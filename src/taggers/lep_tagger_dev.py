# Define our skims or (categories) for Electrons, Muons, LowPtElectrons

from .vid_unpacked import *
from .gen_tagger import *
import numpy as np
import awkward as ak

####################################################################
# functions that take the lepton collections, checks if the lepton is baseline, gold, etc., adds a boolean to it if so. Thats it.


def tag_qual(events):
    
    """
    Returns events with new lepton quality fields and combined Leptons collection
    """
    
    # Tag lepton collections
    tagged_electrons = tag_ele_quality(events.Electron)
    tagged_low_pt_electrons = tag_lpte_quality(events.LowPtElectron)
    tagged_muons = tag_muon_quality(events.Muon)
    
    # And now update the events to have the new fields
    events = ak.with_field(events, tagged_electrons, "Electron")
    events = ak.with_field(events, tagged_low_pt_electrons, "LowPtElectron")
    events = ak.with_field(events, tagged_muons, "Muon")
    
    # Concatenate into combined Leptons collection
    leptons = ak.concatenate([
        tagged_electrons,
        tagged_low_pt_electrons,
        tagged_muons
    ], axis=1)
    
    # Add PtEtaPhiMCandidate behavior (for 4-vector operations like .mass and .delta_r)
    leptons = ak.with_name(leptons, "PtEtaPhiMCandidate")
    
    # Add Leptons as a new event field
    events = ak.with_field(events, leptons, "Leptons")
    
    return events


def tag_qual_and_gen(events):
    
    """
    Returns events with new lepton quality fields and gen fields and a combined Leptons collection
    """

    events = tag_gens(events) # Tag the gen-level information, this only works for MC obviously
    
    # Tag lepton collections
    tagged_electrons = tag_ele_quality(events.Electron)
    tagged_low_pt_electrons = tag_lpte_quality(events.LowPtElectron)
    tagged_muons = tag_muon_quality(events.Muon)
    
    # And now update the events to have the new fields
    events = ak.with_field(events, tagged_electrons, "Electron")
    events = ak.with_field(events, tagged_low_pt_electrons, "LowPtElectron")
    events = ak.with_field(events, tagged_muons, "Muon")
    
    # Concatenate into combined Leptons collection
    leptons = ak.concatenate([
        tagged_electrons,
        tagged_low_pt_electrons,
        tagged_muons
    ], axis=1)
    
    # Add PtEtaPhiMCandidate behavior (for 4-vector operations like .mass and .delta_r)
    leptons = ak.with_name(leptons, "PtEtaPhiMCandidate")
    
    # Add Leptons as a new event field
    events = ak.with_field(events, leptons, "Leptons")
    
    return events

    
##################################################
# global values for convenient "switch flipping" #
##################################################

global_sip3d = 3

##################################################
##################################################

def tag_ele_quality(ele):  # use on raw Electron collection (Awkward array)
    
    """
    Add 'isGold', 'isSilver', 'isBronze' boolean field to events.Electron based on cuts.
    Add 'qual_tag' integer field to events.Electron, binary numbers
    """

    # define variables
    abs_eta   = np.abs(ele.eta)
    abs_dxy   = np.abs(ele.dxy)
    abs_dz    = np.abs(ele.dz)
    pt        = ele.pt
    iso03pt   = ele.pfRelIso03_all * pt
    miniIsoPt = ele.miniPFRelIso_all * pt

    # --- Baseline selection ---
    baseline_mask = (
        (pt >= 7)
        & ( ((pt >= 10) & (abs_eta < 2.5)) | ((pt < 10) & (abs_eta < 1.442)) )
        & (ele.sip3d < 6)
        & (abs_dxy < 0.05)
        & (abs_dz  < 0.1)
        & (iso03pt   < (20 + 300/pt))
        & (miniIsoPt < (20 + 300/pt))
        & (ele.lostHits == 0)
        & (ele.convVeto == 1) #added this for regular electrons, see how it works
        & loose_minus_iso_hoe(ele)
    )

    # --- Quality tags ---

    low_pt_pass = (
        (pt < 20)
        & (iso03pt <= 4)
        & (miniIsoPt <= 4)
        & tight_minus_iso_hoe(ele)
    )
    
    high_pt_pass = (
        (pt >= 20)
        & ele.mvaIso_WP90
    )


    gold_silver_mask = baseline_mask & (low_pt_pass | high_pt_pass)
    
    
    gold_mask = (
        gold_silver_mask
        & (ele.sip3d < global_sip3d)
    )

    silver_mask = (
        gold_silver_mask
        & (ele.sip3d >= global_sip3d)
    )

    bronze_mask = baseline_mask & ~gold_silver_mask

    
    # Add quality masks
    ele = ak.with_field(ele, gold_mask, 'isGold')
    ele = ak.with_field(ele, silver_mask, 'isSilver')
    ele = ak.with_field(ele, bronze_mask, 'isBronze')

    # Used to have an 'isBaseline' but its redundant since you can add isGold, isSilver, isBronze to get baseline.
    # Assuming my categories are correctly mutually exclusive
    
    # Create qual_tag int for convenient histogram filling
    qual_tag = ak.full_like(ele.pt, -1, dtype=int)
    qual_tag = ak.where(ele.isBronze, 001, qual_tag)
    qual_tag = ak.where(ele.isSilver, 010, qual_tag)
    qual_tag = ak.where(ele.isGold, 100, qual_tag)
    
    ele = ak.with_field(ele, qual_tag, "qual_tag")

    return ele




####################################################################
# LowPtElectrons


def tag_lpte_quality(lpte): #use on raw lpte collection

    """
    Add 'isGold', 'isSilver', 'isBronze' boolean field to events.LowPtElectron based on cuts.
    Add 'qual_tag' integer field to events.LowPtElectron, binary numbers
    """
    
    # define variables
    abs_eta   = np.abs(lpte.eta)
    sip3d     = _lpte_sip3d(lpte)
    abs_dxy   = np.abs(lpte.dxy)
    abs_dz    = np.abs(lpte.dz)
    central_eta_ID = (
        ((abs_eta >= 0.8) & (abs_eta < 1.442) & (lpte.ID >= 3)) |
        ((abs_eta < 0.8) & (lpte.ID >= 2.3))
    )

    pt        = lpte.pt
    miniIsoPt = lpte.miniPFRelIso_all * pt

    # --- Baseline selection ---
    baseline_mask = (
        ((pt >= 2) & (pt < 7))
        & (abs_eta < 1.442)
        & (sip3d < 6)
        & (abs_dxy < 0.05)
        & (abs_dz  < 0.1)
        & (miniIsoPt < (20 + 300/pt))
        & (lpte.convVeto == 1)
        & (lpte.lostHits == 0)
        & (lpte.ID >= 1.5)
    )

    # --- Quality tags ---

    gold_silver_mask = (
        baseline_mask
        & (miniIsoPt <= 4)
        & central_eta_ID
    )
    
    gold_mask = (
        gold_silver_mask
        & (sip3d < global_sip3d)
    )

    silver_mask = (
        gold_silver_mask
        & (sip3d >= global_sip3d)
    )

    bronze_mask = baseline_mask & ~gold_silver_mask

    lpte = ak.with_field(lpte, gold_mask, 'isGold')
    lpte = ak.with_field(lpte, silver_mask, 'isSilver')
    lpte = ak.with_field(lpte, bronze_mask, 'isBronze')
    
    qual_tag = ak.full_like(lpte.pt, -1, dtype=int)
    qual_tag = ak.where(lpte.isBronze, 001, qual_tag)
    qual_tag = ak.where(lpte.isSilver, 010, qual_tag)
    qual_tag = ak.where(lpte.isGold, 100, qual_tag)
    
    lpte = ak.with_field(lpte, qual_tag, "qual_tag")


    return lpte


####################################################################
# Muons:

def tag_muon_quality(muon): #use on raw muon collection

    """
    Add 'isGold', 'isSilver', 'isBronze' boolean field to events.Muon based on cuts.
    Add 'qual_tag' integer field to events.Muon, binary numbers
    """

    
    # define variables
    abs_eta   = np.abs(muon.eta)
    abs_dxy   = np.abs(muon.dxy)
    abs_dz    = np.abs(muon.dz)
    sip3d     = muon.sip3d
    pt        = muon.pt
    iso03pt   = muon.pfRelIso03_all * pt
    miniIsoPt = muon.miniPFRelIso_all * pt
    pfIsoId   = muon.pfIsoId
    tight     = muon.tightId
    low_pt_not_in_endcap = (pt >= 6) | (abs_eta < 1.2) # if muon has pt < 6, abs_eta < 1.2 allows them to stay in selection 
                                                       # (i.e. removes low pt endcap muons)
    
    
    # --- Baseline selection ---
    baseline_mask = (
        (abs_eta < 2.5)
        & (sip3d < 6)
        & (abs_dxy < 0.05)
        & (abs_dz  < 0.1)
        & (iso03pt   < (20 + 300/pt))
        & (miniIsoPt < (20 + 300/pt))
        & low_pt_not_in_endcap
    )

    # --- Quality tags ---

    gold_silver_mask = ( 
        baseline_mask
        & (iso03pt   <= 4)
        & (miniIsoPt <= 4)
        & tight
    )
    
    gold_mask = (
        gold_silver_mask
        & (sip3d < global_sip3d)
    )
        
    silver_mask = (
        gold_silver_mask
        & (sip3d >= global_sip3d)
    )

    bronze_mask = baseline_mask & ~gold_silver_mask

    muon = ak.with_field(muon, gold_mask, 'isGold')
    muon = ak.with_field(muon, silver_mask, 'isSilver')
    muon = ak.with_field(muon, bronze_mask, 'isBronze')
    
    qual_tag = ak.full_like(muon.pt, -1, dtype=int)
    qual_tag = ak.where(muon.isBronze, 001, qual_tag)
    qual_tag = ak.where(muon.isSilver, 010, qual_tag)
    qual_tag = ak.where(muon.isGold, 100, qual_tag)
    
    muon = ak.with_field(muon, qual_tag, "qual_tag")

    return muon
    


def _lpte_sip3d(lpte):
     
    #approximation or rough calculation based on what Suyash did years ago
    #NOT SIP3D by any means, this is a "SIP3D-like" variable we constructed
    
    dxy = lpte.dxy
    dz = lpte.dz
    dxy_err = lpte.dxyErr
    dz_err = lpte.dzErr
    
    sigma_xy = dxy/dxy_err
    sigma_z = dz/dz_err
    
    SIP3D = np.sqrt(sigma_xy**2 + sigma_z**2)
    
    return SIP3D