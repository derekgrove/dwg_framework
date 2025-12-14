# Define our skims or (categories) for Electrons, Muons, LowPtElectrons

from .vid_unpacked import *
from .gen_tagger import *
import numpy as np
import awkward as ak

####################################################################
# functions that take the lepton collections, checks if the lepton is baseline, gold, etc., adds a boolean to it if so. Thats it.


def tag_qual_no_gen(events):
    
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

    #ele['isBaseline'] = baseline_mask # redundant since isGold, isSilver and isBronze sum to isBaseline, useful for debugging though
    ele['isGold']     = gold_mask
    ele['isSilver']   = silver_mask
    ele['isBronze']   = bronze_mask

    ele["qual_tag"] = -1 # Filling everything with dummy values for now
    #ele["qual_tag"] = ak.where(ele.isBaseline, 0, ele.qual_tag) # also redundant but useful for debugging by making sure gsb are mutually exlusive (no electrons should be left with "qual_tag == 0"
    ele["qual_tag"] = ak.where(ele.isBronze, 1, ele.qual_tag) # similar for 1
    ele["qual_tag"] = ak.where(ele.isSilver, 10, ele.qual_tag) # similar for 1
    ele["qual_tag"] = ak.where(ele.isGold, 100, ele.qual_tag) # similar for 1

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

    #lpte['isBaseline'] = baseline_mask
    lpte['isGold']     = gold_mask
    lpte['isSilver']   = silver_mask
    lpte['isBronze']   = bronze_mask

    lpte["qual_tag"] = -1
    #lpte["qual_tag"] = ak.where(lpte.isBaseline, 0, lpte.qual_tag) #useful for debugging, see tag_electron_quality 
    lpte["qual_tag"] = ak.where(lpte.isBronze, 1, lpte.qual_tag)
    lpte["qual_tag"] = ak.where(lpte.isSilver, 10, lpte.qual_tag)
    lpte["qual_tag"] = ak.where(lpte.isGold, 100, lpte.qual_tag)


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

    #muon['isBaseline'] = baseline_mask
    muon['isGold']     = gold_mask
    muon['isSilver']   = silver_mask
    muon['isBronze']   = bronze_mask

    muon["qual_tag"] = -1
    #muon["qual_tag"] = ak.where(muon.isBaseline, 0, muon.qual_tag) #useful for debugging, see tag_electron_quality 
    muon["qual_tag"] = ak.where(muon.isBronze, 1, muon.qual_tag)
    muon["qual_tag"] = ak.where(muon.isSilver, 10, muon.qual_tag)
    muon["qual_tag"] = ak.where(muon.isGold, 100, muon.qual_tag)

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