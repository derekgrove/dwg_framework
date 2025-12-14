# Define our skims or (categories) for Electrons, Muons, LowPtElectrons
import json

from vid_unpacked import *
from gen_tagger import *
import numpy as np
import awkward as ak

####################################################################
# functions that take the lepton collections, checks if the lepton is baseline, gold, etc., adds a boolean to it if so. Thats it.

def tag_qual(obj, ID):
    
    acceptable_IDs = [
        'ele', 'electron',
        'lpte', 'lowptelectron',
        'mu', 'muon']

    ID_lower = ID.lower()

    if ID_lower not in acceptable_IDs:
        sys.exit(f'ID {ID} not acceptable, must be in {acceptable_IDs}')
    
    if ID_lower in ['ele', 'electron']:
        obj = tag_ele_quality(obj)
    elif ID_lower in ['lpte', 'lowptelectron']: 
        obj = tag_lpte_quality(obj)
    elif ID_lower in ['mu', 'muon']: 
        obj = tag_muon_quality(obj)

    return obj

def tag_ele(ele):
    
    return tag_qual(tag_gen(ele, 'ele'), 'ele')

def tag_lpte(lpte):
    
    return tag_qual(tag_gen(lpte, 'lpte'), 'lpte')

def tag_muon(muon):
    
    return tag_qual(tag_gen(muon, 'muon'), 'muon')

    
def tag_and_combine_ele(electron, lowptelectron):
    
    tagged_ele = tag_ele(electron)
    tagged_lpte = tag_lpte(lowptelectron)

    pt_selected_electron = tagged_ele[tagged_ele.pt >= 7]
    pt_selected_lpte = tagged_lpte[tagged_lpte.pt < 7]

    ele = ak.concatenate([pt_selected_electron, pt_selected_lpte], axis=1)

    return ele


##################################################
# global values for convenient "switch flipping" #
# #################################################

global_sip3d = 3

##################################################
# #################################################

def tag_ele_quality(ele):  # use on raw Electron collection (Awkward/NanoEvents)
    
    """
    Add 'baseline', 'gold', 'silver', 'bronze' boolean fields to each electron.
    (and some testing fields for vid_unpacked)
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
        & ele.mvaFall17V2Iso_WP90 #trying this for now, unsure which one performs best for UL
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

    ele['isBaseline'] = baseline_mask
    ele['isGold']     = gold_mask
    ele['isSilver']   = silver_mask
    ele['isBronze']   = bronze_mask

    ele["qual_tag"] = -1 # Filling everything with dummy values for now
    ele["qual_tag"] = ak.where(ele.isBaseline, 0, ele.qual_tag) #this just always got overwritten by the bronze, silver, or gold
    ele["qual_tag"] = ak.where(ele.isBronze, 1, ele.qual_tag) # similar for 1
    ele["qual_tag"] = ak.where(ele.isSilver, 10, ele.qual_tag) # similar for 1
    ele["qual_tag"] = ak.where(ele.isGold, 100, ele.qual_tag) # similar for 1

    return ele




####################################################################
# LowPtElectrons


def tag_lpte_quality(lpte): #use on raw lpte collection

    """
    Add an 'isBaseline', 'isGold', 'isSilver', 'isBronze' field to each LowPtElectron based on cuts.
    """
    
    # define variables
    pt        = lpte.pt
    abs_eta   = np.abs(lpte.eta)
    sip3d     = lpte_sip3d(lpte)
    abs_dxy   = np.abs(lpte.dxy)
    abs_dz    = np.abs(lpte.dz)
    miniIsoPt = lpte.miniPFRelIso_all * pt
    #central_eta_ID = (
    #    ((abs_eta >= 0.8) & (abs_eta < 1.442) & (lpte.ID >= 3)) |
    #    ((abs_eta < 0.8) & (lpte.ID >= 2.3))
        #((abs_eta < 0.8) & (lpte.ID >= 2.5))
    #)
    
    central_eta_ID = (
        ((pt < 4) &
        (((abs_eta >= 0.8) & (abs_eta < 1.442) & (lpte.ID >= 3)) |
        ((abs_eta < 0.8) & (lpte.ID >= 2.6)))) |
        ((pt >= 4) &
        (((abs_eta >= 0.8) & (abs_eta < 1.442) & (lpte.ID >= 3.2)) |
        ((abs_eta < 0.8) & (lpte.ID >= 2.8))))
    )
   

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
        & (lpte.ID >= 2)
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


    lpte['isBaseline'] = baseline_mask
    lpte['isGold']     = gold_mask
    lpte['isSilver']   = silver_mask
    lpte['isBronze']   = bronze_mask

    lpte["qual_tag"] = -1 # Filling everything with dummy values for now
    lpte["qual_tag"] = ak.where(lpte.isBaseline, 0, lpte.qual_tag)
    lpte["qual_tag"] = ak.where(lpte.isBronze, 1, lpte.qual_tag)
    lpte["qual_tag"] = ak.where(lpte.isSilver, 10, lpte.qual_tag)
    lpte["qual_tag"] = ak.where(lpte.isGold, 100, lpte.qual_tag)
    
    return lpte


####################################################################
# Muons:

def tag_muon_quality(muon): #use on raw muon collection

    """
    Add 'baseline', 'gold', 'silver', 'bronze' boolean fields to each electron.
    (and some testing fields for vid_unpacked)
    """

    
    # define variables
    abs_eta   = np.abs(muon.eta)
    abs_dxy   = np.abs(muon.dxy)
    abs_dz    = np.abs(muon.dz)
    sip3d = muon.sip3d
    pt        = muon.pt
    iso03pt   = muon.pfRelIso03_all * pt
    miniIsoPt = muon.miniPFRelIso_all * pt
    pfIsoId = muon.pfIsoId
    tight = muon.tightId
    low_pt_not_in_endcap = (muon.pt < 6) & (abs_eta <= 1.2)
    
    
    # --- Baseline selection ---
    baseline_mask = (
        (abs_eta < 2.5)
        & low_pt_not_in_endcap
        & (sip3d < 6)
        & (abs_dxy < 0.05)
        & (abs_dz  < 0.1)
        & (iso03pt   < (20 + 300/pt))
        & (miniIsoPt < (20 + 300/pt))
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

    muon['isBaseline'] = baseline_mask
    muon['isGold']     = gold_mask
    muon['isSilver']   = silver_mask
    muon['isBronze']   = bronze_mask

    muon["qual_tag"] = -1 # Filling everything with dummy values for now
    muon["qual_tag"] = ak.where(muon.isBaseline, 0, muon.qual_tag)
    muon["qual_tag"] = ak.where(muon.isBronze, 1, muon.qual_tag)
    muon["qual_tag"] = ak.where(muon.isSilver, 10, muon.qual_tag)
    muon["qual_tag"] = ak.where(muon.isGold, 100, muon.qual_tag)

    return muon
    


def lpte_sip3d(lpte):
     
    #approximation or rough calculation based on what Suyash did years ago
    
    dxy = lpte.dxy
    dz = lpte.dz
    dxy_err = lpte.dxyErr
    dz_err = lpte.dzErr
    
    sigma_xy = dxy/dxy_err
    sigma_z = dz/dz_err
    
    SIP3D = np.sqrt(sigma_xy**2 + sigma_z**2)
    
    return SIP3D
