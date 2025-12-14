import json

from .gen_filter import *


def tag_lepton_gen(lepton, gen_mask_func):
    
    """
    Generic function to tag leptons with gen-level information
    
    Args:
        lepton: Electron, LowPtElectron, or Muon collection
        gen_mask_func: Function that applies lepton-specific gen mask
    """
    
    signal_mask = (
        (lepton.genPartFlav == 1) 
        & parent_mask(lepton) #Defined for SUSY parents and W/Z bosons
        & gen_mask_func(lepton)
    )
    
    light_fake_mask = (lepton.genPartFlav == 0)
    heavy_decay_mask = ((lepton.genPartFlav == 4) | (lepton.genPartFlav == 5))
    tau_decay_mask = (lepton.genPartFlav == 15)
    
    
    gen_tag = ak.full_like(lepton.pt, -1, dtype=int) #initialize new
    
    # Now reassign -1 to a different binary int
    gen_tag = ak.where(tau_decay_mask, 1000, gen_tag)
    gen_tag = ak.where(heavy_decay_mask, 0100, gen_tag)
    gen_tag = ak.where(light_fake_mask, 0010, gen_tag)
    gen_tag = ak.where(signal_mask, 0001, gen_tag)

    # Any lepton that does not satisfy the above will remain -1, can use that to keep track of how many gen particles we don't explore
    
    lepton = ak.with_field(lepton, gen_tag, "gen_tag")
    
    return lepton

def tag_gens(events):
    
    """Add gen_tag field to all lepton collections"""
    
    events = ak.with_field(
        events, 
        tag_lepton_gen(events.Electron, ele_gen_mask), 
        "Electron"
    )
    
    events = ak.with_field(
        events, 
        tag_lepton_gen(events.LowPtElectron, lpte_gen_mask), 
        "LowPtElectron"
    )
    
    events = ak.with_field(
        events, 
        tag_lepton_gen(events.Muon, muon_gen_mask), 
        "Muon"
    )
    
    return events
    