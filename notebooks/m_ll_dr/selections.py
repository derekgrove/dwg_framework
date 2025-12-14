
# Now do the custom imports
from coffea.analysis_tools import PackedSelection
import lep_tagger as tagger_preUL
import lep_tagger_UL as tagger_UL
import awkward as ak
import numpy as np

import hist

def concat_leps(my_events, is_UL):

    tagger = tagger_UL if is_UL else tagger_preUL
    tag_qual = tagger.tag_qual
    
    ele = tag_qual(my_events.Electron, 'ele')
    mu = tag_qual(my_events.Muon, 'mu')
    return ak.with_name(ak.concatenate([ele, mu], axis=1), "PtEtaPhiMCandidate")

import hist

def make_dr_hist(my_events, is_UL):
    leptons = concat_leps(my_events, is_UL)
    
    # Select only events with exactly 2 gold leptons
    dilepton_mask = ak.num(leptons) == 2
    leptons = leptons[dilepton_mask]
    
    # Get the two leptons
    lead = leptons[:, 0]
    trail = leptons[:, 1]
    
    # Compute ΔR using the delta_r method
    dr = lead.delta_r(trail)
    
    # Calculate invariant mass
    mass = (lead + trail).mass

    qual_sum = ak.sum(leptons.qual_tag, axis=1)

    # Check unique values
    #unique_qual_tag = np.unique(ak.flatten(leptons.qual_tag))
    #unique_qual_sum = np.unique(qual_sum)

    
    # Use regular hist.Hist (works for both dask and virtual)
    h_dr = hist.Hist(
        hist.axis.Regular(160, 0, 8, name="dr", label=r"$\Delta R$"),
        hist.axis.Regular(300, 0, 30, name="mass", label=r"$m_{\ell\ell}$ [GeV]"),
        hist.axis.IntCategory([-2, 0, 2, 9, 11, 20, 99, 101, 110, 200], name="qual_sum")   
    )
    h_dr.fill(dr=dr, mass=mass, qual_sum=qual_sum)

    
    
    return h_dr


def make_dr_div_mll_hist(my_events, is_UL):
    leptons = concat_leps(my_events, is_UL)
    
    # Select only events with exactly 2 gold leptons
    dilepton_mask = ak.num(leptons) == 2
    leptons = leptons[dilepton_mask]
    
    # Get the two leptons
    lead = leptons[:, 0]
    trail = leptons[:, 1]
    
    # Compute ΔR using the delta_r method
    dr = lead.delta_r(trail)
    
    # Calculate invariant mass
    mass = (lead + trail).mass

    dr_div_m = dr/mass

    qual_sum = ak.sum(leptons.qual_tag, axis=1)


    
    # Use regular hist.Hist (works for both dask and virtual)
    h = hist.Hist(
        hist.axis.Regular(150, 0, 1.5, name="dr_div_m", label=r"$\frac{\Delta R}{m_{\ell\ell}}$"),
        hist.axis.IntCategory([-2, 0, 2, 9, 11, 20, 99, 101, 110, 200], name="qual_sum")
        
    )
    h.fill(dr_div_m=dr_div_m, qual_sum=qual_sum)

    
    
    return h

def get_selection(events, is_UL):
    # Create a NEW selection for each call
    
    selection = PackedSelection()
    
    leptons = concat_leps(events, is_UL)

    dilepton_mask = ak.num(leptons) == 2
    
    padded_leptons = ak.pad_none(leptons, 2, axis=1)
    lead = padded_leptons[:, 0]
    trail = padded_leptons[:, 1]

    # Calculate invariant mass (will be None where lead or trail is None)
    mass = (lead + trail).mass
    
    # Compute ΔR (will be None where lead or trail is None)
    dr = lead.delta_r(trail)
   
    # Flavor selections
    selection.add_multiple(
        {
            # Dilepton events (any combination of 2 leptons)
            "dilep": dilepton_mask,
            
            # Same flavor
            "mumu": (ak.num(events.Muon) == 2) & (ak.num(events.Electron) == 0),
            "ee": (ak.num(events.Electron) == 2) & (ak.num(events.Muon) == 0),
            
            # Opposite flavor
            "emu": (ak.num(events.Electron) == 1) & (ak.num(events.Muon) == 1),
            
            # Quality selections - both leptons are baseline
            "bl_leps": ak.all(leptons.qual_tag != -1, axis=1),  # Fixed: added axis=1
            "gold_leps": ak.all(leptons.qual_tag == 3, axis=1),  # Fixed: added axis=1
    
            "ss": (ak.sum(leptons.charge, axis=1) != 0),  # must be used with dilep
            "os": (ak.sum(leptons.charge, axis=1) == 0),  # must be used with dilep
        }
    )

    return selection

def get_event_combos(events, is_UL):
    selection = get_selection(events, is_UL)

    return {
        #"events": events,
        "dilep": events[selection.all("dilep")],
        
        #"mumu": events[selection.all("mumu")],
        #"ee": events[selection.all("ee")],
        
        #"dilep_ss": events[selection.all("dilep", "ss")],
        #"dilep_os": events[selection.all("dilep", "os")],
        #"dilep_sf": events[selection.any("ee", "mumu")],
        #"dilep_of": events[selection.all("emu")],
        
        "dilep_ss_of": events[selection.all("dilep", "ss", "emu")],
        "dilep_os_of": events[selection.all("dilep", "os", "emu")],
        "dilep_ss_sf": events[selection.all("dilep", "ss") & selection.any("ee", "mumu")],
        "dilep_os_sf": events[selection.all("dilep", "os") & selection.any("ee", "mumu")],
        
        "dilep_ss_ee": events[selection.all("dilep", "ss", "ee")],
        "dilep_os_ee": events[selection.all("dilep", "os", "ee")],
        
        "dilep_ss_mumu": events[selection.all("dilep", "ss", "mumu")],
        "dilep_os_mumu": events[selection.all("dilep", "os", "mumu")],
    }
