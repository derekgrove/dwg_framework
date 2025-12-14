import json

from analysis_tools.taggers.gen_filter import *


def tag_gen(obj, obj_name): 

    """Adds a new field to objs (Electron, LowPtElectron, or Muon collections) named 'gen_tag', can access it like Electron.gen_tag (after you run this function and reassign the electron collection to this one thats returned). 

    The resulting integers correspond like this:

    obj.gen_tag == 0, SIGNAL
    obj.gen_tag == 1, UNMATCHED TO PV
    obj.gen_tag == 2, DECAY FROM B OR C
    obj.gen_tag == 3, DECAY FROM TAU
    obj.gen_tag == 999, OTHER DECAY (any other genflav than the above)
    
    """

    if obj_name.lower() in ["ele", "electron"]:

        signal_mask = (
            (obj.genPartFlav == 1) 
            & parent_mask(obj)
            & ele_gen_mask(obj)
        )
        
    elif obj_name.lower() in ["lpte", "lowptelectron"]:
        
        signal_mask = (
            (obj.genPartFlav == 1) 
            & parent_mask(obj)
            & lpte_gen_mask(obj)
        )
        
    elif obj_name.lower() in ["mu", "muon"]:
        
        signal_mask = (
            (obj.genPartFlav == 1) 
            & parent_mask(obj)
            & muon_gen_mask(obj)
        )
        
    else:
        sys.exit(f"invalid obj_name: {obj_name}")
        
    light_fake_mask = (obj.genPartFlav == 0)
    heavy_decay_mask = ( (obj.genPartFlav == 4) | (obj.genPartFlav == 5) )
    tau_decay_mask = (obj.genPartFlav == 15)
    

    
    obj["gen_tag"] = -10 # Filling everything with dummy values for now
    obj["gen_tag"] = ak.where(signal_mask, 10, obj.gen_tag)
    obj["gen_tag"] = ak.where(light_fake_mask, 11, obj.gen_tag)
    obj["gen_tag"] = ak.where(heavy_decay_mask, 12, obj.gen_tag)
    obj["gen_tag"] = ak.where(tau_decay_mask, 13, obj.gen_tag)

#    obj["gen_tag"] = 'other_gen' # Filling everything with dummy values for now
#    obj["gen_tag"] = ak.where(signal_mask, 'signal', obj.gen_tag)
#    obj["gen_tag"] = ak.where(light_fake_mask, 'light_fake', obj.gen_tag)
#    obj["gen_tag"] = ak.where(heavy_decay_mask, 'heavy_decay', obj.gen_tag)
#    obj["gen_tag"] = ak.where(tau_decay_mask, 'tau_decay', obj.gen_tag)
    
    return obj



