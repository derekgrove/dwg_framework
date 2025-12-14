import numpy
import awkward as ak
import sys
import hist.dask as dah
from pathlib import Path

current_dir = Path.cwd()

# taggers = current_dir.parent / 'taggers'

from .plotting.histers import *
#from ../taggers.lep_tagger import *
#from .taggers.lep_tagger import *
#from .taggers.gen_tagger import *

class default_analysis:

    def __init__(self, name, hist=None):
        
        self.name = name
        self.hist = hist or (
        dah.Hist.new
        .Regular(500, 0, 500, name="pt") #500 (1 GeV) bins between 0 and 500
        .Double()
        )
        #self.pt = pt
        #self.pt = getattr(self, 'pt')

"""

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
"""




def poster_plots_dict(obj):
    
    """
    structure will be a dict with lots of hists of various configurations
    """
    
    ##############
    # fill hists #
    #############

    gens = [-10, 10, 11, 12, 13] #I added 1 in front for gens so I know I don't accidentally get it mixed up with qual
    quals = [-1,1,2,3]
    
    pt_eta_gen_qual = make_2d2d_hist_cat(
        obj,                           
        [2,3,4,5,7,10,20,45,75,1000], 
        [0,0.8,1.4442,1.556,2.5],
        cat1_binning = gens,
        cat2_binning = quals,
        var1_name = 'pt',
        var2_name = 'eta',
        cat1_name='gen_tag',
        cat2_name='qual_tag',
        var2_abs=True
       )

    pt_gen_qual = make_1d2d_hist_reg_cat(
    obj,
    [80,0,20], #like [100, 0, 100] 100 bins between 0 and 100
    cat1_binning = gens,
    cat2_binning = quals,
    var_name="pt",
    cat1_name="gen_tag",
    cat2_name="qual_tag",
    var_abs=False,
    )

    ################
    # fill results #
    ###############

    results = {
        "pt_eta_gen_qual_hist": pt_eta_gen_qual,
        "pt_gen_qual_hist": pt_gen_qual,
    }

    return results




