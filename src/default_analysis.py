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




def lep_analysis_dict(obj):
    
    """
    structure will be a dict with lots of hists of various configurations
    """
    
    ##############
    # fill hists #
    #############

    gens = [-10, 10, 11, 12, 13] #I added 1 in front for gens so I know I don't accidentally get it mixed up with qual
    quals = [-1,1,2,3]
    
    pt_eta_hist = make_2d2d_hist_cat(
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
    
    pt_eta_hist_muon_v1 = make_2d2d_hist_cat(
        obj,                           
        [2,3,4,5,7,10,20,45,75,1000], 
        [0,0.8,1.556,2.5],
        cat1_binning = gens,
        cat2_binning = quals,
        var1_name = 'pt',
        var2_name = 'eta',
        cat1_name='gen_tag',
        cat2_name='qual_tag',
        var2_abs=True
       )

    pt_eta_hist_muon_v2 = make_2d2d_hist_cat(
        obj,                           
        [2,3,4,5,7,10,20,45,75,1000], 
        [0,0.9,1.2,2.1,2.4],
        cat1_binning = gens,
        cat2_binning = quals,
        var1_name = 'pt',
        var2_name = 'eta',
        cat1_name='gen_tag',
        cat2_name='qual_tag',
        var2_abs=True
       )

    pt_AN_v1 = make_1d2d_hist_reg_cat(
    obj,
    [40,0,100], #like [100, 0, 100] 100 bins between 0 and 100
    cat1_binning = gens,
    cat2_binning = quals,
    var_name="pt",
    cat1_name="gen_tag",
    cat2_name="qual_tag",
    var_abs=False,
    )

    pt_bins_v2 = [2,3,4,5,7,10,12.5,15.0,17.5,20.0,
               22.5,25.0,27.5,30.0,32.5,35.0,
               37.5,40.0,42.5,45.0,47.5,50.0,
               52.5,55.0,57.5,60.0,62.5,65.0,
               67.5,70.0,72.5,75.0,77.5,80.0,
               82.5,85.0,87.5,90.0,92.5,95.0,
               97.5,100.0]
    
    pt_AN_v2 = make_1d2d_hist_var_cat(
    obj,
    pt_bins_v2, 
    cat1_binning = gens,
    cat2_binning = quals,
    var_name="pt",
    cat1_name="gen_tag",
    cat2_name="qual_tag",
    var_abs=False,
    )

    pt_bins_v3 = [2,5,7.5,10,12.5,15.0,17.5,20.0,
               22.5,25.0,27.5,30.0,32.5,35.0,
               37.5,40.0,42.5,45.0,47.5,50.0,
               52.5,55.0,57.5,60.0,62.5,65.0,
               67.5,70.0,72.5,75.0,77.5,80.0,
               82.5,85.0,87.5,90.0,92.5,95.0,
               97.5,100.0]
    
    pt_AN_v3 = make_1d2d_hist_var_cat(
    obj,
    pt_bins_v3, 
    cat1_binning = gens,
    cat2_binning = quals,
    var_name="pt",
    cat1_name="gen_tag",
    cat2_name="qual_tag",
    var_abs=False,
    )

    pt_bins_muon = [3,5,7.5,10,12.5,15.0,17.5,20.0,
               22.5,25.0,27.5,30.0,32.5,35.0,
               37.5,40.0,42.5,45.0,47.5,50.0,
               52.5,55.0,57.5,60.0,62.5,65.0,
               67.5,70.0,72.5,75.0,77.5,80.0,
               82.5,85.0,87.5,90.0,92.5,95.0,
               97.5,100.0]
    
    pt_AN_muon = make_1d2d_hist_var_cat(
    obj,
    pt_bins_muon, 
    cat1_binning = gens,
    cat2_binning = quals,
    var_name="pt",
    cat1_name="gen_tag",
    cat2_name="qual_tag",
    var_abs=False,
    )

    pt_bins_v4 = [1,2,3,4,5,6,7,8,9,10,
                  12.5,15.0,17.5,20.0,
               22.5,25.0,27.5,30.0,32.5,35.0,
               37.5,40.0,42.5,45.0,47.5,50.0,
               52.5,55.0,57.5,60.0,62.5,65.0,
               67.5,70.0,72.5,75.0,77.5,80.0,
               82.5,85.0,87.5,90.0,92.5,95.0,
               97.5,100.0]
    
    pt_AN_v4 = make_1d2d_hist_var_cat(
    obj,
    pt_bins_v4, 
    cat1_binning = gens,
    cat2_binning = quals,
    var_name="pt",
    cat1_name="gen_tag",
    cat2_name="qual_tag",
    var_abs=False,
    )

    

# (obj, var_binning, cat_binning, var_name = "pt", cat_name="genPartFlav")

    ################
    # fill results #
    ###############

    results = {
        "pt_eta_hist": pt_eta_hist,
        "pt_eta_hist_muon_v1": pt_eta_hist_muon_v1,
        "pt_eta_hist_muon_v2": pt_eta_hist_muon_v2,
        "pt_AN_hist_v1": pt_AN_v1,
        "pt_AN_hist_v2": pt_AN_v2,
        "pt_AN_hist_v3": pt_AN_v3,
        "pt_AN_hist_muon": pt_AN_muon,
        "pt_AN_hist_v4": pt_AN_v4,
    }

    return results


def lpte_analysis_dict(obj): # has lpte specific variables, do not run on a collection without these variables
    
    """
    structure will be a dict with lots of hists of various configurations
    """
    
    ##############
    # fill hists #
    #############

    gens = [-10, 10, 11, 12, 13] #I added 1 in front for gens so I know I don't accidentally get it mixed up with qual
    quals = [-1,1,2,3]
    
    pt_eta_hist = make_2d2d_hist_cat(
        obj,                           
        [1,2,3,4,5,7,10,20,45,75,1000], 
        [0,0.8,1.4442,1.556,2.5],
        cat1_binning = gens,
        cat2_binning = quals,
        var1_name = 'pt',
        var2_name = 'eta',
        cat1_name='gen_tag',
        cat2_name='qual_tag',
        var2_abs=True
       )

    pt_ID_hist = make_2d2d_hist_cat_2reg(
        obj,                           
        [20,0,20], 
        [100,0,10],
        cat1_binning = gens,
        cat2_binning = quals,
        var1_name = 'pt',
        var2_name = 'ID',
        cat1_name='gen_tag',
        cat2_name='qual_tag',
        var2_abs=True
       )  

    ID_eta_hist = make_2d2d_hist_cat_1reg(
        obj,              
        [0,0.8,1.4442,1.556,2.5],
        [100,0,10], 
        cat1_binning = gens,
        cat2_binning = quals,
        var1_name = 'eta',
        var2_name = 'ID',
        cat1_name='gen_tag',
        cat2_name='qual_tag',
        var1_abs=True
       )  

    ################
    # fill results #
    ###############

    results = {
        "pt_eta_hist": pt_eta_hist,
        "pt_ID_hist": pt_ID_hist,
        "ID_eta_hist": ID_eta_hist,
    }

    return results

