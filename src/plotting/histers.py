#file for defining generic yet precise histogramming functions

import numpy as np
import awkward as ak

import hist.dask as dah

"""
Histograms are written to be multi-dimensional, i.e. they will have many axis that you must then later "project out" the axis you are interested in. For example, a pt hist (1 axis is pt) and then 4 GenFlav axis. It is a "5D" hist, but we project out only "pt and genflav=1" to see the pt distribution of genflav1 whatever.

this is useful because do you want to call a function 3 times to fill the exact same hist and then juggle 3 different hist objects or do you want just one hist object thats like a "folder" of all the hists of the same type you're interested in? I know my answer.
"""


###################
## 1D Histograms ##
###################


def make_1d_pt_hist(obj): #I intend to have the user rebin this later
        
    hist = (
        dah.Hist.new
        .Regular(500, 0, 500, name="pt") #500 (1 GeV) bins between 0 and 500
        .Double()
        .fill(
            pt=ak.flatten(obj.pt),
            )
        )
    
    return hist
    

def make_1d_pt_hist_var(obj, pt_bins):
 
    #pt_bins = [2,3,4,5,7,8,10,15,20,30,45,60,75,500]

    hist = (
        dah.Hist.new
        .Variable(pt_bins, name="pt") 
        .Double()
        .fill(
            pt=np.abs(ak.flatten(obj.eta)),
            )
        )
    
    return hist

    

def make_1d_eta_hist(obj): #also intended to be rebinned later
        
    hist = (
        dah.Hist.new
        .Regular(300, 0, 3, name="eta") #300 (0.01) bins between 0 and 3 |η|
        .Double()
        .fill(
            eta=np.abs(ak.flatten(obj.eta)),
            )
        )
    
    return hist
    

def make_1d_eta_hist_var(obj, eta_bins):

    #eta_bins = [0, 0.8, 1.442, 2.5, 2.8]
    
    hist = (
        dah.Hist.new
        .Variable(eta_bins, name="eta") 
        .Double()
        .fill(
            eta=np.abs(ak.flatten(obj.eta)),
            )
        )
    
    return hist


###################
## 2D Histograms ##
###################

def make_2d_hist(obj, var1_binning, var2_binning, var1_name = "pt", var2_name = "eta"):

    obj_var1 = getattr(obj, var1_name)
    obj_var2 = getattr(obj, var2_name)

    #var1_binning = [ 1, 2, 3, 4, 5,
    #            6, 7, 8, 9, 10,
    #            11, 12, 13, 14, 15,
    #            16, 17, 18, 19, 20]
    
    # can rebin later like:
    # hist[:, :10] → selects the first 10 bins in the second axis
    
    hist = (
        dah.Hist.new
        .Variable(var1_binning, name=var1_name, label=var1_name)
        .Variable(var2_binning, name=var2_name, label=var2_name)
        .Double()
        .fill(
            #eta=np.abs(ak.flatten(obj.eta)),
            **{var1_name: ak.flatten(obj_var1),
               var2_name: ak.flatten(obj_var2)}
            )
        )
    
    return hist


############################
## 1D Category Histograms ##
############################

def make_1d_hist_var_cat(obj, var_binning, cat_binning, var_name = "pt", cat_name="genPartFlav"):

    obj_var = getattr(obj, var_name)
    obj_cat = getattr(obj, cat_name)

    #var1_binning = [ 1, 2, 3, 4, 5,
    #            6, 7, 8, 9, 10,
    #            11, 12, 13, 14, 15,
    #            16, 17, 18, 19, 20]
    
    # can rebin later like:
    # hist[:, :10] → selects the first 10 bins in the second axis
    
    hist = (
        dah.Hist.new
        .Variable(var_binning, name=var_name, label=var_name)
        .IntCat(cat_binning, name=cat_name)
        .Double()
        .fill(
            #eta=np.abs(ak.flatten(obj.eta)),
            **{var_name: ak.flatten(obj_var),
               cat_name:ak.flatten(obj_cat)
            }
            )
        )
    
    return hist
    
"""
These are an extension of the above 1D histograms, they have addtional axis known as 'categories'.

categories are very general, can be whatever. If you want to make a selection of your objects based on some variable like, for example, GenFlav (whether your object was from the PV, unmatched to PV, from a tau, etc.) then you can fill these hists under the same binning but keep the histograms separated by this category. Later, you "project out" the histogram of interest, and you could add them together to get combinations like "heavy decay" = GenFlav15 + GenFlav4 + GenFlav5. Enjoy
"""

############################
## 2D Category Histograms ##
############################



def make_2d_hist_cat(obj, var1_binning, var2_binning, cat_binning, var1_name = "pt", var2_name = "eta", cat_name="genPartFlav"):

    obj_var1 = getattr(obj, var1_name)
    obj_var2 = getattr(obj, var2_name)
    obj_cat = getattr(obj, cat_name)

    #var1_binning = [ 1, 2, 3, 4, 5,
    #            6, 7, 8, 9, 10,
    #            11, 12, 13, 14, 15,
    #            16, 17, 18, 19, 20]
    
    # can rebin later like:
    # hist[:, :10] → selects the first 10 bins in the second axis
    
    hist = (
        dah.Hist.new
        .Variable(var1_binning, name=var1_name, label=var1_name)
        .Variable(var2_binning, name=var2_name, label=var2_name)
        .IntCat(cat_binning, name=cat_name)
        .Double()
        .fill(
            #eta=np.abs(ak.flatten(obj.eta)),
            **{var1_name: ak.flatten(obj_var1),
               var2_name: ak.flatten(obj_var2),
               cat_name:ak.flatten(obj_cat)
            }
            )
        )
    
    return hist


###############################
## 1Dx2D Category Histograms ##
###############################

def make_1d2d_hist_reg_cat(
    obj,
    reg_binning, #like [100, 0, 100] 100 bins between 0 and 100
    cat1_binning,
    cat2_binning,
    var_name="pt",
    cat1_name="gen_tag",
    cat2_name="qual_tag",
    var_abs=False,
    ):
    
    """
    Create a 1D x 2D histogram with one numeric variable (and regular ((interpolated)) binning) and two integer categories (typically gen and quality).
    """

    # Extract variables and categories from the object
    obj_var = getattr(obj, var_name)
    
    obj_cat1 = getattr(obj, cat1_name)
    obj_cat2 = getattr(obj, cat2_name)

    # take absolute value if needed (like for eta or phi, maybe)
    if var_abs:
        obj_var = np.abs(obj_var)

    # flatten arrays, store in dict to be deconstructed later in the fill method
    flat_vars = {
        var_name: ak.flatten(obj_var),
        cat1_name: ak.flatten(obj_cat1),
        cat2_name: ak.flatten(obj_cat2)
    }

    bins = reg_binning[0]
    start = reg_binning[1]
    end = reg_binning[2]

    bins, start, stop = reg_binning

    hist = (
        dah.Hist.new
        .Reg(bins, start, stop, name=var_name, label=var_name)
        .IntCat(cat1_binning, name=cat1_name)
        .IntCat(cat2_binning, name=cat2_name)
        .Double()
    )

    # deconstruct dict, it just works
    hist.fill(**flat_vars) 

    return hist

def make_1d2d_hist_reg_cat_weighted(
    obj,
    reg_binning, #like [100, 0, 100] 100 bins between 0 and 100
    cat1_binning,
    cat2_binning,
    var_name="pt",
    cat1_name="gen_tag",
    cat2_name="qual_tag",
    var_abs=False,
    ):
    
    """
    Create a 1D x 2D histogram with one numeric variable (and regular ((interpolated)) binning) and two integer categories (typically gen and quality).
    """

    # Extract variables and categories from the object
    obj_var = getattr(obj, var_name)
    
    obj_cat1 = getattr(obj, cat1_name)
    obj_cat2 = getattr(obj, cat2_name)
    weight = getattr(obj, "Weight")

    # take absolute value if needed (like for eta or phi, maybe)
    if var_abs:
        obj_var = np.abs(obj_var)

    # flatten arrays, store in dict to be deconstructed later in the fill method
    flat_vars = {
        var_name: ak.flatten(obj_var),
        cat1_name: ak.flatten(obj_cat1),
        cat2_name: ak.flatten(obj_cat2),
        weight:    ak.flatten(weight)
    }

    bins = reg_binning[0]
    start = reg_binning[1]
    end = reg_binning[2]

    bins, start, stop = reg_binning

    hist = (
        dah.Hist.new
        .Reg(bins, start, stop, name=var_name, label=var_name)
        .IntCat(cat1_binning, name=cat1_name)
        .IntCat(cat2_binning, name=cat2_name)
        .Weight(obj_weight)
        .Double()
    )

    # deconstruct dict, it just works
    hist.fill(**flat_vars) 

    return hist

    
def make_1d2d_hist_var_cat(
    obj,
    var_binning,
    cat1_binning,
    cat2_binning,
    var_name="pt",
    cat1_name="gen_tag",
    cat2_name="qual_tag",
    var_abs=False,
    ):
    
    """
    Create a 1D x 2D histogram with one numeric variable and two integer categories (typically gen and quality).
    """

    # Extract variables and categories from the object
    obj_var = getattr(obj, var_name)
    
    obj_cat1 = getattr(obj, cat1_name)
    obj_cat2 = getattr(obj, cat2_name)

    # take absolute value if needed (like for eta or phi, maybe)
    if var_abs:
        obj_var = np.abs(obj_var)

    # flatten arrays, store in dict to be deconstructed later in the fill method
    flat_vars = {
        var_name: ak.flatten(obj_var),
        cat1_name: ak.flatten(obj_cat1),
        cat2_name: ak.flatten(obj_cat2)
    }

    hist = (
        dah.Hist.new
        .Variable(var_binning, name=var_name, label=var_name)
        .IntCat(cat1_binning, name=cat1_name)
        .IntCat(cat2_binning, name=cat2_name)
        .Double()
    )

    # deconstruct dict, it just works
    hist.fill(**flat_vars) 

    return hist

###############################
## 2Dx2D Category Histograms ##
###############################

def make_2d2d_hist_cat(
    obj,
    var1_binning,
    var2_binning,
    cat1_binning,
    cat2_binning,
    var1_name="pt",
    var2_name="eta",
    cat1_name="gen_tag",
    cat2_name="qual_tag",
    var1_abs=False,
    var2_abs=False
    ):
    
    """
    Create a 2D x 2D histogram with two numeric variables and two integer categories (typically gen and quality).
    """

    # Extract variables and categories from the object
    obj_var1 = getattr(obj, var1_name)
    obj_var2 = getattr(obj, var2_name)
    obj_cat1 = getattr(obj, cat1_name)
    obj_cat2 = getattr(obj, cat2_name)

    # take absolute value if needed (like for eta or phi, maybe)
    if var1_abs:
        obj_var1 = np.abs(obj_var1)
    if var2_abs:
        obj_var2 = np.abs(obj_var2)

    # flatten arrays, store in dict to be deconstructed later in the fill method
    flat_vars = {
        var1_name: ak.flatten(obj_var1),
        var2_name: ak.flatten(obj_var2),
        cat1_name: ak.flatten(obj_cat1),
        cat2_name: ak.flatten(obj_cat2)
    }

    hist = (
        dah.Hist.new
        .Variable(var1_binning, name=var1_name, label=var1_name)
        .Variable(var2_binning, name=var2_name, label=var2_name)
        .IntCat(cat1_binning, name=cat1_name)
        .IntCat(cat2_binning, name=cat2_name)
        .Double()
    )

    # deconstruct dict, it just works
    hist.fill(**flat_vars) 

    return hist

def make_2d2d_hist_cat_1reg(
    obj,
    var1_binning,
    var2_binning,
    cat1_binning,
    cat2_binning,
    var1_name="pt",
    var2_name="eta",
    cat1_name="gen_tag",
    cat2_name="qual_tag",
    var1_abs=False,
    var2_abs=False
    ):
    
    """
    Create a 2D x 2D histogram with two numeric variables and two integer categories (typically gen and quality).
    """

    # Extract variables and categories from the object
    obj_var1 = getattr(obj, var1_name)
    obj_var2 = getattr(obj, var2_name)
    obj_cat1 = getattr(obj, cat1_name)
    obj_cat2 = getattr(obj, cat2_name)

    # take absolute value if needed (like for eta or phi, maybe)
    if var1_abs:
        obj_var1 = np.abs(obj_var1)
    if var2_abs:
        obj_var2 = np.abs(obj_var2)

    # flatten arrays, store in dict to be deconstructed later in the fill method
    flat_vars = {
        var1_name: ak.flatten(obj_var1),
        var2_name: ak.flatten(obj_var2),
        cat1_name: ak.flatten(obj_cat1),
        cat2_name: ak.flatten(obj_cat2)
    }

    bins, start, stop = var2_binning
    hist = (
        dah.Hist.new
        .Variable(var1_binning, name=var1_name, label=var1_name)
        .Reg(bins, start, stop, name=var2_name, label=var2_name)
        .IntCat(cat1_binning, name=cat1_name)
        .IntCat(cat2_binning, name=cat2_name)
        .Double()
    )

    # deconstruct dict, it just works
    hist.fill(**flat_vars)

    return hist

def make_2d2d_hist_cat_2reg(
    obj,
    var1_binning,
    var2_binning,
    cat1_binning,
    cat2_binning,
    var1_name="pt",
    var2_name="eta",
    cat1_name="gen_tag",
    cat2_name="qual_tag",
    var1_abs=False,
    var2_abs=False
    ):
    
    """
    Create a 2D x 2D histogram with two numeric variables and two integer categories (typically gen and quality).
    """

    # Extract variables and categories from the object
    obj_var1 = getattr(obj, var1_name)
    obj_var2 = getattr(obj, var2_name)
    obj_cat1 = getattr(obj, cat1_name)
    obj_cat2 = getattr(obj, cat2_name)

    # take absolute value if needed (like for eta or phi, maybe)
    if var1_abs:
        obj_var1 = np.abs(obj_var1)
    if var2_abs:
        obj_var2 = np.abs(obj_var2)

    # flatten arrays, store in dict to be deconstructed later in the fill method
    flat_vars = {
        var1_name: ak.flatten(obj_var1),
        var2_name: ak.flatten(obj_var2),
        cat1_name: ak.flatten(obj_cat1),
        cat2_name: ak.flatten(obj_cat2)
    }

    bins1, start1, stop1 = var1_binning
    bins2, start2, stop2 = var2_binning
    
    hist = (
        dah.Hist.new
        .Reg(bins1, start1, stop1, name=var1_name, label=var1_name)
        .Reg(bins2, start2, stop2, name=var2_name, label=var2_name)
        .IntCat(cat1_binning, name=cat1_name)
        .IntCat(cat2_binning, name=cat2_name)
        .Double()
    )

    # deconstruct dict, it just works
    hist.fill(**flat_vars)

    return hist