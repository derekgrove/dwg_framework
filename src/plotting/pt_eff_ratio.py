import matplotlib.pyplot as plt
import mplhep
import numpy as np
import os
from matplotlib.lines import Line2D


def make_1d_pt_eff_ratio(
    hist1,
    hist2,
    name="default_name",
    title=None,
    plot_txt=None,
    savefig=False,
    plot_directory='plots',
    lep_type='Electron',
    com=13.6, 
    max_pt=100, 
    log_scale=False,
    marker_size = 14,
    cap_size = 11,
    hist1_label='hist1',
    hist2_label='hist2'):


    fig, ax = plt.subplots(figsize=(25, 12))

    hist1_int = hist1.integrate("pt").integrate("eta")
    hist2_int = hist2.integrate("pt").integrate("eta")

    # This block is to remove 0-efficiency points, essentially empty values in the histograms
    min_pt = 0
    
    if lep_type == "Electron":
        min_pt = 2 # the binning of my hists conveniently lines up with the pt we are removing
    if lep_type == "Muon":
        min_pt = 3

    hist1_int  = hist1_int[min_pt:, :]
    hist2_int = hist2_int[min_pt:, :]

    hist1_gold = hist1_int[:, 100j]
    hist1_silver = hist1_int[:, 10j]
    hist1_bronze = hist1_int[:, 1j]
    
    hist2_gold = hist2_int[:, 100j]
    hist2_silver = hist2_int[:, 10j]
    hist2_bronze = hist2_int[:, 1j]

    hist1_baseline = (hist1_gold + hist1_silver + hist1_bronze)
    hist2_baseline = (hist2_gold + hist2_silver + hist2_bronze)
    
    hist1_gold_eff_err = calc_eff_err(hist1_gold, hist1_baseline) 
    hist1_silver_eff_err = calc_eff_err(hist1_silver, hist1_baseline)
    hist1_bronze_eff_err = calc_eff_err(hist1_bronze, hist1_baseline)
    
    hist2_gold_eff_err = calc_eff_err(hist2_gold, hist2_baseline)
    hist2_silver_eff_err = calc_eff_err(hist2_silver, hist2_baseline)
    hist2_bronze_eff_err = calc_eff_err(hist2_bronze, hist2_baseline)

    

    gold_ratio, gold_ratio_err = ratio_and_err(
        hist1_gold_eff_err[0], hist1_gold_eff_err[1],
        hist2_gold_eff_err[0], hist2_gold_eff_err[1],
    )
    
    silver_ratio, silver_ratio_err = ratio_and_err(
        hist1_silver_eff_err[0], hist1_silver_eff_err[1],
        hist2_silver_eff_err[0], hist2_silver_eff_err[1],
    )
    
    bronze_ratio, bronze_ratio_err = ratio_and_err(
        hist1_bronze_eff_err[0], hist1_bronze_eff_err[1],
        hist2_bronze_eff_err[0], hist2_bronze_eff_err[1],
    )



    edges = hist1_int.axes[0].edges
    xes = (edges[:-1] + edges[1:]) * 0.5
    
    plt.errorbar(xes, y=gold_ratio, yerr=gold_ratio_err, fmt='o', markersize=marker_size, capsize=cap_size, color='darkorange')
    plt.errorbar(xes, y=silver_ratio, yerr=gold_ratio_err, fmt='o', markersize=marker_size, capsize=cap_size, color='dodgerblue')
    plt.errorbar(xes, y=bronze_ratio, yerr=gold_ratio_err, fmt='o', markersize=marker_size, capsize=cap_size, color='firebrick')
    
    # X-axis settings (always linear)
    plt.xlim(0, max_pt)
    # Adjust xticks based on max_pt
    if max_pt <= 30:
        plt.xticks([5, 10, 15, 20, 25, 30], fontsize=40)
    elif max_pt <= 50:
        plt.xticks([10, 20, 30, 40, 50], fontsize=40)
    else:
        plt.xticks([20, 40, 60, 80, 100], fontsize=40)
    
    plt.xlabel(f"{lep_type} $p_T$ [GeV]", fontsize=40)
    
    # Y-axis settings - handle log scale vs linear scale
    if log_scale:
        plt.yscale('log')
        plt.ylim(0.3, 2.5)  # Log scale needs to avoid 0
    else:
        plt.ylim(0.3, 2.5)
    
    #plt.yticks(fontsize=40)
    plt.ylabel(f"FullSim Efficiencies / FastSim Efficiencies", fontsize=40, labelpad = 20)
    plt.tick_params(axis='both', which='both', right=True, labelright=True, size=30)
    plt.yticks(fontsize=30)
    plt.grid(visible=None, which='major', axis='y', linewidth=2.0)
    
    if lep_type.lower() == "electron":
        plt.axvline(x=7, linestyle='dotted', color='black', linewidth=1.5)
        plt.axvline(x=20, linestyle='dotted', color='black', linewidth=1.5)
    
    mplhep.cms.label(loc=0, fontsize=40, com=com)
    #mplhep.cms.text("Work in Progress", fontsize=40, loc=0)
    
    if plot_txt is not None:
        # Use axes coordinates for positioning
        plt.text(0.7, 1.055, plot_txt,
                 fontsize=35,
                 color='black',
                 ha='center',
                 va='top',
                 transform=ax.transAxes,
                 rotation=0,
                 weight='bold'
                )
    
    if title is not None:
        plt.title(title, fontsize=45, pad=60)

    """
    handles = [
        Line2D([0], [0], marker='s', color='firebrick', label='Bronze', markersize=25, linestyle=''),
        Line2D([0], [0], marker='s', color='dodgerblue', label='Silver', markersize=25, linestyle=''),
        Line2D([0], [0], marker='s', color='darkorange', label='Gold', markersize=25, linestyle=''), 
    ]
    
    fig.legend(
        handles=handles,
        loc='center left',
        bbox_to_anchor=(0.95, 0.5),
        fontsize=30,
        ncol=1,
        frameon=True,
        facecolor='white',
        edgecolor='black',
        framealpha=1
    )
    """
    if savefig:
        os.makedirs(plot_directory, exist_ok=True)
        suffix = "_log" if log_scale else ""
        suffix = suffix + "_" + str(max_pt) + "GeV"
        plt.savefig(f"{plot_directory}/eff_ratio_plot_{name}{suffix}.pdf", bbox_inches='tight')
        
    plt.close()
    return ax



def calc_eff_err(hist_1, hist_2): #hist_2 is the denominator of the efficiency

    
    num = hist_1.values()
    denom = hist_2.values()
    
    eff = num/denom
    err = np.sqrt(eff * (1 - eff)/ denom)

    filtered_eff = np.nan_to_num(eff, nan=0).tolist()
    filtered_err = np.nan_to_num(err, nan=0).tolist()
        
    return filtered_eff, filtered_err

def ratio_and_err(eff1, err1, eff2, err2):
    
    """
    Compute ratio = eff1 / eff2 and its propagated uncertainty.

    All inputs are numpy arrays of the same shape.
    Assumes eff1 and eff2 are independent.
    """

    eff1 = np.array(eff1)
    err1 = np.array(err1)
    eff2 = np.array(eff2)
    err2 = np.array(err2)

    ratio = eff1 / eff2

    ratio_err = ratio * np.sqrt(
        (err1 / eff1)**2 +
        (err2 / eff2)**2
    )

    return ratio, ratio_err