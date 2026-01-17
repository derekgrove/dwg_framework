# 2D efficiency plot for lepton pt and eta

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import mplhep
import os


def plot_eta_pt_eff(
    hist_num,
    hist_denom,
    source,
    sample_name="test",
    save_name="test",
    color="cividis",
    save_fig=False,
    plot_directory = "plots",
    vmin=0, vmax=1,
    event_count=None,
    muon_hist = False,
    com=13.6
):

    """
    takes 2 2D histograms of same binning and calculates their efficiencies/errs, plots it and labels each bin.

    This function uses whatever bins and edges are stored in the histogram, so if it looks messed up you probably have too many bins

    I had to change the functionality slightly for muons, we skip the overlap region of tracker and ecal endcap (?) in the middle of the two bins, maybe there is a clever way to have it do this automatically, idk

    example binning for electrons:

    'pt_eta_hist': Hist(
    Variable([2, 3, 4, 5, 7, 10, 20, 45, 75, 1000], name='pt'),
    Variable([0, 0.8, 1.4442, 1.556, 2.5], name='eta'),
    IntCategory([-10, 1, 10, 100, 1000], name='gen_tag'),
    IntCategory([-1, 1, 10, 100], name='qual_tag'),
    storage=Double()) # Sum: 435009.0 (480359.0 with flow)

    example binning muons:

    'pt_eta_hist_muon': Hist(
    Variable([2, 3, 4, 5, 7, 10, 20, 45, 75, 1000], name='pt'),
    Variable([0, 0.9, 1.2, 2.1, 2.4], name='eta'),
    IntCategory([-10, 1, 10, 100, 1000], name='gen_tag'),
    IntCategory([-1, 1, 10, 100], name='qual_tag'),
    storage=Double()) # Sum: 430443.0 (480359.0 with flow)

    but then select the ones you want to put into it like:
    

    numerator = hist[:, :, 1j, 100j]
    
    denominator = hist.integrate("qual", [1j, 10j, 100j])[:, :, 1j]

    sample_name = sample_display_names.get(key, key[:20])
    
    file_name = sample_file_names.get(key, key[:20])

    
    plot_eta_pt_eff(
        numerator, 
        denominator,
        "",
        sample_name=sample_name,
        save_name=f"{file_name}_signal_ele",
        color='plasma',
        vmin=0,
        vmax=1,
        event_count=None,
        com=13 
    )
    
    """
    
    
    
    fig, ax = plt.subplots(figsize=(20, 20))
    mplhep.style.use(mplhep.style.CMS)
    #mplhep.cms.label(loc=0, fontsize=35, com=com)
    #mplhep.cms.label(loc=0, fontsize=35, com=com)
    mplhep.cms.text("Work in Progress", fontsize=30, loc=0)
    eff = two_d_eff_err(hist_num, hist_denom)
    print(eff)
    pt_edges = hist_num.axes['pt'].edges
    pt_labels = [fr"{pt_edges[i]:.3g} - {pt_edges[i+1]:.4g}" for i in range(len(pt_edges) - 1)]
    pt_range = range(len(pt_edges))
    print(pt_edges)
    print(pt_labels)
    print(pt_range)
    
    eta_edges = hist_num.axes['eta'].edges
    
    if muon_hist:
        # Don't skip any bins for muon histograms
        eta_labels = [fr"{eta_edges[i]:.3g} $< |\eta| \leq$ {eta_edges[i+1]:.3g}" 
                      for i in range(len(eta_edges) - 1)]
        eta_indices = list(range(len(eta_edges) - 1))
        eff_filtered = eff[0]
        errors_filtered = eff[1]
    else:
        # Filter out the 1.4442-1.556 bin (index 2)
        skip_bin = 2
        eta_labels = [fr"{eta_edges[i]:.3g} $< |\eta| \leq$ {eta_edges[i+1]:.3g}" 
                      for i in range(len(eta_edges) - 1) if i != skip_bin]
        eta_indices = [i for i in range(len(eta_edges) - 1) if i != skip_bin]
        
        # Filter the efficiency array to exclude the skipped bin
        eff_filtered = np.delete(eff[0], skip_bin, axis=1)
        errors_filtered = np.delete(eff[1], skip_bin, axis=1)
    
    print(eta_edges)
    print(eta_labels)
    print(eta_indices)
    
    plt.imshow(
        eff_filtered,
        aspect='auto',
        cmap=color,
        origin='lower',
        vmin=vmin,
        vmax=vmax
    )
    
    plt.xticks(ticks=range(len(eta_labels)), labels=eta_labels, fontsize=30)
    plt.yticks(ticks=range(len(pt_labels)), labels=pt_labels, fontsize=30)
    
    # Get the efficiency values and errors (filtered)
    values = np.nan_to_num(eff_filtered, nan=0.0)
    errors_filtered = np.nan_to_num(errors_filtered, nan=0.0)
    
    for i in range(len(pt_labels)):
        for j in range(len(eta_labels)):
            val = values[i, j]
            err = errors_filtered[i, j]
            plt.text(j, i, 
                    f"{val:.3f} ± {err:.3f}", 
                    ha='center', va='center', 
                    color='white', fontsize=38,
                    path_effects=[
                    path_effects.Stroke(linewidth=4, foreground='black'),
                    path_effects.Normal()
            ])
    
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    plt.tick_params(axis='x', pad=15)
    plt.tick_params(axis='y', pad=15)
    plt.ylabel("$p_T$ (GeV)", fontsize=50)
    if event_count is None:
        plt.text(-1, 2.4, f"Sample: {sample_name} {source}", ha='center', rotation=90, va='center', color='black', fontsize=35)
    if event_count is not None:
        plt.text(-1.3, 2.4, f"Sample: {sample_name} {source}", ha='center', rotation=90, va='center', color='black', fontsize=35)
        plt.text(-1, 2.4, f"num Events: {event_count}", ha='center', rotation=90, va='center', color='black', fontsize=35)
    
    plt.colorbar()

    if save_fig:
        os.makedirs(plot_directory, exist_ok=True)
        plt.savefig(f"{plot_directory}/{save_name}.pdf", bbox_inches="tight")
    
    plt.show()

    return

def two_d_eff_err(h_num, h_denom): # h_num and h_denom must have same binning in both dimensions

    eff_h = h_num/h_denom # Creating a histogram
    
    err_h = np.sqrt(eff_h.values() * (1 - eff_h.values())/ h_denom.values())

    return eff_h, err_h