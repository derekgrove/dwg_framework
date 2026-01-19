import matplotlib.pyplot as plt
import mplhep
import numpy as np
import os


def make_AN_1d_pt_eff(signal, 
                      fakes, 
                      name="default_name", 
                      title=None, 
                      plot_txt=None, 
                      savefig=False,
                      plot_directory = "plots",
                      lepton_type="Electron",
                      com=13.6, 
                      max_pt=100, 
                      log_scale=False,
                      marker_size = 14,
                      cap_size = 11,
                      fake_label="Misidentified"):
    
    from matplotlib.lines import Line2D
    fig, ax = plt.subplots(figsize=(25, 12))
    # signal and fakes are 4-axis: [pt, pt_long, eta, qual]
    # Integrate over pt and eta to get [pt_long, qual] for plotting
    
    sig_integrated = signal.integrate("pt").integrate("eta")
    fake_integrated = fakes.integrate("pt").integrate("eta")

    
    # This block is to remove 0-efficiency points, essentially empty values in the histograms
    min_pt = 0
    
    if lepton_type == "Electron":
        min_pt = 2 # the binning of my hists conveniently lines up with the pt we are removing
    if lepton_type == "Muon":
        min_pt = 3

    sig_integrated  = sig_integrated[min_pt:, :]
    fake_integrated = fake_integrated[min_pt:, :]
    
    #end block of removing 0 eff leptons
        
    
    sig_gold = sig_integrated[:, 100j]
    sig_silver = sig_integrated[:, 10j]
    sig_bronze = sig_integrated[:, 1j]
    
    fake_gold = fake_integrated[:, 100j]
    fake_silver = fake_integrated[:, 10j]
    fake_bronze = fake_integrated[:, 1j]
    
    sig_baseline = (sig_bronze + sig_silver + sig_gold)
    fake_baseline = (fake_bronze + fake_silver + fake_gold)
    
    sig_gold_eff_err = calc_eff_err(sig_gold, sig_baseline) 
    sig_silver_eff_err = calc_eff_err(sig_silver, sig_baseline)
    sig_bronze_eff_err = calc_eff_err(sig_bronze, sig_baseline)
    fake_gold_eff_err = calc_eff_err(fake_gold, fake_baseline)
    fake_silver_eff_err = calc_eff_err(fake_silver, fake_baseline)
    fake_bronze_eff_err = calc_eff_err(fake_bronze, fake_baseline)
    
    edges = sig_integrated.axes[0].edges
    xes = (edges[:-1] + edges[1:]) * 0.5
    
    plt.errorbar(xes, y=sig_gold_eff_err[0], yerr=sig_gold_eff_err[1], fmt='o', markersize=marker_size, capsize=cap_size, color='darkorange')
    plt.errorbar(xes, y=sig_silver_eff_err[0], yerr=sig_silver_eff_err[1], fmt='o', markersize=marker_size, capsize=cap_size, color='dodgerblue')
    plt.errorbar(xes, y=sig_bronze_eff_err[0], yerr=sig_bronze_eff_err[1], fmt='o', markersize=marker_size, capsize=cap_size, color='firebrick')
    
    plt.errorbar(xes, y=fake_gold_eff_err[0], yerr=fake_gold_eff_err[1], fmt='s', mfc='none', markersize=marker_size, capsize=cap_size, color='darkorange')
    plt.errorbar(xes, y=fake_silver_eff_err[0], yerr=fake_silver_eff_err[1], fmt='s', mfc='none', markersize=marker_size, capsize=cap_size, color='dodgerblue')
    plt.errorbar(xes, y=fake_bronze_eff_err[0], yerr=fake_bronze_eff_err[1], fmt='s', mfc='none', markersize=marker_size, capsize=cap_size, color='firebrick')
    
    # X-axis settings (always linear)
    plt.xlim(0, max_pt)
    # Adjust xticks based on max_pt
    if max_pt <= 30:
        plt.xticks([5, 10, 15, 20, 25, 30], fontsize=40)
    elif max_pt <= 50:
        plt.xticks([10, 20, 30, 40, 50], fontsize=40)
    else:
        plt.xticks([20, 40, 60, 80, 100], fontsize=40)
    
    plt.xlabel(f"{lepton_type} $p_T$ [GeV]", fontsize=40)
    
    # Y-axis settings - handle log scale vs linear scale
    if log_scale:
        plt.yscale('log')
        plt.ylim(0.001, 1.0)  # Log scale needs to avoid 0
    else:
        plt.ylim(0, 1.0)
    
    #plt.yticks(fontsize=40)
    plt.ylabel(f"Quality Assignment Efficiency", fontsize=40, labelpad = 20)
    plt.tick_params(axis='both', which='both', right=True, labelright=True, size=20)
    plt.yticks(fontsize=30)
    plt.grid(visible=None, which='major', axis='y', linewidth=2.0)
    
    if lepton_type.lower() == "electron":
        plt.axvline(x=7, linestyle='dotted', color='black', linewidth=1.5)
    
    #mplhep.cms.label(loc=0, fontsize=30, com=com)
    mplhep.cms.text("Work in Progress", fontsize=30, loc=0)
    
    if plot_txt is not None:
        # Use axes coordinates for positioning
        plt.text(0.7, 0.95, plot_txt,
                 fontsize=30,
                 color='black',
                 ha='center',
                 va='top',
                 transform=ax.transAxes,
                 rotation=0)
    
    if title is not None:
        plt.title(title, fontsize=45, pad=40)
    
    handles = [
        Line2D([0], [0], marker='s', color='brown', label='Bronze', markersize=25, linestyle=''),
        Line2D([0], [0], marker='s', color='dodgerblue', label='Silver', markersize=25, linestyle=''),
        Line2D([0], [0], marker='s', color='darkorange', label='Gold', markersize=25, linestyle=''),
        Line2D([0], [0], marker='o', color='black', label='Prompt', markersize=25, linestyle=''),
        Line2D([0], [0], marker='s', color='black', label=fake_label, mfc='none', markersize=25, linestyle='')
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
    
    if savefig:
        os.makedirs(plot_directory, exist_ok=True)
        suffix = "_log" if log_scale else ""
        suffix = suffix + "_" + str(max_pt) + "GeV"
        plt.savefig(f"{plot_directory}/eff_plot_{name}{suffix}.pdf", bbox_inches='tight')
    return fig, ax



def calc_eff_err(hist_1, hist_2): #hist_2 is the denominator of the efficiency

    
    num = hist_1.values()
    denom = hist_2.values()
    
    eff = num/denom
    err = np.sqrt(eff * (1 - eff)/ denom)

    filtered_eff = np.nan_to_num(eff, nan=0).tolist()
    filtered_err = np.nan_to_num(err, nan=0).tolist()
        
    return filtered_eff, filtered_err


"""

This one now has problems with sqrt of a negative number

def calc_eff_err(hist_1, hist_2):  # hist_2 is the denominator of the efficiency
    
    num = hist_1.values()
    denom = hist_2.values()
    
    # Avoid division by zero warning
    eff = np.divide(num, denom, out=np.zeros_like(num, dtype=float), where=denom!=0)
    
    # Calculate error only where denom > 0 to avoid sqrt of negative or division by zero
    err = np.where(denom > 0, np.sqrt(eff * (1 - eff) / denom), 0)
    
    filtered_eff = eff.tolist()
    filtered_err = err.tolist()
        
    return filtered_eff, filtered_err
"""


    