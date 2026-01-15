import matplotlib.pyplot as plt
import mplhep
import numpy as np


def make_AN_1d_pt_eff(signal, fakes, name="default_name", title=None, plot_txt=None, savefig=False):

    from matplotlib.lines import Line2D

    fig, ax = plt.subplots(figsize=(20, 12))

    # signal and fakes are 4-axis: [pt, pt_long, eta, qual]
    # Integrate over pt and eta to get [pt_long, qual] for plotting
    
    # Integrate over pt (axis 0) and eta (axis 2)
    sig_integrated = signal.integrate("pt").integrate("eta")  # Now [pt_long, qual]
    fake_integrated = fakes.integrate("pt").integrate("eta")  # Now [pt_long, qual]
    
    # Select quality levels
    sig_gold = sig_integrated[:, 100j]     # [pt_long]
    sig_silver = sig_integrated[:, 10j]    # [pt_long]
    sig_bronze = sig_integrated[:, 1j]     # [pt_long]
    
    fake_gold = fake_integrated[:, 100j]   # [pt_long]
    fake_silver = fake_integrated[:, 10j]  # [pt_long]
    fake_bronze = fake_integrated[:, 1j]   # [pt_long]
    
    # Baseline = sum of bronze, silver, gold
    sig_baseline = (sig_bronze + sig_silver + sig_gold)
    fake_baseline = (fake_bronze + fake_silver + fake_gold)
    
    sig_gold_eff_err = calc_eff_err(sig_gold, sig_baseline) 
    sig_silver_eff_err = calc_eff_err(sig_silver, sig_baseline)
    sig_bronze_eff_err = calc_eff_err(sig_bronze, sig_baseline)

    fake_gold_eff_err = calc_eff_err(fake_gold, fake_baseline)
    fake_silver_eff_err = calc_eff_err(fake_silver, fake_baseline)
    fake_bronze_eff_err = calc_eff_err(fake_bronze, fake_baseline)
    
    ax.set_ylim(0, 1)

    # Get edges from pt_long axis (now axis 0 after integration)
    edges = sig_integrated.axes[0].edges
    
    xes = (edges[:-1] + edges[1:]) * 0.5
    
    my_ms = 14
    my_cs = 11
    
    plt.errorbar(xes, y=sig_gold_eff_err[0], yerr=sig_gold_eff_err[1], fmt='o', markersize=my_ms, capsize=my_cs, color='darkorange')
    plt.errorbar(xes, y=sig_silver_eff_err[0], yerr=sig_silver_eff_err[1], fmt='o', markersize=my_ms, capsize=my_cs, color='dodgerblue')
    plt.errorbar(xes, y=sig_bronze_eff_err[0], yerr=sig_bronze_eff_err[1], fmt='o', markersize=my_ms, capsize=my_cs, color='firebrick')
    
    plt.errorbar(xes, y=fake_gold_eff_err[0], yerr=fake_gold_eff_err[1], fmt='s', mfc='none', markersize=my_ms, capsize=my_cs, color='darkorange')
    plt.errorbar(xes, y=fake_silver_eff_err[0], yerr=fake_silver_eff_err[1], fmt='s', mfc='none', markersize=my_ms, capsize=my_cs, color='dodgerblue')
    plt.errorbar(xes, y=fake_bronze_eff_err[0], yerr=fake_bronze_eff_err[1], fmt='s', mfc='none', markersize=my_ms, capsize=my_cs, color='firebrick')
    
    plt.xlim(0, 100)
    plt.xticks([20, 40, 60, 80], fontsize=40)
    plt.xlabel("Electron $p_T$ [GeV]", fontsize=40)

    plt.tick_params(axis='y', which='both', right=True, labelright=True, size=20)
    
    plt.ylim(0, 1.46)
    plt.yticks(fontsize=40)
    plt.ylabel("Efficiency", fontsize=40)

    plt.grid(visible=None, which='major', axis='both', linewidth=2.0)
    
    mplhep.cms.label(loc=0, fontsize=30, com=13.6)

    if plot_txt is not None:
        plt.text(65, 1.2, plot_txt,
                 fontsize=30,
                 color='black',
                 ha='center',
                 va='bottom',
                 rotation=0)

    if title is not None:
        plt.title(title, fontsize=50, pad=40)

    handles = [
        Line2D([0], [0], marker='s', color='brown', label='Bronze', markersize=25, linestyle=''),
        Line2D([0], [0], marker='s', color='dodgerblue', label='Silver', markersize=25, linestyle=''),
        Line2D([0], [0], marker='s', color='darkorange', label='Gold', markersize=25, linestyle=''),
        Line2D([0], [0], marker='o', color='black', label='Prompt', markersize=my_ms, linestyle=''),
        Line2D([0], [0], marker='s', color='black', label='Misidentified', mfc='none', markersize=my_ms, linestyle='')
    ]
    
    fig.legend(
        handles=handles,
        loc='upper left',
        bbox_to_anchor=(0.2, 0.85),
        fontsize=30,
        ncol=2,
        frameon=True,
        facecolor='white',
        edgecolor='black',
        framealpha=1
    )

    if savefig:
        plt.savefig(f"eff_plot_{name}.pdf", bbox_inches='tight')
    plt.show()
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


    