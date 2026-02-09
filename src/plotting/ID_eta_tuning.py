from matplotlib import colors
import matplotlib.pyplot as plt
import mplhep
import numpy as np
import os
import hist


def make_AN_ID_eta_plot(
    signal_hist,
    background_hist,
    name="default_name", 
    title=None, 
    plot_txt=None, 
    savefig=False,
    plot_directory = "plots",
    wip=True,
    com=13.6,
    log_color_scale=True,
):
    
    """
    This function should be handed two histograms (of same shape) like:
    (Regular(40, 1, 5, name='ID', label='LowPtElectron ID'),
    Variable([-2.5, -1.4442, -0.8, 0, 0.8, 1.4442, 2.5], name='eta', label='LowPtElectron $\\eta$'))
    """

    norm_sig = signal_hist / signal_hist.sum(flow=False)

    norm_bg = background_hist / background_hist.sum(flow=False)
    
    # Compute ratio safely
    s_o_b = norm_sig / norm_bg
    
    # Extract values
    vals = s_o_b.values(flow=False)
    
    # Clean NaN / inf / negative (log-safe)
    vals = np.nan_to_num(vals, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Rebuild histogram with same axes
    s_o_b_clean = hist.Hist(*s_o_b.axes, storage=s_o_b.storage_type())
    s_o_b_clean[...] = vals
    
    fig, ax = plt.subplots(figsize=(16, 12))
    
    if wip:
        mplhep.cms.text("Work in Progress", fontsize=40, loc=0)
    else:
        mplhep.cms.label(loc=0, fontsize=30, com=com)

    if log_color_scale:
        s_o_b_clean.project("eta", "ID").plot(norm=colors.LogNorm())
    else:
        s_o_b_clean.project("eta", "ID").plot()

    ax.set_xlabel("LowPtElectron $\eta$", fontsize=45, labelpad=15)
    ax.set_ylabel("LowPtElectron ID", fontsize=45, labelpad=20)
    
    ax.tick_params(
    axis="y",
    which="major",
    labelsize=40,
    length=10,
    width=2
    )

    ax.tick_params(
    axis="x",
    which="major",
    labelsize=35,
    length=10,
    width=2,
    rotation=35
    )


    if title is not None:
        ax.text(
            0.85, 1.01, title,
            transform=ax.transAxes,
            fontsize=45,
            ha="center",
            va="bottom"
        )
    
    cbar = plt.gca().collections[0].colorbar
    cbar.set_label("Signal / Background", fontsize=45, labelpad=15)
    cbar.ax.tick_params(labelsize=35)

    if savefig:
        os.makedirs(plot_directory, exist_ok=True)
        plt.savefig(f"{plot_directory}/{name}.pdf",  bbox_inches='tight')

    plt.close(fig)
    
    return ax