#some efficiency tools and plotting tools for said efficiency tools

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import mplhep
import numpy as np
import awkward as ak

import hist.dask as dah



def init_plt():
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import mplhep
    mplhep.style.use(mplhep.style.CMS)
    plt.figure()
    mplhep.style.use(mplhep.style.CMS)


def calc_eff_err_old(hist_1, hist_2): #hist_2 is the denominator of the efficiency

    
    num = hist_1.values()
    denom = hist_2.values()
    
    eff = num/denom
    err = np.sqrt(eff * (1 - eff)/ denom)
    
        
    return eff, err


def calc_eff_err(hist_1, hist_2): #hist_2 is the denominator of the efficiency

    
    num = hist_1.values()
    denom = hist_2.values()
    
    eff = num/denom
    err = np.sqrt(eff * (1 - eff)/ denom)

    filtered_eff = np.nan_to_num(eff, nan=0).tolist()
    filtered_err = np.nan_to_num(err, nan=0).tolist()
        
    return filtered_eff, filtered_err



#### 2D EFFs HERE ###

def two_d_eff_err(h_num, h_denom): # h_num and h_denom must have same binning in both dimensions

    eff_h = h_num/h_denom # Creating a histogram
    
    err_h = np.sqrt(eff_h.values() * (1 - eff_h.values())/ h_denom.values())

    return eff_h, err_h



### Plotting Effs ####


def make_1d_eff_plot_even_binning(eff_errs, pt_bins, label, title="Default title", ymin=-0.1, ymax=1.1):

    fig, ax = plt.subplots(figsize=(20, 12))


    """
    makes plot with even-spaced binning, doesn't matter where the real-space between bins is
    """
    
    
    #pt_bins = [2,3,4,5,7,8,10,15,20,30,45,60,75,500]
    
    str_names = [str(num) for num in range(len(pt_bins))]
    
    all_colors = [
    'black', 'red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive',
    'cyan', 'magenta', 'lime', 'teal', 'navy', 'maroon', 'coral', 'gold', 'indigo', 'violet',
    'turquoise', 'crimson', 'chocolate', 'darkgreen', 'darkblue', 'darkred', 'darkorange',
    'darkviolet', 'darkcyan', 'darkmagenta', 'darkgray', 'lightgray', 'lightblue', 'lightgreen',
    'lightcoral', 'lightpink', 'lightsalmon', 'lightseagreen', 'lightskyblue', 'lightsteelblue',
    'mediumblue', 'mediumseagreen', 'mediumslateblue', 'mediumturquoise', 'mediumvioletred',
    'midnightblue', 'orangered', 'orchid', 'palegreen', 'paleturquoise', 'palevioletred',
    'peru', 'plum', 'rosybrown', 'royalblue', 'saddlebrown', 'salmon', 'sandybrown',
    'seagreen', 'sienna', 'skyblue', 'slateblue', 'springgreen', 'steelblue', 'tan', 'tomato'
    ]

    colors = ['darkmagenta','darkblue','lightsteelblue','sienna', 'plum', 'darkcyan']
    #x = np.arange(len(data)) #literally don't see the advantage to using this, range works just fine


    for i, eff_err in enumerate(eff_errs):
        data = eff_err[0]
        errs = eff_err[1]
        #xs = range(len(data))
        xs = np.arange(len(data))
        
        ax.errorbar(
            xs, data,
            xerr=0.5,
            fmt='o',
            #capsize=10,
            elinewidth=3,
            color=colors[i],
            label=labels[i]
        )
    
        ax.errorbar(
            xs, data,
            yerr=errs,
            fmt='none', #if you would like no point, just the error bar
            capsize=10,
            elinewidth=3,
            color=colors[i]
        )
    
    edge_ticks = np.arange(-0.5, len(data)+0.5)
    #ax.set_xticks([i + 1 for i in x]) #no longer needed for step if I use where='mid'
    ax.set_xticks(edge_ticks)
    ax.set_xticklabels([str(n) for n in pt_bins])
    ax.set_title(title, pad=20, fontweight="bold")
    
    ax.tick_params(axis='x', labelsize=30)
    ax.tick_params(axis='y', labelsize=30)
    
    ax.xaxis.grid(True, which='major', linestyle='--', alpha=0.5)
    ax.legend()
    
    ax.set_xlabel("$p_T$ (GeV)")
    ax.set_ylabel("Efficiency")

    ax.axhline(y=0, linewidth=1, linestyle='--', color='0.5')  # dashed grey
    ax.axhline(y=1, linewidth=1, linestyle='--', color='0.5')  # dashed grey

    ax.set_ylim(ymin, ymax)
    #ax.set_xlim(ymin, ymax)
    
    fig.show()

def make_1d_eff_plot_even_test(
    eff_errs,
    pt_bins,
    title="Default title",
    labels=None,
    ymin=-0.1,
    ymax=1.1,
    fig=None,
    ax=None,
):
    """
    Makes a plot with even-spaced binning. If fig/ax are passed, reuse them;
    otherwise create a new figure.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    # Create fig/ax if not provided
    if fig is None or ax is None:
        fig, ax = plt.subplots(figsize=(20, 12))

    if labels is None:
        labels = [f"Series {i}" for i in range(len(eff_errs))]

    colors = [
        'darkmagenta', 'darkblue', 'lightsteelblue',
        'sienna', 'plum', 'darkcyan'
    ]

    for i, eff_err in enumerate(eff_errs):
        data = eff_err[0]
        errs = eff_err[1]
        xs = np.arange(len(data))

        ax.errorbar(
            xs, data,
            xerr=0.5,
            fmt='o',
            elinewidth=3,
            color=colors[i % len(colors)],
            label=labels[i]
        )
        ax.errorbar(
            xs, data,
            yerr=errs,
            fmt='none',
            capsize=10,
            elinewidth=3,
            color=colors[i % len(colors)]
        )

    # Configure only once (otherwise it gets re-set each call)
    if ax.get_title() == "":
        edge_ticks = np.arange(-0.5, len(pt_bins) + 0.5)
        ax.set_xticks(edge_ticks)
        ax.set_xticklabels([str(n) for n in pt_bins])
        ax.set_title(title, pad=20, fontweight="bold")
        ax.tick_params(axis='x', labelsize=30)
        ax.tick_params(axis='y', labelsize=30)
        ax.xaxis.grid(True, which='major', linestyle='--', alpha=0.5)
        ax.set_xlabel("$p_T$ (GeV)")
        ax.set_ylabel("Efficiency")
        ax.axhline(y=0, linewidth=1, linestyle='--', color='0.5')
        ax.axhline(y=1, linewidth=1, linestyle='--', color='0.5')
        ax.set_ylim(ymin, ymax)

    return fig, ax



"""

examples:

from NM's code:
efficiencyinfo = (
    hist.Hist.new
    .Reg(20, 40, 300, name="pt")
    .Reg(4, 0, 2.5, name="abseta")
    .IntCat([0, 4, 5], name="flavor")
    .Bool(name="passWP")
    .Double()
    .fill(
        pt=jets.pt,
        abseta=abs(jets.eta),
        flavor=jets.hadronFlavour,
        passWP=jets.btagDeepFlavB > 0.2783, # UL 2018 medium WP
    )
)

from Hist documentation:
h = Hist(
    axis.IntCategory([3, 1, 2], label="Number"),
    axis.StrCategory(["Teacher", "Police", "Artist"], label="Profession"),
)
# Sort Number axis increasing and Profession axis decreasing
h1 = h.sort("Number").sort("Profession", reverse=True)
"""

"""
def make_plot_pt_ID(hist, title, text, filename):
    plt.title(title, pad=35, fontsize=35)
    vmin = 0.0001
    vmax = 0.01
    

    # Normalize bin contents
    view = hist.view(flow=False)
    norm = view.sum()
    if norm > 0:
        hist.view(flow=False)[:] = view / norm  # overwrite in place

    mplhep.cms.label(loc=0, fontsize=20, com=13.6)
    plt.figtext(0.1, 0.04, text, fontsize=20)
    
    #hist.plot2d(norm=colors.LogNorm(vmin=vmin, vmax=vmax))
    hist.plot2d(norm=colors.LogNorm(vmin=vmin, vmax=vmax))
    #hist.plot2d()

    plt.yticks(range(int(plt.ylim()[0]), int(plt.ylim()[1]) + 1, 2))
    plt.axhline(y=2.2, color='red', linestyle='--', linewidth=3)
    plt.savefig(filename)
    plt.show()

hist = SlepSnu_results['lpte_hists']['two_d']['pt_ID']['blp_test']
make_plot_pt_ID(hist, title="SlepSnu Signal (Baseline w/o ID cut)", text="SlepSnu signal", filename = "SlepSnu_tighter_ID_pt")




import matplotlib.pyplot as plt
import mplhep
from matplotlib.lines import Line2D

# Load histograms
pt_histograms = ele_hists['one_d']['pt']
eta_histograms = ele_hists['one_d']['eta']
pt_hist_items = list(pt_histograms.items())
eta_hist_items = list(eta_histograms.items())

# Format
selected_indices = [0, 1, 3]
custom_labels = ["raw", "gen filtered", "baseline"]
custom_colors = ["purple", "red", "green"]
total_entries = results["total_entries"]

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10), sharey=True)

# Shared title
fig.suptitle(r"Electron Cutflow (SlepSnu Signal)", fontsize=45, y=1.03)

# -------------------- $p_T$ Plot --------------------
ax1.set_yscale('log')
ax1.tick_params(axis='y', which='both', right=True, labelright=True)
ax1.set_xlabel('Electron $p_T$ (GeV)', fontsize=35)
ax1.set_ylabel('Counts', fontsize=35)
ax1.set_title(r"$p_T$", fontsize=40, pad=40)
mplhep.cms.label(loc=0, fontsize=25, com=13.6, ax=ax1)

for i, label, color in zip(selected_indices, custom_labels, custom_colors):
    name, plot = pt_hist_items[i]
    plot.plot(ax=ax1, label=label, linewidth=3, color=color)

# -------------------- $\eta$ Plot --------------------
ax2.set_yscale('log')
ax2.tick_params(axis='y', which='both', right=True, labelright=True)
ax2.set_xlabel('Electron $\eta$', fontsize=35)
ax2.set_title(r"$\eta$ (integrated $p_T$)", fontsize=40, pad=40)
mplhep.cms.label(loc=0, fontsize=25, com=13.6, ax=ax2)

# Add vertical η lines
for x in [-1.479, 1.479, -0.8, 0.8]:
    ax2.axvline(x=x, color='black', linestyle=':', linewidth=1)

for i, label, color in zip(selected_indices, custom_labels, custom_colors):
    name, plot = eta_hist_items[i]
    plot.plot(ax=ax2, label=label, linewidth=3, color=color)

# -------------------- Shared Legend --------------------
handles = [
    Line2D([0], [0], color=color, lw=3, label=label)
    for label, color in zip(custom_labels, custom_colors)
]
fig.legend(handles=handles, loc='lower center', fontsize=30, ncol=3, bbox_to_anchor=(0.5, -0.05))

# Total events text
fig.text(0.75, 0, f"Total events: {total_entries}", fontsize=30, ha='left')

# Layout
plt.tight_layout()
plt.subplots_adjust(top=0.85, bottom=0.2)
plt.savefig("final_presentation_electron_pt_eta_cutflow_shared", bbox_inches="tight")
plt.show()










def test_yield_plot(h_signal, h_fakes):
    h_yield = h_signal.copy()
    
    # Compute the raw ratio
    h_yield.view()[:] = h_signal.view() / h_fakes.view()

    #print(h_signal.view())
    #print(h_fakes.view())
    #print(h_signal.view() / h_fakes.view())
    array_clean = np.nan_to_num(h_yield, nan=0.0, posinf=0.0, neginf=0.0)
    #print(array_clean)

    total = array_clean.view().sum()
    total_num = h_signal.view().sum()
    total_denom = h_fakes.view().sum()
    c = total_num/total_denom
    #c = total_denom/total_num
    print(c)
    print(1/c)
    print(total)
    print(1/total)
    # Normalize the entire histogram to have unit integral
    #total = h_yield.view().sum()
    
    array_clean.view()[:] /= c    
    print("Sum after normalization:", total)  # should be 1.0

    
    return h_yield

    
#plot_key = ["2_3", "3_4", "4_5", "5_6", "6_7", "7_8", "5_8", "2_5"]
plot_key = ["2_3", "3_4", "4_5", "5_6", "6_7", "7_8"]
os.makedirs("sob_id_pt", exist_ok=True)

for key in plot_key:
    print(f"S/B for baseline_{key}")
    test = test_yield_plot(ss_eta_ID_hists[f"baseline_{key}"], tt_eta_ID_hists[f"baseline_{key}"])
    plt.title(f"S/B for Baseline $p_T$ {key} GeV", pad = 35, fontsize = 35)
    vmin = 1e-1
    vmax = 1e1
    #tot_e = hist.values(flow=True).sum()
    mplhep.cms.label(loc=0, fontsize=20, com=13.6)
    plt.figtext(0.1, 0.04, r"SlepSnu Signal, $t\bar{t}$ Background", fontsize = 20)
    test.plot2d(norm=colors.LogNorm(vmin=vmin, vmax=vmax))
    #test.plot2d()
    
    plt.savefig(f"sob_id_pt/{key}_ID_pT.png")
    plt.show()

"""


