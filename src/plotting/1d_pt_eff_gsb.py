import matplotlib.pyplot as plt
import mplhep
import numpy as np


def plot_1d_pt_gsb_eff(
    hist_sig_num,
    hist_sig_denom,
    hist_bg_num,
    hist_bg_denom,
    source,
    title,
    xlabel="Electron $p_T$ (GeV)",
    sample_name="test",
    save_name="test",
    ymin=0, ymax=1,
    event_count=None,
):

    """
    written for a histogram of shape:

    'pt_eta_hist': Hist(
    Variable([2, 3, 4, 5, 7, 10, 20, 45, 75, 1000], name='pt'),
    Variable([0, 0.8, 1.4442, 1.556, 2.5], name='eta'),
    IntCategory([-10, 1, 10, 100, 1000], name='gen_tag'),
    IntCategory([-1, 1, 10, 100], name='qual_tag'),
    storage=Double()) # Sum: 435009.0 (480359.0 with flow)
    """
    fig, ax = plt.subplots(figsize=(30, 12))
    
    indx = 100
    
    xes = np.arange(len(sig_gold_eff[0][:indx])) + 0.5
    
    trim = 2
    
    ax.set_ylim(ymin, ymax)
    
    plt.xlim(0, 100)
    plt.xticks(np.arange(0, 100, 5))
    
    my_ms = 8
    my_cs = 5

    sig_eff, sig_err = calc_eff_err(hist_sig_num, hist_sig_denom)
    fake_eff, fake_err = calc_eff_err(hist_bg_num, hist_bg_denom)

    
    plt.errorbar(x=xes[trim:], y=sig_eff[trim:indx], yerr=sig_err[trim:indx], fmt='^', markersize=my_ms, capsize=my_cs, color='orange', label="Signal gold")
    plt.errorbar(x=xes[trim:], y=sig_silver_eff[0][trim:indx], yerr=sig_silver_eff[1][trim:indx], fmt='^', markersize=my_ms, capsize=my_cs, color='blue', label="Signal silver")
    plt.errorbar(x=xes[trim:], y=sig_bronze_eff[0][trim:indx], yerr=sig_bronze_eff[1][trim:indx], fmt='^', markersize=my_ms, capsize=my_cs, color='brown', label="Signal bronze")
    
    plt.errorbar(x=xes[trim:], y=lf_gold_eff[0][trim:indx], yerr=lf_gold_eff[1][trim:indx], fmt='o', markersize=my_ms, capsize=my_cs, color='orange', label="Light Fakes gold")
    plt.errorbar(x=xes[trim:], y=lf_silver_eff[0][trim:indx], yerr=lf_silver_eff[1][trim:indx], fmt='o', markersize=my_ms, capsize=my_cs, color='blue', label="Light Fakes silver")
    plt.errorbar(x=xes[trim:], y=lf_bronze_eff[0][trim:indx], yerr=lf_bronze_eff[1][trim:indx], fmt='o', markersize=my_ms, capsize=my_cs, color='brown', label="Light Fakes bronze")
    
    
    plt.xlabel(xlabel, fontsize=40)
    plt.ylabel("Efficiency", fontsize=40)
    plt.ylim(0, 1)
    #plt.legend(fontsize=25)
    plt.axvline(x=7, color='gray', linestyle='--', linewidth=2)
    mplhep.cms.label(loc=0, fontsize=30, com=13.6)
    plt.title(title, fontsize=50, pad=40)
    
    # Metadata annotations
    entries = signal_entries
    
    plt.tick_params(axis='y', which='both', right=True, labelright=True)
    plt.grid(visible=None, which='major', axis='both')
    
    #plt.text(1, -0.09, f"Signal objects: {slepsnu_entries}", fontsize=25) #I don't think these are necessary?
    #plt.text(1, -0.14, fr"Fake objects: {ttbar_entries}", fontsize=25) # number of electrons across all events isn't useful information, other than to prove I have high stats, but that is also evident in my error bars
    fig.legend(loc='lower center', fontsize=25, ncol=3, bbox_to_anchor=(0.5, -0.05))
    
    # Save and show
    plt.savefig(savename, bbox_inches='tight')
    plt.show()

    return
    


def calc_eff_err(hist_1, hist_2): #hist_2 is the denominator of the efficiency

    
    num = hist_1.values()
    denom = hist_2.values()
    
    eff = num/denom
    err = np.sqrt(eff * (1 - eff)/ denom)

    filtered_eff = np.nan_to_num(eff, nan=0).tolist()
    filtered_err = np.nan_to_num(err, nan=0).tolist()
        
    return filtered_eff, filtered_err




    