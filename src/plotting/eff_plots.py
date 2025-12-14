from .efficiencies import *

from .histers import *



def make_1d_eff_plot(eff_errs, pt_bins, title="Default title", label=labels, ymin=-0.1, ymax=1.1):

    fig, ax = plt.subplots(figsize=(20, 12))
    
    #pt_bins = [2,3,4,5,7,8,10,15,20,30,45,60,75,500]
    
    str_names = [str(num) for num in range(len(pt_bins))]
    
    colors = [
    #'black', 'red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive',
    'cyan', 'magenta', 'lime', 'teal', 'navy', 'maroon', 'coral', 'gold', 'indigo', 'violet',
    'turquoise', 'crimson', 'chocolate', 'darkgreen', 'darkblue', 'darkred', 'darkorange',
    'darkviolet', 'darkcyan', 'darkmagenta', 'darkgray', 'lightgray', 'lightblue', 'lightgreen',
    'lightcoral', 'lightpink', 'lightsalmon', 'lightseagreen', 'lightskyblue', 'lightsteelblue',
    'mediumblue', 'mediumseagreen', 'mediumslateblue', 'mediumturquoise', 'mediumvioletred',
    'midnightblue', 'orangered', 'orchid', 'palegreen', 'paleturquoise', 'palevioletred',
    'peru', 'plum', 'rosybrown', 'royalblue', 'saddlebrown', 'salmon', 'sandybrown',
    'seagreen', 'sienna', 'skyblue', 'slateblue', 'springgreen', 'steelblue', 'tan', 'tomato'
    ]
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