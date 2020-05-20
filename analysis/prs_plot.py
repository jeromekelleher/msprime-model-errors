# Figure (A) plot demography with no demographic events, (B) Martin figure 5B,
# biases across populations, (C)-(E) Martin figure 5C, with same parameters
# as shown there. 
# To do: see if I have rerun the pipeline with that paper's incorrect demography,
# and plot side by side in (C)-(E) in dashed or alpha < 1, for direct comparison.

import models
import demography
import math
import numpy as np
import matplotlib.pylab as plt
import matplotlib
import matplotlib.font_manager as fm
from bokeh import palettes
import pandas, pickle
import seaborn as sns

cb_colors = palettes.Colorblind8

matplotlib.rcParams['text.usetex'] = True
plt.rcParams['legend.title_fontsize'] = 'xx-small'
matplotlib.rc('xtick', labelsize=6) 
matplotlib.rc('ytick', labelsize=6) 
matplotlib.rc('axes', labelsize=8) 
matplotlib.rc('axes', titlesize=8) 
matplotlib.rc('legend', fontsize=8)


# from https://matplotlib.org/3.1.1/gallery/statistics/customized_violin.html
# to help with plotting violins

def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels) + 1))
    ax.set_xticklabels(labels)
    ax.set_xlim(0.25, len(labels) + 0.75)
    ax.set_xlabel('Sample name')


if __name__ == "__main__":
    
    pc, mm, de, dg = models.martin_model(merge_time=110000)

    # prs simulation data
    fname = f'data/example_sim_nhaps_400000_400000_400000_h2_0.67_m_1000_alpha_0.0.prs.gz'
    data = pandas.read_csv(fname, sep='\t')
    AFR_true = data[data['Pop']=='AFR']['PRS_true']
    AFR_infer = data[data['Pop']=='AFR']['PRS_infer']
    EUR_data = data[data['Pop']=='EUR']  
    EUR_data = EUR_data[EUR_data.isnull().any(axis=1)]
    EUR_true = EUR_data[EUR_data['Pop']=='EUR']['PRS_true']
    EUR_infer = EUR_data[EUR_data['Pop']=='EUR']['PRS_infer']
    EAS_true = data[data['Pop']=='EAS']['PRS_true']
    EAS_infer = data[data['Pop']=='EAS']['PRS_infer']

    # standardize PRSs
    true_mean = np.mean( np.concatenate(( AFR_true, EUR_true, EAS_true )) )
    true_std = np.std( np.concatenate(( AFR_true, EUR_true, EAS_true )) )
    infer_mean = np.mean( np.concatenate(( AFR_infer, EUR_infer, EAS_infer )) )
    infer_std = np.std( np.concatenate(( AFR_infer, EUR_infer, EAS_infer )) )

    AFR_true = (AFR_true - true_mean) / true_std
    EUR_true = (EUR_true - true_mean) / true_std
    EAS_true = (EAS_true - true_mean) / true_std
    AFR_infer = (AFR_infer - infer_mean) / infer_std
    EUR_infer = (EUR_infer - infer_mean) / infer_std
    EAS_infer = (EAS_infer - infer_mean) / infer_std

    correlations = pickle.load(open('data/correlations_alpha_0.0.bp', 'rb'))
    distributions = pickle.load(open('data/distributions_alpha_0.0.bp', 'rb'))
    
    
    fig = plt.figure(2, figsize=(6.5, 4))
    fig.clf()
    
    ax1 = plt.subplot2grid((4, 9), (0, 0), rowspan=2, colspan=4)
    demography.plotting.plot_demography(dg, leaf_order=['YRI', 'CEU', 'CHB'], ax=ax1,
                                        gen=25, stacked=[], flipped=['CEU'],
                                        padding=3)
    ax1.set_ylim(top=100)
    ax1.spines['left'].set_bounds(0, ax1.get_ylim()[1])
    
    ## plot distributions over 100 replicates
    ax2 = plt.subplot2grid((4, 9), (0, 4), colspan=2)
    sns.distplot(distributions['true'][.67][1000]['EUR']['mean'],
                 color=cb_colors[1], hist = False, kde = True,
                 kde_kws = {'shade': True, 'linewidth': 1},
                 ax=ax2)
    sns.distplot(distributions['true'][.67][1000]['EAS']['mean'],
                 color=cb_colors[3], hist = False, kde = True,
                 kde_kws = {'shade': True, 'linewidth': 1},
                 ax=ax2)
    sns.distplot(distributions['true'][.67][1000]['AFR']['mean'],
                 color=cb_colors[0], hist = False, kde = True,
                 kde_kws = {'shade': True, 'linewidth': 1},
                 ax=ax2)
    
    ax2.set_ylabel('Density')
    ax2.set_xlabel('True PRS')
    ax2.set_yticks([0, 0.05])
    ax2.set_xticks([-15, 0, 15])
    
    ax3 = plt.subplot2grid((4, 9), (1, 4), colspan=2)
    sns.distplot(distributions['infer'][.67][1000]['EUR']['mean'],
                 color=cb_colors[1], hist = False, kde = True,
                 kde_kws = {'shade': True, 'linewidth': 1},
                 ax=ax3)
    sns.distplot(distributions['infer'][.67][1000]['EAS']['mean'],
                 color=cb_colors[3], hist = False, kde = True,
                 kde_kws = {'shade': True, 'linewidth': 1},
                 ax=ax3)
    sns.distplot(distributions['infer'][.67][1000]['AFR']['mean'],
                 color=cb_colors[0], hist = False, kde = True,
                 kde_kws = {'shade': True, 'linewidth': 1},
                 ax=ax3)
    
    ax3.set_ylabel('Density')
    ax3.set_xlabel('Inferred PRS')
    ax3.set_yticks([0, 0.02, 0.04])
    ax3.set_xticks([-100, 0, 100])
    
    ## plot example of True vs Inferred individual PRS
    ax4 = plt.subplot2grid((4, 9), (0, 6), rowspan=2, colspan=3)
    ax4.plot(EUR_true[:10000], EUR_infer[:10000], 'o', markersize=.5,
             c=cb_colors[1], label='EUR')
    ax4.plot(EAS_true[:10000], EAS_infer[:10000], 'o', markersize=.5,
             c=cb_colors[3], label='EAS')
    ax4.plot(AFR_true[:10000], AFR_infer[:10000], 'o', markersize=.5,
             c=cb_colors[0], label='AFR')
    
    ax4.set_xlabel(r'\emph{Z} true')
    ax4.set_ylabel(r'\emph{Z} inferred')
    ax4.legend(fontsize=6, frameon=False)
    
    ## violin plots
    h2 = 0.67
    violin_colors = [cb_colors[0], cb_colors[3], cb_colors[1]]

    ax5 = plt.subplot2grid((4, 9), (2, 0), rowspan=2, colspan=3)
    m = 200
    data = [correlations[h2][m][k] for k in ['AFR','EAS','EUR']]
    
    for y in [0.25, 0.5, 0.75]:
        ax5.axhline(y, ls='--', lw=0.5, color='black')
    
    parts = ax5.violinplot(data,showmeans=False, showmedians=False,
                           showextrema=False)
    
    for ii, pc in enumerate(parts['bodies']):
        pc.set_facecolor(violin_colors[ii])
        pc.set_edgecolor('black')
        pc.set_linewidth(1)
        pc.set_alpha(0.6)

    quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=1)
    inds = np.arange(1, len(medians) + 1)
    ax5.scatter(inds, medians, marker=".", color='white', s=15, zorder=3)
    ax5.vlines(inds, quartile1, quartile3, color='black', linestyle='-', lw=5)
    
    quartile1, medians, quartile3 = np.percentile(data, [5, 50, 95], axis=1)
    ax5.vlines(inds, quartile1, quartile3, color='black', linestyle='-', lw=1)

    labels = ['AFR','EAS','EUR']
    set_axis_style(ax5, labels)

    ax5.set_ylim([0,1])
    ax5.set_title(r'${0}$ causal variants'.format(m))
    ax5.set_ylabel("Pearson's correlation")
    ax5.set_xlabel("")
    ax5.set_yticks([0, 0.25, 0.5, 0.75, 1])
    
    # m = 500
    ax6 = plt.subplot2grid((4, 9), (2, 3), rowspan=2, colspan=3)
    m = 500
    data = [correlations[h2][m][k] for k in ['AFR','EAS','EUR']]
    
    for y in [0.25, 0.5, 0.75]:
        ax6.axhline(y, ls='--', lw=0.5, color='black')
    
    parts = ax6.violinplot(data,showmeans=False, showmedians=False,
                           showextrema=False)
    
    for ii, pc in enumerate(parts['bodies']):
        pc.set_facecolor(violin_colors[ii])
        pc.set_edgecolor('black')
        pc.set_linewidth(1)
        pc.set_alpha(0.6)

    quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=1)
    inds = np.arange(1, len(medians) + 1)
    ax6.scatter(inds, medians, marker=".", color='white', s=15, zorder=3)
    ax6.vlines(inds, quartile1, quartile3, color='black', linestyle='-', lw=5)
    
    quartile1, medians, quartile3 = np.percentile(data, [5, 50, 95], axis=1)
    ax6.vlines(inds, quartile1, quartile3, color='black', linestyle='-', lw=1)

    labels = ['AFR','EAS','EUR']
    set_axis_style(ax6, labels)

    ax6.set_ylim([0,1])
    ax6.set_title(r'${0}$ causal variants'.format(m))
    ax6.set_xlabel("")
    ax6.set_yticks([0, 0.25, 0.5, 0.75, 1])

    # m = 1000
    ax7 = plt.subplot2grid((4, 9), (2, 6), rowspan=2, colspan=3)
    m = 1000
    data = [correlations[h2][m][k] for k in ['AFR','EAS','EUR']]
    
    for y in [0.25, 0.5, 0.75]:
        ax7.axhline(y, ls='--', lw=0.5, color='black')
    
    parts = ax7.violinplot(data,showmeans=False, showmedians=False,
                           showextrema=False)
    
    for ii, pc in enumerate(parts['bodies']):
        pc.set_facecolor(violin_colors[ii])
        pc.set_edgecolor('black')
        pc.set_linewidth(1)
        pc.set_alpha(0.6)

    quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=1)
    inds = np.arange(1, len(medians) + 1)
    ax7.scatter(inds, medians, marker=".", color='white', s=15, zorder=3)
    ax7.vlines(inds, quartile1, quartile3, color='black', linestyle='-', lw=5)
    
    quartile1, medians, quartile3 = np.percentile(data, [5, 50, 95], axis=1)
    ax7.vlines(inds, quartile1, quartile3, color='black', linestyle='-', lw=1)

    labels = ['AFR','EAS','EUR']
    set_axis_style(ax7, labels)

    ax7.set_ylim([0,1])
    ax7.set_title(r'${0}$ causal variants'.format(m))
    ax7.set_xlabel("")
    ax7.set_yticks([0, 0.25, 0.5, 0.75, 1])


    fig.text(0.025, 0.98, 'A', fontsize=9, ha='center', va='center')
    fig.text(0.45, 0.98, 'B', fontsize=9, ha='center', va='center')
    fig.text(0.45, 0.75, 'C', fontsize=9, ha='center', va='center')
    fig.text(0.7, 0.98, 'D', fontsize=9, ha='center', va='center')
    fig.text(0.025, 0.45, 'E', fontsize=9, ha='center', va='center')
    fig.text(0.36, 0.45, 'F', fontsize=9, ha='center', va='center')
    fig.text(0.685, 0.45, 'G', fontsize=9, ha='center', va='center')

    fig.tight_layout()
    plt.savefig('../figures/prs_fig.pdf')
    plt.show()
