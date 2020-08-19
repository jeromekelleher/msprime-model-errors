# Here we compute the AFS and LD decay curves for each population. To compute these
# statistics, we use moments and use demography to define the OOA demographic models,
# which also wraps the moments functions to compute the AFS or LD decay curves.
# We also compute the coalescent rate trajectories for each model and compare to the
# true population size history under the model. Finally, we print Fst for each pair
# of populations under both the correct and incorrect model, as well as expected
# heterozygosites in each population. Those values are summaries of the joint AFS.


import msprime
import models
import demography
import networkx as nx
import math
import moments, moments.LD as mold
import numpy as np
import matplotlib.pylab as plt
import matplotlib
import matplotlib.font_manager as fm
from bokeh import palettes

cb_colors = palettes.Colorblind8

matplotlib.rcParams['text.usetex'] = True
plt.rcParams['legend.title_fontsize'] = 'xx-small'
matplotlib.rc('xtick', labelsize=6) 
matplotlib.rc('ytick', labelsize=6) 
matplotlib.rc('axes', labelsize=8) 
matplotlib.rc('axes', titlesize=8) 
matplotlib.rc('legend', fontsize=8)


if __name__ == "__main__":
    dg = models.ooa_dg()
    dg_mig = models.ooa_mig_dg() # for plotting. we compute using a more stable function

    ## get allele frequency spectra using moments (Jouganous et al 2017)
    num_sam = 100 # number of samples per population
    try:
        fs = moments.Spectrum.from_file(f'data/ooa.{num_sam}.fs')
    except IOError:
        fs = dg.SFS(['YRI', 'CEU', 'CHB'], [num_sam, num_sam, num_sam])
        fs.to_file(f'data/ooa.{num_sam}.fs')
    try:
        fs_mig = moments.Spectrum.from_file(f'data/ooa.ancient_mig.{num_sam}.fs')
    except IOError:
        fs_mig = models.get_sfs_ancient_migration(num_sam)
        fs_mig.to_file(f'data/ooa.ancient_mig.{num_sam}.fs')
    
    # print Heterozygosities and Fst
    print(f"Statistics computed with {num_sam} sampled chromosomes from each population")
    print("Heterozygosity (theta=1):")
    print("pop\tcorrect\tincorrect")
    print(f"YRI\t{fs.marginalize([1,2]).pi():.4}\t{fs_mig.marginalize([1,2]).pi():.4}")
    print(f"CEU\t{fs.marginalize([0,2]).pi():.4}\t{fs_mig.marginalize([0,2]).pi():.4}")
    print(f"CHB\t{fs.marginalize([0,1]).pi():.4}\t{fs_mig.marginalize([0,1]).pi():.4}")
    print("Fst:")
    print("pops\tcorrect\tincorrect")
    print(f"YRI-CEU\t{fs.marginalize([2]).Fst():.4}\t{fs_mig.marginalize([2]).Fst():.4}")
    print(f"YRI-CHB\t{fs.marginalize([1]).Fst():.4}\t{fs_mig.marginalize([1]).Fst():.4}")
    print(f"CEU-CHB\t{fs.marginalize([0]).Fst():.4}\t{fs_mig.marginalize([0]).Fst():.4}")
    
    ## get LD statistics using moments.LD (Hill Robertson 1968, Ragsdale & Gravel 2019)
    rs = np.logspace(-6, -3, 21)
    rhos = 4 * dg.Ne * rs
    y = dg.LD(['YRI', 'CEU', 'CHB'], rho=rhos)
    y_mig = models.get_ld_ancient_migration(rho=rhos)

    D2_yri = [y[ii][0] for ii,_ in enumerate(rs)]
    D2_chb = [y[ii][5] for ii,_ in enumerate(rs)]
    D2_mig_yri = [y_mig[ii][0] for ii,_ in enumerate(rs)]
    D2_mig_chb = [y_mig[ii][5] for ii,_ in enumerate(rs)]

    ## get coalescence rates
    generation_time = 25
    T_AF = 220e3 / generation_time
    T_B = 140e3 / generation_time
    T_EU_AS = 21.2e3 / generation_time

    gens = int(1e6/generation_time)

    coal_rates = dg.coal_rates(['YRI', 'CEU', 'CHB'], gens=gens)
    coal_rates_mig = dg_mig.coal_rates(['YRI', 'CEU', 'CHB'], gens=gens)

    model = models.out_of_africa()
    model_mig = models.out_of_africa_ancient_migration()
    dd = msprime.DemographyDebugger(
        population_configurations=model[0],
        migration_matrix=model[1],
        demographic_events=model[2]
    )
    dd_mig = msprime.DemographyDebugger(
        population_configurations=model_mig[0],
        migration_matrix=model_mig[1],
        demographic_events=model_mig[2]
    )

    steps = np.linspace(0, gens-1, gens)
    true_sizes = dd.population_size_trajectory(steps = steps)
    yri_size = true_sizes[:,0]
    ceu_size = np.concatenate((
        true_sizes[:int(T_B), 1],
        true_sizes[int(T_B):, 0]
    ))
    chb_size = np.concatenate((
        true_sizes[:int(T_EU_AS), 2],
        true_sizes[int(T_EU_AS):int(T_B), 1],
        true_sizes[int(T_B):, 0]
    ))

    # population colors, to match Martin et al
    yri_color = palettes.Paired[10][5]
    ceu_color = palettes.Paired[10][1]
    chb_color = palettes.Paired[10][9]

    fig = plt.figure(1, figsize=(6.5, 6))
    fig.clf()

    # plot the two models to make sure they look right
    ax1 = plt.subplot2grid((3, 6), (0, 0), colspan=3)
    demography.plotting.plot_demography(dg, leaf_order=['YRI', 'CEU', 'CHB'], ax=ax1,
                                        gen=25, stacked=[('A', 'YRI')], flipped=['CEU'],
                                        root_length=2, padding=3)
    ax1.set_ylim([0, 300])
    ax1.spines['left'].set_bounds(0, ax1.get_ylim()[1])

    ax2 = plt.subplot2grid((3, 6), (0, 3), colspan=3)
    demography.plotting.plot_demography(dg_mig, leaf_order=['YRI', 'CEU', 'CHB'], ax=ax2,
                                        gen=25, stacked=[('A', 'YRI')], flipped=['CEU'],
                                        padding=3)
    ax2.set_ylim([0, 300])
    ax2.spines['left'].set_bounds(0, ax2.get_ylim()[1])
    ax2.set_ylabel('')

    ax3 = plt.subplot2grid((3, 6), (1, 0), colspan=3)
    ax3.plot(fs.marginalize([1,2]),
             '-', lw=1, color=yri_color, label='YRI (model A)')
    ax3.plot(fs_mig.marginalize([1,2]),
             '--', lw=1, color=yri_color, label='YRI (model B)')
    
    ax3.plot(fs.marginalize([0,2]),
             '-', lw=1, color=ceu_color, label='CEU (A)')
    ax3.plot(fs_mig.marginalize([0,2]),
             '--', lw=1, color=ceu_color, label='CEU (B)')

    ax3.set_yscale('log')
    ax3.set_xlabel('Allele frequency')
    ax3.set_ylabel(r'Count (with $\theta=1$)')
    ax3.legend(frameon=False, fontsize=6)

    ax4 = plt.subplot2grid((3, 6), (1, 3), colspan=3)
    ax4.plot(1e2*rs, D2_yri, '-', color=yri_color, lw=1, label='YRI (A)')
    ax4.plot(1e2*rs, D2_mig_yri, '--', color=yri_color, lw=1, label='YRI (B)')
    ax4.plot(1e2*rs, D2_chb, '-', color=chb_color, lw=1, label='CHB (A)')
    ax4.plot(1e2*rs, D2_mig_chb, '--', color=chb_color, lw=1, label='CHB (B)')

    ax4.set_yscale('log')
    ax4.set_ylabel('$E[D^2]$')
    ax4.set_xlabel('cM')
    ax4.set_xscale('log')
    ax4.legend(frameon=False, fontsize=6)

    ax5 = plt.subplot2grid((3, 6), (2, 0), colspan=2)
    
    ax5.plot(steps*generation_time, yri_size, ':', color='gray', label='YRI size')
    ax5.plot(steps*generation_time, 1 / 2 / coal_rates[('YRI', 'YRI')],
             '-', lw=1, color=yri_color, label='Model A')
    ax5.plot(steps*generation_time, 1 / 2 / coal_rates_mig[('YRI', 'YRI')],
             '--', lw=1, color=yri_color, label='Model B')
    
    ax5.set_xscale('log')
    ax5.set_ylabel('Inverse coalalescence rate')
    ax5.set_xlim([5000, 1e6])
    ax5.legend(frameon=False, fontsize=6)

    ax6 = plt.subplot2grid((3, 6), (2, 2), colspan=2)
    
    ax6.plot(steps*generation_time, ceu_size, ':', color='gray', label='CEU size')
    ax6.plot(steps*generation_time, 1 / 2 / coal_rates[('CEU', 'CEU')],
             '-', lw=1, color=ceu_color, label='Model A')
    ax6.plot(steps*generation_time, 1 / 2 / coal_rates_mig[('CEU', 'CEU')],
             '--', lw=1, color=ceu_color, label='Model B')
    
    ax6.set_xscale('log')
    ax6.set_xlim([5000, 1e6])
    ax6.legend(frameon=False, fontsize=6)

    ax7 = plt.subplot2grid((3, 6), (2, 4), colspan=2)
    
    ax7.plot(steps*generation_time, chb_size, ':', color='gray', label='CHB size')
    ax7.plot(steps*generation_time, 1 / 2 / coal_rates[('CHB', 'CHB')],
             '-', lw=1, color=chb_color, label='Model A')
    ax7.plot(steps*generation_time, 1 / 2 / coal_rates_mig[('CHB', 'CHB')],
             '--', lw=1, color=chb_color, label='Model B')
    
    ax7.set_xscale('log')
    ax7.set_xlim([5000, 1e6])
    ax7.legend(frameon=False, fontsize=6)

    ax5.set_ylim([0, 20e3])
    ax6.set_ylim([0, 20e3])
    ax7.set_ylim([0, 20e3])
    ax5.set_xlabel('Time in past (years)')
    ax6.set_xlabel('Time in past (years)')
    ax7.set_xlabel('Time in past (years)')

    fig.text(0.025, 0.98, 'A', fontsize=9, ha='center', va='center')
    fig.text(0.5, 0.98, 'B', fontsize=9, ha='center', va='center')
    fig.text(0.025, 0.65, 'C', fontsize=9, ha='center', va='center')
    fig.text(0.525, 0.65, 'D', fontsize=9, ha='center', va='center')
    fig.text(0.025, 0.32, 'E', fontsize=9, ha='center', va='center')
    fig.text(0.35, 0.32, 'F', fontsize=9, ha='center', va='center')
    fig.text(0.675, 0.32, 'G', fontsize=9, ha='center', va='center')

    #fig.tight_layout()
    fig.subplots_adjust(wspace=1.5, hspace=0.4, right=0.975, left=0.1, top=0.975, bottom=.08)
    plt.savefig('../figures/ooa_expected_stats.pdf', dpi=300)
