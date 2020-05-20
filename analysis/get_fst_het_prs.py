# This script runs some number of simulations under the correct
# and incorrect demographic model parameterized in Gravel et al (2011).
# We use branch statistics to compute heterozygosity in each population
# and pairwise Fst between pops.

# Usage: python get_fst_het_prs.py num_reps num_samples_per_population

import models
import msprime
import tskit
import sys


if __name__ == "__main__":
    reps, nsam = int(sys.argv[1]), int(sys.argv[2])

    pc, mm, de = models.out_of_africa_gravel()
    for i in range(3):
        pc[i].sample_size = nsam

    # run without demographic events
    H0_inc = 0
    H1_inc = 0
    H2_inc = 0
    F01_inc = 0
    F02_inc = 0
    F12_inc = 0

    for _ in range(reps):
        ts = msprime.simulate(
            population_configurations=pc,
            migration_matrix=mm,
            recombination_rate=0.01,
            length=1,
        ) 
        H0_inc += ts.diversity(range(nsam), mode='branch') 
        H1_inc += ts.diversity(range(nsam, 2*nsam), mode='branch') 
        H2_inc += ts.diversity(range(2*nsam, 3*nsam), mode='branch')
        F01_inc += ts.Fst((range(nsam), range(nsam, 2*nsam)), mode='branch')
        F02_inc += ts.Fst((range(nsam), range(2*nsam, 3*nsam)), mode='branch')
        F12_inc += ts.Fst((range(nsam, 2*nsam), range(2*nsam, 3*nsam)), mode='branch')

    # run with demographic events
    H0_cor = 0
    H1_cor = 0
    H2_cor = 0
    F01_cor = 0
    F02_cor = 0
    F12_cor = 0

    for _ in range(reps):
        ts = msprime.simulate(
            population_configurations=pc,
            migration_matrix=mm,
            recombination_rate=0.01,
            length=1,
            demographic_events=de
        ) 
        H0_cor += ts.diversity(range(nsam), mode='branch') 
        H1_cor += ts.diversity(range(nsam, 2*nsam), mode='branch') 
        H2_cor += ts.diversity(range(2*nsam, 3*nsam), mode='branch')
        F01_cor += ts.Fst((range(nsam), range(nsam, 2*nsam)), mode='branch')
        F02_cor += ts.Fst((range(nsam), range(2*nsam, 3*nsam)), mode='branch')
        F12_cor += ts.Fst((range(nsam, 2*nsam), range(2*nsam, 3*nsam)), mode='branch')

    print("Average Fst:")
    print("Incorrect model:")
    print(f"YRI-CEU:\t{F01_inc/reps}")
    print(f"YRI-CHB:\t{F02_inc/reps}")
    print(f"CEU-CHB:\t{F12_inc/reps}")
    print("Correct model:")
    print(f"YRI-CEU:\t{F01_cor/reps}")
    print(f"YRI-CHB:\t{F02_cor/reps}")
    print(f"CEU-CHB:\t{F12_cor/reps}")
    print()
    print("Ratio of heterozygosity (incorrect/correct)")
    print(f"YRI:\t{H0_inc/H0_cor}")
    print(f"CEU:\t{H1_inc/H1_cor}")
    print(f"CHB:\t{H2_inc/H2_cor}")
