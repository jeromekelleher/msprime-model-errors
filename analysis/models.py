# Here we store the demographic models used in the analysis as functions that
# return population configurations, demographic events, and migration matrices.
# These are taken from the msprime documentation, before and after the mistake
# in the model was found. They are copied exactly, except to set the choice
# of printing the demography debugger as an argument.


import msprime
import math
import demography
import networkx as nx
import moments, moments.LD as mold


def out_of_africa_ancient_migration(print_debugger=False):
    """
    The out of africa model that does not turn off migration between demes 0 and 1
    after the merger of the African and Eurasian populations.

    print_debugger: If True, prints msprime demography debugger. Defaults to not print.

    Returns population_configurations, migration_matrix, and demographic_events.
    """
    # First we set out the maximum likelihood values of the various parameters
    # given in Table 1.
    N_A = 7300
    N_B = 2100
    N_AF = 12300
    N_EU0 = 1000
    N_AS0 = 510
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_AF = 220e3 / generation_time
    T_B = 140e3 / generation_time
    T_EU_AS = 21.2e3 / generation_time
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.004
    r_AS = 0.0055
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.
    m_AF_B = 25e-5
    m_AF_EU = 3e-5
    m_AF_AS = 1.9e-5
    m_EU_AS = 9.6e-5
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=0, initial_size=N_AF),
        msprime.PopulationConfiguration(
            sample_size=1, initial_size=N_EU, growth_rate=r_EU),
        msprime.PopulationConfiguration(
            sample_size=1, initial_size=N_AS, growth_rate=r_AS)
    ]
    migration_matrix = [
        [      0, m_AF_EU, m_AF_AS],
        [m_AF_EU,       0, m_EU_AS],
        [m_AF_AS, m_EU_AS,       0],
    ]
    demographic_events = [
        # CEU and CHB merge into B with rate changes at T_EU_AS
        msprime.MassMigration(
            time=T_EU_AS, source=2, destination=1, proportion=1.0),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=0),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
        msprime.PopulationParametersChange(
            time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
        # Population B merges into YRI at T_B
        msprime.MassMigration(
            time=T_B, source=1, destination=0, proportion=1.0),
        # Size changes to N_A at T_AF
        msprime.PopulationParametersChange(
            time=T_AF, initial_size=N_A, population_id=0),
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    if print_debugger is True:
        dd = msprime.DemographyDebugger(
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,
            demographic_events=demographic_events)
        dd.print_history()

    return population_configurations, migration_matrix, demographic_events


def out_of_africa(print_debugger=False):
    """
    The out of africa model with ancient migration turned off (i.e. the correct model).

    print_debugger: If True, prints msprime demography debugger. Defaults to not print.

    Returns population_configurations, migration_matrix, and demographic_events.
    """
    # First we set out the maximum likelihood values of the various parameters
    # given in Table 1.
    N_A = 7300
    N_B = 2100
    N_AF = 12300
    N_EU0 = 1000
    N_AS0 = 510
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_AF = 220e3 / generation_time
    T_B = 140e3 / generation_time
    T_EU_AS = 21.2e3 / generation_time
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.004
    r_AS = 0.0055
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.
    m_AF_B = 25e-5
    m_AF_EU = 3e-5
    m_AF_AS = 1.9e-5
    m_EU_AS = 9.6e-5
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=0, initial_size=N_AF),
        msprime.PopulationConfiguration(
            sample_size=1, initial_size=N_EU, growth_rate=r_EU),
        msprime.PopulationConfiguration(
            sample_size=1, initial_size=N_AS, growth_rate=r_AS)
    ]
    migration_matrix = [
        [      0, m_AF_EU, m_AF_AS],
        [m_AF_EU,       0, m_EU_AS],
        [m_AF_AS, m_EU_AS,       0],
    ]
    demographic_events = [
        # CEU and CHB merge into B with rate changes at T_EU_AS
        msprime.MassMigration(
            time=T_EU_AS, source=2, destination=1, proportion=1.0),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=0),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
        msprime.PopulationParametersChange(
            time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
        # Population B merges into YRI at T_B
        msprime.MassMigration(
            time=T_B, source=1, destination=0, proportion=1.0),
        msprime.MigrationRateChange(
            time=T_B, rate=0),
        # Size changes to N_A at T_AF
        msprime.PopulationParametersChange(
            time=T_AF, initial_size=N_A, population_id=0),
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    if print_debugger is True:
        dd = msprime.DemographyDebugger(
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,
            demographic_events=demographic_events)
        dd.print_history()

    return population_configurations, migration_matrix, demographic_events


def ooa_dg():
    N_A = 7300
    N_B = 2100
    N_AF = 12300
    N_EU0 = 1000
    N_AS0 = 510
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_AF = 220e3 / generation_time
    T_B = 140e3 / generation_time
    T_EU_AS = 21.2e3 / generation_time
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.004
    r_AS = 0.0055
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.
    m_AF_B = 25e-5
    m_AF_EU = 3e-5
    m_AF_AS = 1.9e-5
    m_EU_AS = 9.6e-5
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.

    G = nx.DiGraph()
    G.add_node('root', nu=1, T=0)
    G.add_node('A', nu=N_AF/N_A, T=(T_AF-T_B)/2/N_A)
    G.add_node('B', nu=N_B / N_A, T=(T_B-T_EU_AS)/2/N_A,
               m={'YRI': 2 * N_A * m_AF_B})
    G.add_node('YRI', nu=N_AF / N_A, T=T_B / 2 / N_A,
               m={'B': 2 * N_A * m_AF_B, 'CEU': 2 * N_A * m_AF_EU, 'CHB': 2 * N_A * m_AF_AS})
    G.add_node('CEU', nu0=N_EU0/N_A, nuF=N_EU/N_A, T=T_EU_AS/2/N_A,
               m={'YRI': 2 * N_A * m_AF_EU, 'CHB': 2 * N_A * m_EU_AS})
    G.add_node('CHB', nu0=N_AS0/N_A, nuF=N_AS/N_A, T=T_EU_AS/2/N_A,
               m={'YRI': 2 * N_A * m_AF_AS, 'CEU': 2 * N_A * m_EU_AS})

    G.add_edges_from(
        [('root', 'A'),
         ('A', 'B'),
         ('A', 'YRI'),
         ('B', 'CEU'),
         ('B', 'CHB')
         ]
    )

    dg = demography.DemoGraph(G, Ne=N_A)
    return dg


def ooa_mig_dg():
    N_A = 7300
    N_B = 2100
    N_AF = 12300
    N_EU0 = 1000
    N_AS0 = 510
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_AF = 220e3 / generation_time
    T_B = 140e3 / generation_time
    T_EU_AS = 21.2e3 / generation_time
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.004
    r_AS = 0.0055
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.
    m_AF_B = 25e-5
    m_AF_EU = 3e-5
    m_AF_AS = 1.9e-5
    m_EU_AS = 9.6e-5
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.

    G = nx.DiGraph()
    G.add_node('root', nu=1, T=0)
    # from the root, add two pops with the ancient migration that exists for 40N gens
    # A0 and B0 exchanges migrants with each other, and A and B0 also exchange migrants
    G.add_node('A0', nu=1, T=40,
               m={'B': 2 * N_A * m_AF_B})
    G.add_node('A', nu=N_AF/N_A, T=(T_AF-T_B)/2/N_A,
               m={'B': 2 * N_A * m_AF_B},
               pulse={('B', 1, 1)})
    G.add_node('B', nu=N_B / N_A, T=40 + (T_AF-T_EU_AS)/2/N_A,
               m={'A0': 2 * N_A * m_AF_B, 'A': 2 * N_A * m_AF_B, 'YRI': 2 * N_A * m_AF_B})
    G.add_node('YRI', nu=N_AF / N_A, T=T_B / 2 / N_A,
               m={'B': 2 * N_A * m_AF_B, 'CEU': 2 * N_A * m_AF_EU, 'CHB': 2 * N_A * m_AF_AS})
    G.add_node('CEU', nu0=N_EU0/N_A, nuF=N_EU/N_A, T=T_EU_AS/2/N_A,
               m={'YRI': 2 * N_A * m_AF_EU, 'CHB': 2 * N_A * m_EU_AS})
    G.add_node('CHB', nu0=N_AS0/N_A, nuF=N_AS/N_A, T=T_EU_AS/2/N_A,
               m={'YRI': 2 * N_A * m_AF_AS, 'CEU': 2 * N_A * m_EU_AS})

    G.add_edges_from(
        [('root', 'A0'),
         ('root', 'B'),
         ('A0', 'A'),
         ('A', 'YRI'),
         ('B', 'CEU'),
         ('B', 'CHB')
         ]
    )

    dg = demography.DemoGraph(G, Ne=N_A)
    return dg


def get_sfs_ancient_migration(ns):
    N_A = 7300
    N_B = 2100
    N_AF = 12300
    N_EU0 = 1000
    N_AS0 = 510
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_AF = 220e3 / generation_time
    T_B = 140e3 / generation_time
    T_EU_AS = 21.2e3 / generation_time
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.004
    r_AS = 0.0055
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.
    m_AF_B = 25e-5
    m_AF_EU = 3e-5
    m_AF_AS = 1.9e-5
    m_EU_AS = 9.6e-5
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.

    fs = moments.Demographics1D.snm([4*ns])
    fs = moments.Manips.split_1D_to_2D(fs, 3*ns, ns)
    fs.integrate([1, N_B / N_A], 40,
                 m=[[0, 2 * N_A * m_AF_B], [2 * N_A * m_AF_B, 0]])
    fs.integrate([N_AF / N_A, N_B / N_A], (T_AF - T_B) / 2 / N_A,
                 m=[[0, 2 * N_A * m_AF_B], [2 * N_A * m_AF_B, 0]])
    fs = fs.marginalize([1])
    fs = moments.Manips.split_1D_to_2D(fs, ns, 2*ns)
    fs.integrate([N_AF / N_A, N_B / N_A], (T_B - T_EU_AS) / 2 / N_A,
                 m=[[0, 2 * N_A * m_AF_B], [2 * N_A * m_AF_B, 0]])
    fs = moments.Manips.split_2D_to_3D_2(fs, ns, ns)
    nu_func = lambda t: [N_AF / N_A,
                         N_EU0 / N_A * math.exp(math.log(N_EU / N_EU0) * t/(T_EU_AS / 2 / N_A)),
                         N_AS0 / N_A * math.exp(math.log(N_AS / N_AS0) * t/(T_EU_AS / 2 / N_A))]
    fs.integrate(nu_func, T_EU_AS / 2 / N_A,
                 m=[[0, 2 * N_A * m_AF_EU, 2 * N_A * m_AF_AS],
                    [2 * N_A * m_AF_EU, 0, 2 * N_A * m_EU_AS],
                    [2 * N_A * m_AF_AS, 2 * N_A * m_EU_AS, 0]])

    fs.pop_ids = ['YRI', 'CEU', 'CHB']
    return fs


def get_ld_ancient_migration(rho=None, theta=1):
    N_A = 7300
    N_B = 2100
    N_AF = 12300
    N_EU0 = 1000
    N_AS0 = 510
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_AF = 220e3 / generation_time
    T_B = 140e3 / generation_time
    T_EU_AS = 21.2e3 / generation_time
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.004
    r_AS = 0.0055
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.
    m_AF_B = 25e-5
    m_AF_EU = 3e-5
    m_AF_AS = 1.9e-5
    m_EU_AS = 9.6e-5
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.

    y = mold.Demographics1D.snm(rho=rho, theta=1)
    y = y.split(1)
    y.integrate([1, N_B / N_A], 40,
                m=[[0, 2 * N_A * m_AF_B], [2 * N_A * m_AF_B, 0]],
                rho=rho, theta=1)
    y.integrate([N_AF / N_A, N_B / N_A], (T_AF - T_B) / 2 / N_A,
                m=[[0, 2 * N_A * m_AF_B], [2 * N_A * m_AF_B, 0]],
                rho=rho, theta=1)
    y = y.marginalize([2])
    y = y.split(1)
    y.integrate([N_AF / N_A, N_B / N_A], (T_B - T_EU_AS) / 2 / N_A,
                m=[[0, 2 * N_A * m_AF_B], [2 * N_A * m_AF_B, 0]],
                rho=rho, theta=1)
    y = y.split(2)
    nu_func = lambda t: [N_AF / N_A,
                         N_EU0 / N_A * math.exp(math.log(N_EU / N_EU0) * t/(T_EU_AS / 2 / N_A)),
                         N_AS0 / N_A * math.exp(math.log(N_AS / N_AS0) * t/(T_EU_AS / 2 / N_A))]
    y.integrate(nu_func, T_EU_AS / 2 / N_A,
                m=[[0, 2 * N_A * m_AF_EU, 2 * N_A * m_AF_AS],
                   [2 * N_A * m_AF_EU, 0, 2 * N_A * m_EU_AS],
                   [2 * N_A * m_AF_AS, 2 * N_A * m_EU_AS, 0]],
                rho=rho, theta=1)

    y.pop_ids = ['YRI', 'CEU', 'CHB']
    return y


def out_of_africa_gravel(print_debugger=False):
    """
    """
    # First we set out the maximum likelihood values of the various parameters
    # given in Gravel et al, 2011 Table 2.
    N_A = 7300
    N_B = 1861
    N_AF = 14474
    N_EU0 = 1032
    N_AS0 = 554
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_AF = 148e3 / generation_time
    T_B = 51e3 / generation_time
    T_EU_AS = 23e3 / generation_time
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.0038
    r_AS = 0.0048
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.
    m_AF_B = 15e-5
    m_AF_EU = 2.5e-5
    m_AF_AS = 0.78e-5
    m_EU_AS = 3.11e-5
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=0, initial_size=N_AF),
        msprime.PopulationConfiguration(
            sample_size=1, initial_size=N_EU, growth_rate=r_EU),
        msprime.PopulationConfiguration(
            sample_size=1, initial_size=N_AS, growth_rate=r_AS)
    ]
    migration_matrix = [
        [      0, m_AF_EU, m_AF_AS],
        [m_AF_EU,       0, m_EU_AS],
        [m_AF_AS, m_EU_AS,       0],
    ]
    demographic_events = [
        # CEU and CHB merge into B with rate changes at T_EU_AS
        msprime.MassMigration(
            time=T_EU_AS, source=2, destination=1, proportion=1.0),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=0),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
        msprime.PopulationParametersChange(
            time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
        # Population B merges into YRI at T_B
        msprime.MassMigration(
            time=T_B, source=1, destination=0, proportion=1.0),
        msprime.MigrationRateChange(
            time=T_B, rate=0),
        # Size changes to N_A at T_AF
        msprime.PopulationParametersChange(
            time=T_AF, initial_size=N_A, population_id=0),
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    if print_debugger is True:
        dd = msprime.DemographyDebugger(
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,
            demographic_events=demographic_events)
        dd.print_history()

    return population_configurations, migration_matrix, demographic_events

def martin_model(merge_time=None, print_debugger=False):
    """
    """
    generation_time = 25
    pc, mm, de = out_of_africa_gravel()
    if merge_time is not None:
        de = [
            msprime.MassMigration(
                time=merge_time / generation_time,
                source=2, destination=0, proportion=1.0),
            msprime.MassMigration(
                time=merge_time / generation_time,
                source=1, destination=0, proportion=1.0),
            msprime.MigrationRateChange(
                time=merge_time / generation_time,
                rate=0),
        ]
    else:
        de = []

    tol = 1e-5
    if merge_time is None:
        merge_time = 3e5

    merge_time /= generation_time

    # First we set out the maximum likelihood values of the various parameters
    # given in Gravel et al, 2011 Table 2.
    N_A = 7300
    N_B = 1861
    N_AF = 14474
    N_EU0 = 1032
    N_AS0 = 554
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_AF = 148e3 / generation_time
    T_B = 51e3 / generation_time
    T_EU_AS = 23e3 / generation_time
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.0038
    r_AS = 0.0048
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    
    N_EU0 = N_EU * math.exp(-r_EU * merge_time)
    N_AS0 = N_AS * math.exp(-r_EU * merge_time)
    # Migration rates during the various epochs.
    m_AF_B = 15e-5
    m_AF_EU = 2.5e-5
    m_AF_AS = 0.78e-5
    m_EU_AS = 3.11e-5

    G = nx.DiGraph()
    G.add_node('root', nu=1, T=0)
    G.add_node('YRI', nu=N_AF/N_A, T=merge_time/2/N_A,
               m={'CEU':2*N_A*m_AF_EU, 'CHB':2*N_A*m_AF_AS})
    G.add_node('B', nu=1, T=tol)
    G.add_node('CEU', nu0=N_EU0/N_A, nuF=N_EU/N_A, T=merge_time/2/N_A - tol,
               m={'YRI':2*N_A*m_AF_EU, 'CHB':2*N_A*m_EU_AS})
    G.add_node('CHB', nu0=N_AS0/N_A, nuF=N_AS/N_A, T=merge_time/2/N_A - tol,
               m={'YRI':2*N_A*m_AF_AS, 'CEU':2*N_A*m_EU_AS})
    G.add_edges_from([('root','YRI'), ('root','B'), ('B','CEU'), ('B','CHB')])
    dg = demography.DemoGraph(G, Ne=N_A)
    
    return pc, mm, de, dg
