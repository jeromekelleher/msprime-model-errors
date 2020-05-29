# msprime-model-errors

This repo collects together information about two errors in 
[msprime](https://msprime.readthedocs.io/) simulations. The 
[manuscript](https://github.com/jeromekelleher/msprime-model-errors/blob/master/paper.pdf)
describes the problems and background in detail.

In short, there are two distinct problems here:

1. The three population Out-of-Africa model given as an example in the msprime 
   documentation was not an accurate description of the 
   [true model](https://doi.org/10.1371/journal.pgen.1000695). In the most 
   ancient time period, migration was allowed to continue between ancestral
   African and European populations. Fortunately, the difference is a subtle
   one, and the differences in expected diversity measures between the models
   is small. However, this code has been extensively copied --- see below.
   
2. The [simulation pipeline](https://github.com/armartin/ancestry_pipeline) for the 
   analysis for the influential [Martin et al paper](https://doi.org/10.1016/j.ajhg.2017.03.004) 
   contained an error, leading to the simulations being of a substantially 
   different model from what was expected.

See the manuscript for more information and analyis.

## How did this happen?

The error present in the msprime documentation was found as part of the 
quality control process for [stdpopsim](https://stdpopsim.readthedocs.io/en/latest/),
as described in the [preprint](https://www.biorxiv.org/content/10.1101/2019.12.20.885129v2).

See these issues for more details:

- https://github.com/popsim-consortium/stdpopsim/pull/496#issuecomment-620408186
- https://github.com/popsim-consortium/stdpopsim/pull/496
- https://github.com/popsim-consortium/stdpopsim/issues/516

## What do I need to do?

If you have copied incorrect code, you have two basic options to fix it:

### Use stdpopsim

Defining demographic models is hard and error-prone. In an attempt to 
reduce the duplicated effort involved in reimplementing published models
multiple times, the [PopSim Consortium](https://github.com/popsim-consortium)
developed [stdpopsim](https://stdpopsim.readthedocs.io/en/latest/), a 
standard library of population genetic models. The correct version of the 
three population Out-of-Africa model is 
[defined](https://stdpopsim.readthedocs.io/en/latest/catalog.html#sec_catalog_homsap_models_outofafrica_3g09)
and can be run as simply as

```
   $ stdpopsim HomSap -d OutOfAfrica_3G09 -c chr22 10 -o ooa.trees
```

There is also a [Python API](https://stdpopsim.readthedocs.io/en/latest/api.html)
which can plug directly into your existing pipeline and significantly 
simplify your code.

### Fix your model code

The problem with the OOA model is that migration is allowed between two 
ancestral populations until the indefinite past, when we should only
have a single ancestral population. The solution is to add 
`MigrationRateChange` events to ensure that this erroneous migration
isn't happening.

Here is the correct model with a single randomly mating ancestral population:
```python
import math
import msprime
def out_of_africa():
    # First we set out the maximum likelihood values of the various parameters
    # given in Table 1.
    N_A = 7300
    N_AF = 12300
    N_B = 2100
    N_EU0 = 1000
    N_AS0 = 510
    r_EU = 0.004   # 0.4% EU growth
    r_AS = 0.0055  # 0.55% AS growth
    # Migration rates during the various epochs.
    m_AF_B = 25e-5
    m_AF_EU = 3e-5
    m_AF_AS = 1.9e-5
    m_EU_AS = 9.6e-5
    # Times in Table 1 are provided in years, calculated on the assumption
    # of 25 years per generation: we need to convert back into generations.
    generation_time = 25
    T_AF = 220e3 / generation_time
    T_B = 140e3 / generation_time
    T_EU_AS = 21.2e3 / generation_time
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
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
        msprime.MigrationRateChange(time=T_EU_AS, rate=0),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
        msprime.PopulationParametersChange(
            time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
        # Population B merges into YRI at T_B
        msprime.MassMigration(
            time=T_B, source=1, destination=0, proportion=1.0),
        msprime.MigrationRateChange(time=T_B, rate=0),  
        # Size changes to N_A at T_AF
        msprime.PopulationParametersChange(
            time=T_AF, initial_size=N_A, population_id=0)
    ]
    return {
        'population_configurations':population_configurations,
        'migration_matrix':migration_matrix,
        'demographic_events':demographic_events}
```

## Copies of Out-of-Africa example

These GitHub repos have a copy of the faulty Out-of-Africa model that was in the 
msprime documentation. Each link is to a file containing either a direct copy of the 
code, or code that is obviously derived from it.

The list is probably not exhaustive.

- [DomNelson/wf_coalescent](https://github.com/DomNelson/wf_coalescent/blob/842a3f22c075b6499b13f214adfb752b80c4e4a4/scripts/simulate_ooa.py)
- [Ephraim-usc/egrm](https://github.com/Ephraim-usc/egrm/blob/3baf5009aaf1519ebf175074e46494004849bbc7/egrm/simulation.py)
- [OasisYE/MsprimeSimul](https://github.com/OasisYE/MsprimeSimul/blob/77181059bc2d7f6d5cd970a64a56192f26eccc95/Gutenkunst-out-of-Africa.py)
- [PrincetonUniversity/msprime_scripts](https://github.com/PrincetonUniversity/msprime_scripts/blob/892506fa28af98ed80d76a3db558adcbe8cf34e9/src/Demography_Models.py)
- [TishkoffLab/data_simulation](https://github.com/TishkoffLab/data_simulation/blob/105e137881646b5e5ad054de0eecff072d2d8bbd/generate_simulated_phenogeno.py)
- [YingZhou001/POPdemog](https://github.com/YingZhou001/POPdemog/blob/40c78a7a26c93c6755ebd8a74061b424d49042c3/doc/demo1.py)
- [abwolf/msprime_scripts](https://github.com/abwolf/msprime_scripts/blob/f4a383831b3e4156eb9d732eccc4e3c192453709/src/Demography_Models.py)
- [armartin/ancestry_pipeline](https://github.com/armartin/ancestry_pipeline/blob/2e83e68bb5f32858a95046b4048c49899948ab1d/simulate_prs.py)
- [arundurvasula/migration](https://github.com/arundurvasula/migration/blob/017b7b355e4e3b16e199c046568abe129272930a/migration3.py)
- [astheeggeggs/msprime_sim](https://github.com/astheeggeggs/msprime_sim/blob/8ec1945290fcfd2889dbb2a677e21012162fbc89/src/msprime_sim_scenarios.py)
- [awohns/relative-allele-age](https://github.com/awohns/relative-allele-age/blob/7c22e19917207d4e31aebad73ebdb86bacf553df/out_of_africa_fig_transfer/ooa_sim.py)
- [awohns/tsdate_paper](https://github.com/awohns/tsdate_paper/blob/d6e0cee1393f3dc1cf1112a6a7e543a6c2e1a0cb/src/evaluation.py)
- [carjed/primeval](https://github.com/carjed/primeval/blob/9fa2442ad7a28cbf75af83aa73e24007c3d17abb/primeval.py)
- [carjed/topmed_singleton_clusters](https://github.com/carjed/topmed_singleton_clusters/blob/65e5de669af5499703c50aebeeadd50dabc4b96e/scripts/.ipynb_checkpoints/simulate_ext_branches-checkpoint.ipynb)
- [cran/POPdemog](https://github.com/cran/POPdemog/blob/c58939d20c253d2cf18ef30397b4351bbb7ed1bd/inst/doc/popdemog_tutorial.Rmd)
- [dmctong/rv_imp](https://github.com/dmctong/rv_imp/blob/a73fe5074f8ea630e7cd672a1294f24930861fba/s2018-10-25.pipeline1.AFR.py)
- [fbaumdicker/AIMsetfinder](https://github.com/fbaumdicker/AIMsetfinder/blob/17737eb149d3927a705285c76a4c5fc6b4783f68/code/simulate4biogeo.py)
- [isaacovercast/gimmeSAD](https://github.com/isaacovercast/gimmeSAD/blob/1e6dfce63c30997a90051775cff4a5fd4e0ace96/ipython-notebooks/msprime-debugging.ipynb)
- [jiahuanglin/GSoC2019](https://github.com/jiahuanglin/GSoC2019/blob/15659f5c4bb686f127f8e2f354c4d554618e2418/tut/out_of_africa.ipynb)
- [jshleap/Cotagging_playground](https://github.com/jshleap/Cotagging_playground/blob/7700af78408a38a114149a17b1f134d7481c5682/Simulate_PRS.py)
- [mathii/spectrum](https://github.com/mathii/spectrum/blob/955c9d56c435227a58a149f2bb976b479038dfbc/simulate_demography.py)
- [mccoy-lab/sim_introgression](https://github.com/mccoy-lab/sim_introgression/blob/79165bdf59f35f3be1226cb1140fcc9308c13064/sim_introgression.py)
- [mcveanlab/treeseq-inference](https://github.com/mcveanlab/treeseq-inference/blob/697faec29b61b6ff46a1dc8bfc32f8c32e0ba56a/src/evaluation.py)
- [mcveanlab/tskit-workshop](https://github.com/mcveanlab/tskit-workshop/blob/4785d230526083710ec091f3128a302e868b3d6d/ts_workshop_part2.ipynb)
- [messDiv/MESS](https://github.com/messDiv/MESS/blob/8e96c5f68dbf1ccd34ff4b7a6dd18e08f419474e/jupyter-notebooks/_arch/msprime-stuff.ipynb)
- [nikbaya/risk_gradients](https://github.com/nikbaya/risk_gradients/blob/cf1ad95bc8249be0275034c357193bbf46c8d73f/python/msprime_prs.py)
- [pkalbers/ScriptCollection](https://github.com/pkalbers/ScriptCollection/blob/148943ef77afe39bcd4be713e10bcee92d34ce55/demography/simulate.py)
- [pmckenz1/quartet_proj](https://github.com/pmckenz1/quartet_proj/blob/8a29b1a3d4cbce47d3d9bfb1d9281b81f590fdb8/simulate_introgression.ipynb)
- [popgengent/pipeline](https://github.com/popgengent/pipeline/blob/735cdcc5cb240a4bb3f8911fc5b65bec4cc09003/simulate_prs.py)
- [slowkoni/out-of-africa-model](https://github.com/slowkoni/out-of-africa-model/blob/587e10fac40ea29e35659409698c67eadb75e8a8/msprime-out-of-africa-3-pops.py)
- [tszfungc/scripts](https://github.com/tszfungc/scripts/blob/af30ebc6a5862550a0ac715f4ac6bc38d8e0c16c/simulation/sim_demography.py)


## Copies of Martin et al pipeline

The analysis for the [Martin et al paper](https://doi.org/10.1016/j.ajhg.2017.03.004) is 
define in [armartin/ancestry_pipeline](https://github.com/armartin/ancestry_pipeline/blob/2e83e68bb5f32858a95046b4048c49899948ab1d/simulate_prs.py).
We found the following repos that have code derived from it:

- [nikbaya/risk_gradients](https://github.com/nikbaya/risk_gradients/blob/cf1ad95bc8249be0275034c357193bbf46c8d73f/python/msprime_prs.py)
- [jshleap/Cotagging_playground](https://github.com/jshleap/Cotagging_playground/blob/7700af78408a38a114149a17b1f134d7481c5682/Simulate_PRS.py)
- [popgengent/pipeline](https://github.com/popgengent/pipeline/blob/735cdcc5cb240a4bb3f8911fc5b65bec4cc09003/simulate_prs.py)


## Publications affected

By searching for the GitHub repository URLs above, we were able to identify 
a number of papers that may be affected by the erroneous model.

- [Dating genomic variants and shared ancestry in population-scale sequencing
  data](https://doi.org/10.1371/journal.pbio.3000586). The OOA model was 
  used as an example of a complicated demography when evaluating the 
  accuracy of allele age estimates. The precise details of the 
  model are not important and it is highly unlikely the incorrect 
  model specification has any impact.
- [Inferring whole-genome histories in large population datasets](https://doi.org/10.1038/s41588-019-0483-y).
  The OOA model was used here as an example of a more complex demographic
  history, used to test ancestry inference methods. The specifics of the
  demography is not important, and the model does not affect the conclusions
  of the paper in any way.
- [Population genetic simulation study of power in association testing across
  genetic architectures and study designs](https://doi.org/10.1002/gepi.22264)
  The authors use an implementation of the [Tennessen
model](https://stdpopsim.readthedocs.io/en/latest/catalog.html#sec_catalog_homsap_models_outofafrica_2t12)
  that is based on the incorrect msprime OOA example. It would 
  appear that this implemented model also does not switch off 
  migration in the most ancient time period. However, the method 
  is not concerned with detecting detailed population structure, and 
  so the details of the model used are unlikely to be significant.
- [An integrated model of population genetics and community ecology](https://doi.org/10.1111/jbi.13541)
  The OOA model appears in the GitHub repository associated with this paper
  ([isaacovercast/gimmeSAD](https://github.com/isaacovercast/gimmeSAD)), but it
  appears to only have been used as a temporary debugging example.
- [POPdemog: Visualizing Population Demographic History From Simulation
  Scripts](https://doi.org/10.1093/bioinformatics/bty184) POPDemog is a method 
  for visualising demographic histories as described by a number of population
  genetic tools. The OOA example is included as an example of how they 
  can convert msprime input into ms compatible demography descriptions,
  which they then process.
- [How to choose sets of ancestry informative markers: A supervised feature
  selection approach](https://doi.org/10.1016/j.fsigen.2020.102259) In this
  paper the OOA model was used to evaluate a new method for choosing ancestry
  informative markers. Given the very subtle effect of the incorrect
  model on demography (and the fact the method was evaluated using other simulations and real data), 
  it seems unlikely that the model specified will have any effect on the conclusions of the paper.
