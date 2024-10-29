# Test-and-treat posts to reduce malaria in border regions :world_map: ---- :mosquito: :construction: ---- :mosquito: 

This repo uses a [metapopulation version](https://github.com/mrc-ide/malariasimulation/tree/feat/metapopulation_model) of the [Imperial College Individual Malaria Model](https://github.com/mrc-ide/malariasimulation).


## Objectives
:one: Estimate the effectiveness of border posts on reducing cases of malaria along international borders
  
:two: Identify characteristics of the sites most amenable to implementation of border posts

:three: Determine the level of coverage required to make border posts an effective intervention


## Directory

```
.
├── 01_data                                         # Data files
|   ├── Marshall et al. 2016                          # Travel data from Marshall et al. 2016
|   ├── Marshall et al. 2018                          # Trip and coordinate data from Marshall et al. 2018
|   ├── ssa_demography_2021.csv                       # Demography
├── 02_code                                         # R code files
|   ├── Figures                                       # Code for paper figures 
|   ├── HPC_functions                                 # Functions to run model and process results
|   |   ├── function_mixing_matrices.R                  # Helper-function for mixing matrices
|   |   ├── function_run_metapopulation_casestudies.R   # Helper-function for PfPR calibration and case studies
|   |   ├── function_run_metapopulation_mix_examples.R  # Helper-function for examples
|   |   ├── function_run_metapopulation_sims.R          # Helper-function for main metapopulation runs
|   ├── A_shapefiles_centroids.R                      # Download shapfiles and find admin1 centroids
|   ├── B_population.R                                # Extract population sizes 
|   ├── C_clusters.R                                  # Create metapopulation clusters 
|   ├── D_mixing_matrices.R                           # Create mixing matrices 
|   ├── E_parameterization.R                          # Parameterize clusters 
|   ├── F_modelling.R                                 # Run model
|   ├── G_processing.R                                # Process model results
|   ├── MISC_border_mixing_casestudy.R                # Border mixing case study
|   ├── MISC_calibrate.R                              # Calibrate PfPRs to EIRs 
|   ├── MISC_crossing_coverage_casestudy.R            # Case study: border coverage
|   ├── MISC_microscopy_RDT_PCR.R                     # Transforming PCR to RDT prevalence
|   ├── MISC_netz.R                                   # Exploring net distribution vs. use 
|   ├── MISC_PfPR_casestudy.R                         # Case study: PfPR combos 
|   ├── MISC_refitting_gravity_model.R                # Extracting Marshall et al. travel times
|   ├── MISC_three_unit_model.R                       # Toy model
|   ├── MISC_travel_prob.R                            # DHS overnight travel data
|   ├── MISC_trip_duration_travel_time_fits.R         # Trip duration fits
|   ├── packages_paths.R                              # Libraries 
├── 03_output                                       # Figures and .rds files
├── border_elimination.Rproj                        # R.Studio project file
├── sessionInfo.txt                                 # R version and package versions
└── README.md                                       # Project overview

```

## R Package Versions

gdistance: 1.6.4

geodata: 0.5-9

malariaAtlas: 1.5.1

malariasimulation: 1.6.0; [https://github.com/mrc-ide/malariasimulation@ba4bedc](https://github.com/mrc-ide/malariasimulation/commit/ba4bedc2585275abed21c0d6c7fdaf14b31fa7f3)
