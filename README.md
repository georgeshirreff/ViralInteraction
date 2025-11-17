This repository contains the code for our paper modelling interactions between influenza and RSV. 

System requirements:

I have run the code most recently on:
R version 4.4.3 (2025-02-28 ucrt) -- "Trophy Case"
Platform: x86_64-w64-mingw32/x64
Windows 11 Professionel

The following R packages are required (and the version numbers are just those I'm currently using, not necessarily required): 
tidyverse_2.0.0
odin_1.2.7
FME_1.3.6.3
stringi_1.8.4
ggh4x_0.3.0
coda_0.19-4.1
data.table_1.17.0
ggstance_0.3.7
ggpubr_0.6.0

Installation:

All packages are listed in the "0 installation.R" script. This will also run through the installation of the odin model. 

Installation should not take longer than a few minutes.

Any issues directly with installation of odin can be addressed at https://mrc-ide.github.io/odin/index.html


Execution:

The R scripts are numbered and to be run in that order: 

The script "1 run_odin_m5_twomethods_jumps_reps_lrompsi.R" is used to run the analyses. The user can choose different proposal methods: 
- "cyc" which is described in the paper, and cycles through the seasons proposing season-specific parameters each time
- "basic" a simple MCMC with the same proposal each iteration
- "mod" an adaptive MCMC using FME::modMCMC

The user can also choose different proposal files such as "interactiononly", or create their own and specify the name in the "jumpfile". 

Each run should produce two output files "..._posterior.csv" and "_bestParams.csv" which are all saved iterations, and the best parameter set, respectively. 

The output folder contains the results of previously conducted statistical inference, and can be used to recreate all the figures in the paper. It should take less than 5 minutes to run one model for 1000 iterations. 

The scripts beginning "2..." are used to create output tables and figures, and will directly read analysis files already created in the output folder. 

There is a data_public folder which contains the two public datasets (census data and ECDC community surveillance data) as well as toy versions of the Valencia multiplex data. This data is not currently publicly available, hence the data folder in the repository is empty. 

The code can be run locally. It can also be run using SLURM to allow parallelisation of different parameter values, for which I provide a SLURM shell script:
odin5.sh

Running on your own data:

The analysis requires single detection and co-detection data, for which the current versions are shown in data_public/toydata_allcatch_coinc_age.tsv and data_public/toydata_allcatch_coinc_age.tsv
You can run on your own data by reproducing this format, where iso_standard is the week number (relative to the start of the year), Inc is the number of cases and CatchmentPop is the size of the catchment population for that age group. 

You can also choose to include an ILI dataset as found in data_public/allcatch_ili_age2.tsv

