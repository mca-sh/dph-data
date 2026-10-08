This repository contains all data and MATLAB scripts produced for the manuscript **Using Acyclic Discrete Phase-type Distributions to Model Kinetic Heterogeneity in Ensemble Dwell time histograms**.

Folders are structured as follows:

```
dph-data/
|-- dataset<D><S><N>/                            # Numerical experiments
|   |-- 01-mladph/                               # Output of ML-ADPH
|   |   |-- <prefix>_simprm/
|   |       |-- <prefix>_1_mldphres.mat          # Analysis results file
|   |       |-- ...
|   |       \-- <prefix>_10_mldphres.mat
|   |   |-- dataset<D><S><N>_linbin_<F>.png      # Fit distribution plots
|   |   \-- ...
|   |-- 02-iem-crp/                              # Output of iEMM
|   |   \-- <prefix>_simprm/
|   |       |-- <prefix>_1_mldphres.mat          # Analysis results file
|   |       |-- ...
|   |       \-- <prefix>_10_mldphres.mat
|   |   |-- dataset<D><S><N>_linbin_<F>.png      # Fit distribution plots
|   |   \-- ...
|   |-- dataset<D><S><N>_linbin_<F>.png          # Ground truth distribution plots for <N>=0
|   |-- ...
|   |-- <prefix>_simprm.mat                      # Simulation parameters file
|   |-- <prefix>_1_simres.mat                    # Simulated data file
|   |-- ...
|   \-- <prefix>_10_simres.mat
|-- EBS-IBS/                                     # EBS-IBS experiment
|-- scripts/                                     # MATLAB scripts for analysis pipeline
|   |-- PHtest_analysisRoutine.m                 # Entry point of the analysis pipeline
|   \-- ...
\-- readme.txt
```

Conventions for directory and file naming are as follows:

### DATASET FOLDERS: `dataset<D><S><N>/`
   Represents a distinct simulation dataset configuration.
   - `<D>` : Aggregate size.
   - `<S>` : Reaction scheme:
   		   - 0 = all non-redundant (or,  for `<D>`=4, a subset of 100) reaction schemes
   		   - 1 = fully coupled reaction scheme, with $w_{ba}=0$ to $1$ and $\frac{tau_{b2}}{tau_{b1}}=1$ to $10$
   - `<N>` : Sample size index:
           - 0 = 500 samples
           - 1 = 5,000 samples
           - 2 = 50,000 samples

### SIMULATION PARAMETERS & DATA REPLICATES (in `dataset<D><S><N>/`):
   - `dataset<D><S><N>_linbin_<F>.png`:
       Plots depicting the ground truth distributions for each simulation parameter set.
       - `<F>` : Image index.
   - `<prefix>_simprm.mat`:
       Simulation parameter file. Each defines a distinct `<prefix>`.
   - `<prefix>_<R>_simres.mat`:
       Simulated dwell time set for replicate `<R>`, generated from `<prefix>_simprm.mat`.
       - `<R>` : Replicate index from 1 to 10.

### METHOD OUTPUTS (in `01-mladph/` and `02-iem-crp/`):
   Subdirectories corresponding to the two analysis methods. Both follow 
   the exact same internal structure:
   - `<prefix>/` :
       Folder containing the results for the corresponding parameter set.
   - `<prefix>_<R>_mldphres.mat`:
       Analysis results for replicate `<R>`.
   - `dataset<D><S><N>_linbin_<F>.png`:
       Plots depicting the fitted distributions for all replicates of each simulation parameter set.
