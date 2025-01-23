# hpge-ana
Analysis of high-purity Germanium detector data with Juleana. It is structured as follows: 

- **`/processing_funs/`** 

This folder contains a set of files, following the naming scheme `process_<name-of-process>.jl`. 
Each files contains a function definition. This function automatically performs the proccesing for `<name-of-process>`. Example: the file `process_energy_calibration.jl` contains a generic processing function of the energy calibration. 

- **`/dataflow/`**

This folder contains an example on how a full data analysis chain could look like. The files are numbered and meant to be run in ascending order. 

- **`/utils/`**

Some utility function used for I/O and plotting. 

- **`/scripts/`** 

This directory contains the submodule `LBNL_PP01`. This submodule contains analysis scripts (following the structure in dataflow) of the LBNL teststand data. If you want to add your own test-stand analsyis, you can create and add a submodule here. 

# Data configuration / Accessing Data
This code uses the `LegendDataManagement.jl` package to read and write data (tiers and pars) following the LEGEND-200 logic. 
To use it you have to define the environmental variable `"LEGEND_DATA_CONFIG"`, which points to a `.json` file in which all data paths are set. An example for such a config file can be found on NERSC: `export LEGEND_DATA_CONFIG=""/global/cfs/projectdirs/m2676/data/teststands/lbnl/ppc01/config.json"`. An example on how to use this code can be found in the `dataflow/` folder. 

# Cloning 
To clone this repository with submodlues run: `git clone --recurse-submodules -j8 git@github.com:LisaSchlueter/hpge-ana.git` 

# Questions
For questions or requests contact Lisa Schlueter at `lschlueter@lbl.gov`. 