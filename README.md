# hpge-ana
Analysis of high-purity Germanium detector data with Juleana. 

# Data configuration / Accessing Data
This software uses the `LegendDataManagement.jl` package to read and write data (tiers and pars) following the LEGEND-200 logic. 

You have to define the environmental variable `"LEGEND_DATA_CONFIG"`, which points to a `.json` file in which all data paths are set. An example for such a config file can be found on NERSC: `export LEGEND_DATA_CONFIG=""/global/cfs/projectdirs/m2676/data/teststands/lbnl/ppc01/config.json"`

