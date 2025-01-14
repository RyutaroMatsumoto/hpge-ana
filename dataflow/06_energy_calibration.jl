using LegendDataManagement
using LegendDataManagement: readlprops, writelprops
using LegendDataManagement.LDMUtils
using LegendSpecFits
using LegendSpecFits: get_friedman_diaconis_bin_width 
using LegendHDF5IO
using RadiationSpectra
using PropDicts
using Unitful
using TypedTables
using Statistics, StatsBase
using IntervalSets
using Plots 
using Unitful, Measures
using Measurements: value as mvalue
# set data configuration (where to find data; and where to save results)
ENV["LEGEND_DATA_CONFIG"] = "/global/cfs/projectdirs/m2676/data/teststands/lbnl/ppc01/config.json"

# include relevant functions 
include("$(@__DIR__)/../utils/utils_aux.jl")
include("$(@__DIR__)/../utils/utils_plot.jl")
include("$(@__DIR__)/../processing_funcs/process_energy_calibration.jl")

# inputs
reprocess = true 
asic = LegendData(:ppc01)
period = DataPeriod(1)
run = DataRun(1)
channel = ChannelId(1)
category = :cal 
e_type = :e_trap
ecal_config = asic.metadata.config.energy.energy_config.default

ecal_config.source
data = asic 
# do calibration 
process_energy_calibration(asic, period, run, category, channel; reprocess = reprocess, ecal_config = ecal_config, e_types = [:e_trap, :e_cusp, :e_zac])

# read calibration parameters
asic.par.rpars.ecal[period, run, channel].e_trap



