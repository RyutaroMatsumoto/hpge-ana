
using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDataManagement.LDMUtils
using LegendHDF5IO
using LegendSpecFits
using LegendDSP
using LegendDSP: get_fltpars
using RadiationDetectorSignals
using RadiationDetectorDSP
using Measurements: value as mvalue
using PropDicts
using StatsBase, IntervalSets
using Unitful
using TypedTables
using Plots
using Measures
using Optim
using BSplineKit
using Printf
# set data configuration (where to find data; and where to save results)
ENV["LEGEND_DATA_CONFIG"] = "/global/cfs/projectdirs/m2676/data/teststands/lbnl/l1k65n/config.json"

# include relevant functions 
include("$(@__DIR__)/../src/filteropt_rt_optimization_blnoise.jl")
include("$(@__DIR__)/../utils/utils_plot.jl")
include("$(@__DIR__)/../utils/utils_aux.jl")
include("$(@__DIR__)/../processing_funcs/process_filteropt.jl")

# inputs 
asic = LegendData(:l1k65n)
period = DataPeriod(1)
run = DataRun(1)
channel = ChannelId(1)
category = DataCategory(:cal)

# do optimization 
process_filteropt(asic, period, run, category, channel; reprocess = true, rt_opt_mode = :pickoff)

# read filter optimization pars
fltopt_pars = asic.par.rpars.fltopt[period,run,channel]
