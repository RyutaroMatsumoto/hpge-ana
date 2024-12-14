
using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDataManagement.LDMUtils
using LegendHDF5IO
using LegendSpecFits
using LegendDSP
using LegendDSP: get_fltpars
# using RadiationDetectorSignals
# using RadiationDetectorDSP
using Measurements: value as mvalue
using PropDicts
using StatsBase, IntervalSets
using Unitful
using TypedTables
using Plots
using Measures
using Printf

# set data configuration (where to find data; and where to save results)
ENV["LEGEND_DATA_CONFIG"] = "/global/cfs/projectdirs/m2676/data/teststands/lbnl/l1k65n/config.json"
# include relevant functions 
include("$(@__DIR__)/../src/simple_dsp.jl")
include("$(@__DIR__)/../processing_funcs/process_dsp.jl")

# inputs 
asic = LegendData(:l1k65n)
period = DataPeriod(1)
run = DataRun(1)
channel = ChannelId(1)
category = DataCategory(:cal)

# debug
data = asic 
# do dsp
pars_filter = data.par.rpars.fltopt[period,run,channel] # filter parameter used in dsp 
process_dsp(asic, period, run, category, channel; reprocess = true)

# read dsp pars
dsp_pars = read_ldata(asic, :jldsp, category, period, run, channel);
Table(dsp_pars)
columnnames(Table(dsp_pars))

