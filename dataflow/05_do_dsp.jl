
relPath = relpath(split(@__DIR__, "main")[1], @__DIR__)
using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDataManagement.LDMUtils
using RadiationDetectorDSP
using RadiationDetectorSignals
using ArraysOfArrays
using LegendDSP
using LegendDSP: get_fltpars
using Measurements: value as mvalue
using CSV, DataFrames
using HDF5, LegendHDF5IO
using Dates
using PropDicts
using IntervalSets
using Unitful
using TypedTables
using Statistics
include("$(@__DIR__)/$relPath/src/simple_dsp.jl")
include("$(@__DIR__)/$relPath/processing_funcs/process_dsp.jl")

# inputs 
asic = LegendData(:l1k65n)
period = DataPeriod(1)
run = DataRun(1)
channel = ChannelId(1)
category = DataCategory(:cal)

# do dsp
pars_filter_opt = PropDict() # filter optimization still has to be done! 
process_dsp(asic, period, run, category, channel; pars_filter = pars_filter_opt, reprocess = true)

# read dsp pars
dsp_pars = read_ldata(asic, :jldsp, category, period, run, channel);
Table(dsp_pars)
columnnames(Table(dsp_pars))
