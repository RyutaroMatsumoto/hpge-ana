
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
using LegendHDF5IO
using PropDicts
using IntervalSets
using Unitful
using TypedTables
using Statistics
using Optim
using BSplineKit

# inputs 
asic = LegendData(:l1k65n)
period = DataPeriod(1)
run = DataRun(1)
channel = ChannelId(1)
category = DataCategory(:cal)

data = asic
include("$(@__DIR__)/$relPath/src/filteropt_rt_blnoise.jl")
# include("$(@__DIR__)/$relPath/processing_funcs/process_filteropt.jl")

# TO DO: fix filter optimization --> include QC before. THIS IS NOT WORKING PROPERLY YET! 
process_filteropt(asic, period, run, category, channel; reprocess = true)

# add QC in peak splitting function ? 