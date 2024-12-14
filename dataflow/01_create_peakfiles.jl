
using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDataManagement.LDMUtils
using LegendDSP
using LegendDSP: get_fltpars
using LegendHDF5IO
using RadiationSpectra
using RadiationDetectorDSP
using Measurements: value as mvalue
using PropDicts
using StatsBase, IntervalSets
using Unitful
using TypedTables

# set data configuration (where to find data; and where to save results)
ENV["LEGEND_DATA_CONFIG"] = "/global/cfs/projectdirs/m2676/data/teststands/lbnl/l1k65n/config.json"

# include relevant functions 
include("$(@__DIR__)/../processing_funcs/process_peak_split.jl")
include("$(@__DIR__)/../src/simple_dsp.jl")
include("$(@__DIR__)/../src/apply_qc.jl")

# inputs
asic = LegendData(:l1k65n)
period = DataPeriod(1)
run = DataRun(1)
channel = ChannelId(1)
category = DataCategory(:cal)
source = :co60

# create peak files 
process_peak_split(asic, period, run, category, channel ; source = source, reprocess = true)

# sanity check: open peak files and plot 
peakA = read_ldata(:Co60a, asic, :jlpeaks, category, period, run, channel)
peakB = read_ldata(:Co60b, asic, :jlpeaks, category, period, run, channel)

using Plots
stephist(peakA.e_trap, nbins = 400, label = "Co60a", fill = true, xlabel = "Energy")
stephist!(peakB.e_trap, nbins = 300, label = "Co60b", fill = true)
length(peakA.e_simplecal), length(peakB.e_simplecal)