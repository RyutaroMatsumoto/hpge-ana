relPath = relpath(split(@__DIR__, "main")[1], @__DIR__)
using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDataManagement.LDMUtils
using RadiationDetectorDSP
using RadiationDetectorSignals
using LegendDSP
using LegendDSP: get_fltpars
using Measurements: value as mvalue
using CSV, DataFrames
using HDF5, LegendHDF5IO
using Dates
using PropDicts
using StatsBase, IntervalSets
using Unitful
using TypedTables
using RadiationSpectra
using RadiationDetectorDSP
include("$(@__DIR__)/$relPath/processing_funcs/process_peak_split.jl")
include("$(@__DIR__)/$relPath/src/simple_dsp_qc.jl")
include("$(@__DIR__)/$relPath/src/apply_qc.jl")

# inputs
asic = LegendData(:l1k65n)
period = DataPeriod(1)
run = DataRun(1)
channel = ChannelId(1)
category = DataCategory(:cal)
source = :co60

# create peak files 
process_peak_split(asic, period, run, category, channel ; reprocess = true)

# sanity check: open peak files and plot 
peakA = read_ldata(:Co60a, asic, :jlpeaks, category, period, run, channel)
peakB = read_ldata(:Co60b, asic, :jlpeaks, category, period, run, channel)
using Plots
stephist(peakA.e_trap, nbins = 400, label = "Co60a", fill = true, xlabel = "Energy")
stephist!(peakB.e_trap, nbins = 300, label = "Co60b", fill = true)
length(peakA.e_simplecal), length(peakB.e_simplecal)