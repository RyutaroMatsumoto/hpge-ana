relPath = relpath(split(@__DIR__, "main")[1], @__DIR__)
using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDataManagement.LDMUtils
using LegendSpecFits
using LegendDSP
using RadiationSpectra
using RadiationDetectorDSP
using CSV, DataFrames
using HDF5, LegendHDF5IO
using Dates
using PropDicts
using StatsBase, IntervalSets
using Unitful
using TypedTables
using Plots 
using Measures
include("$(@__DIR__)/$relPath/processing_funcs/process_decaytime.jl")

# inputs
asic = LegendData(:l1k65n)
period = DataPeriod(1)
run = DataRun(1)
channel = ChannelId(1)
category = DataCategory(:cal)


# create decay time pars / test different ways to call the function
# 1. using default values for pz_config and dsp_config
process_decayime(asic, period, run, category, channel; reprocess = true)

# 2. using pz_config from metadata
pz_config = asic.metadata.config.dsp.dsp_config.pz.default
dsp_config = DSPConfig(asic.metadata.config.dsp.dsp_config.default)
process_decayime(asic, period, run, category, channel, pz_config, dsp_config; reprocess = true)

#3. using custom values
min_τ, max_τ = pz_config.min_tau, pz_config.max_tau
nbins        = pz_config.nbins
rel_cut_fit  = pz_config.rel_cut_fit
peak = Symbol(pz_config.peak)
bl_window = dsp_config.bl_window
tail_window = dsp_config.tail_window
process_decayime(asic, period, run, category, channel, min_τ, max_τ, nbins, rel_cut_fit, peak, bl_window, tail_window; reprocess = true)

## sanity check: read decay time pars
pars_pz = asic.par.rpars.pz[period, run, channel]