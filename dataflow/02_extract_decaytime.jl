# Extract decay time of exponential tail of HPGe waveforms using peak files 
# Perform fit (truncated gaussian) to get average decay time, that is later used for pole-zero correction of all waveforms
# Note: Use peak files instead of all waveforms, because they are "good" waveforms and are quality cuts 
using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDataManagement.LDMUtils
using LegendDSP
using LegendDSP: get_fltpars
using LegendHDF5IO
using LegendSpecFits
using RadiationSpectra
using RadiationDetectorDSP
using Measurements: value as mvalue
using PropDicts
using StatsBase, IntervalSets
using Unitful
using TypedTables
using Plots
using Measures

# set data configuration (where to find data; and where to save results)
ENV["LEGEND_DATA_CONFIG"] = "/global/cfs/projectdirs/m2676/data/teststands/lbnl/ppc01/config.json"

# include relevant functions 
include("$(@__DIR__)/../processing_funcs/process_decaytime.jl")
include("$(@__DIR__)/../utils/utils_aux.jl")

# inputs
asic = LegendData(:ppc01)
period = DataPeriod(1)
run = DataRun(1)
channel = ChannelId(1)
category = DataCategory(:cal)

# create decay time pars using default values for pz_config and dsp_config
process_decayime(asic, period, run, category, channel; reprocess = true)

# sanity check: read decay time pars
pars_pz = asic.par.rpars.pz[period, run, channel]

###### Alternative ways to call function: 
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

