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
include("$(@__DIR__)/../processing_funcs/process_hit.jl")

# inputs
reprocess = true
asic = LegendData(:ppc01)
period = DataPeriod(1)
run = DataRun(1)
channel = ChannelId(1)
category = :cal 
e_types = [:e_trap, :e_cusp, :e_zac]

# do processing 
process_hit(asic, period, run, category, channel; reprocess = reprocess, e_types = e_types)

# read hit files and plot -> sanity check 
hit_par = Table(read_ldata(asic, :jlhit, category, period, run))
Plots_theme()
stephist(hit_par.e_trap, fill = true, color = :silver, label = "Before qc", nbins =1000, ylims = [1, 1000], yscale = :log10)
stephist!(hit_par.e_trap[hit_par.qc], fill = true, alpha = 1, color = :orange, label = "After qc", xlabel = "Calibrated energy", ylabel = "Counts", nbins = 1000)

