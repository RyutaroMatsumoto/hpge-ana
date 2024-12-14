using LegendDataManagement
using LegendDataManagement: readlprops, writelprops
using LegendDataManagement.LDMUtils
using LegendHDF5IO
using PropDicts
using Unitful
using TypedTables

# set data configuration (where to find data; and where to save results)
ENV["LEGEND_DATA_CONFIG"] = "/global/cfs/projectdirs/m2676/data/teststands/lbnl/l1k65n/config.json"

# include relevant functions 
include("$(@__DIR__)/../src/apply_qc.jl")
include("$(@__DIR__)/../processing_funcs/process_qualitycuts.jl")

# inputs
asic = LegendData(:l1k65n)
period = DataPeriod(1)
run = DataRun(1)
channel = ChannelId(1)
category = :cal 

# load qc config (optional) and apply 
qc_config = asic.metadata.config.qc.qc_config.default
process_qualitycuts(asic, period, run, category, channel; reprocess = true, qc_config = qc_config);

# read quality cuts from pars 
qc = asic.par.rpars.qc[period][run][channel]


