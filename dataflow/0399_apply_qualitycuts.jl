relPath = relpath(split(@__DIR__, "main")[1], @__DIR__)
using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDataManagement.LDMUtils
using CSV, DataFrames
using HDF5, LegendHDF5IO
using Dates
using PropDicts
using IntervalSets
using Unitful
using TypedTables
using RadiationDetectorDSP
include("$(@__DIR__)/$relPath/utils/utils_IO.jl")

# inputs
asic = LegendData(:l1k65n)
period = DataPeriod(1)
run = DataRun(1)
channel = ChannelId(1)
cat = :cal 

# do qc 
data = asic
qc_config = data.metadata.config.qc.qc_config.default