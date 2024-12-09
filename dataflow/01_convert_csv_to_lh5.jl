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
category = :cal 
ti = 0.0u"µs"..550.0u"µs" # time interval within waveform is truncated and saved (to reduce file size)
csv_folder = "/Users/lisa/Documents/Workspace/LEGEND/LBL_ASIC/ASIC_data/" * "L1k65n_CSA_DATA_Board_A/Board_A_2nd_Run_Detector/sept27samebut10k/"

# convert csv files to lh5 files
csv_to_lh5(asic, period, run, category, channel, csv_folder; csv_heading = 14, nChannels = 1, nwvfmax = NaN, ti = ti)

# read waveforms as sanity check 
filekeys = search_disk(FileKey, asic.tier[DataTier(:raw), category , period, run])
data = read_ldata(asic, DataTier(:raw), filekeys[1], channel)

Table(data)