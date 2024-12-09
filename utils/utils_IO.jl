using HDF5
using CSV, DataFrames
using Unitful
using RadiationDetectorSignals

## CSV conversion functions 
"""
    read_csv_metadata(filepath::String; heading::Int = 17, nChannels::Int = 2, timezone = "PT")
read metadata from a csv file (as taken from oscilloscope)
input:
- filepath::String: file name
- heading::Int: number of lines (beginning at top) in csv file contain metadata 
- nChannels::Int: number of channels in csv file
- timezone::String: timezone to which the timekey in filename refers (standard is PT = Pacific Time --> Berkeley)
"""
function read_csv_metadata(filepath::String; heading::Int = 17, nChannels::Int = 2, timezone = "PT")
    # get time from filename. 
    timestr_file = split(split(filepath,".csv")[1],"_")[end][1:end-3] 
    time_iso = timestr_file[1:8]*"T"*timestr_file[9:end]*timezone
    time_unix = datetime2unix( DateTime(time_iso, "yyyymmddTHHMMSS$(timezone)"))
    MetaData = (timestamp = time_unix, time_iso = time_iso)

    # read header from csv file 
    header_fileds = CSV.read(filepath, DataFrame; delim=',', limit = heading-1, header = false)
    MetaData1_lbls = replace.(collect(header_fileds.Column1), " " => "")
    MetaData1_vals = collect(header_fileds.Column2)
    IdxKeep = map(x -> !ismissing(x), MetaData1_vals)
    MetaData1_nt = NamedTuple{Tuple(Symbol.(MetaData1_lbls[IdxKeep]))}(Tuple(MetaData1_vals[IdxKeep]))
    for (k, val) in pairs(MetaData1_nt)
        try
            MetaData1_nt = merge(MetaData1_nt, (k => parse(Float64, val), ))
        catch
            continue
        end
    end
    if !haskey(MetaData1_nt, :Channel)
        MetaData1_nt = merge(MetaData1_nt, (Channel = MetaData1_nt.TIME, ))
    end 
    MetaData = merge(MetaData, (Ch1 = MetaData1_nt, )) 
    timestep = uconvert(u"µs", MetaData.Ch1.SampleInterval .* uparse(MetaData.Ch1.HorizontalUnits))
  
    if nChannels == 2
        MetaData2_lbls = collect(header_fileds.Column4)
        MetaData2_vals = collect(header_fileds.Column5)
        IdxKeep = map(x -> !ismissing(x), MetaData2_vals)
        MetaData2_nt = NamedTuple{Tuple(Symbol.(MetaData2_lbls[IdxKeep]))}(Tuple(MetaData2_vals[IdxKeep]))
        for (k, val) in pairs(MetaData2_nt)
            try
                MetaData2_nt = merge(MetaData2_nt, (k => parse(Float64, val), ))
                @info k
            catch
                continue
            end
        end
        if !haskey(MetaData2_nt, :Channel)
            MetaData2_nt = merge(MetaData2_nt, (Channel = MetaData2_nt.TIME, ))
        end 
       MetaData = merge(MetaData, (Ch2 = MetaData2_nt, )) 
    end 
    return MetaData, timestep 
end

"""
    read_folder_csv_oscilloscope(csv_folder::String; heading::Int = 17, nwvfmax::Union{Int, Float64, Vector{Int64}} = NaN, nChannels::Int = 2)
read folder with csv files from oscilloscope
input:
- csv_folder::String: absolute csv folder path
- heading::Int: number of lines to skip
- nwvfmax::Union{Int, Float64, Vector{Int64}}: number of waveforms to read OR vector of waveforms indices to read
- nChannels::Int: number of channels in csv files to read (1 or 2)
"""
function read_folder_csv_oscilloscope(csv_folder::String; heading::Int = 17, nwvfmax::Union{Int, Float64, Vector{Int64}} = NaN, nChannels::Int = 2)
    timezone = "PT"
    files = readdir(csv_folder)
    filter!(x -> occursin(".csv", x), files)
    if length(nwvfmax) == 1 && !isnan.(nwvfmax)
        files = files[1:nwvfmax]
    elseif length(nwvfmax) > 1
        files = files[nwvfmax]
    end

    fnames = [joinpath(csv_folder, file) for file in files]
    data = [CSV.read(fname, DataFrame; delim=',', header = heading) for fname in fnames]

    # Get MetaData and timestep
    MetaData, timestep = read_csv_metadata(fnames[1], heading = heading, nChannels = nChannels, timezone = timezone)

   # SOME BASIC DATA CLEANING:
    MvIdx = .![any(map(x-> occursin(MetaData.Ch1.Channel,x), names(df))) .&&  # good files need to have column named with a channel 
            !any(ismissing.(df[!,  MetaData.Ch1.Channel])) .&& # no missing data points
            all(isfinite.(df[!, MetaData.Ch1.Channel])) for df in data] # no infinite or nan data points
    if any(MvIdx)
        @info "Moving $(sum(MvIdx)) waveforms to ignore folder because of: missing data points, infinite or nan data points, missing channel"
        dst_folder = csv_folder * "ignore/"
        if !ispath(dst_folder)
        mkpath(dst_folder)
        end
        mv.(fnames[MvIdx], [joinpath(dst_folder, file) for file in files[MvIdx]])
        fnames = fnames[.!MvIdx]
        data = data[.!MvIdx]
    end

    # add timestamps for all files (not just first)
    timestr_file = [split.(split.(fn,".csv")[1],"_")[end][1:end-3] for fn in fnames]
    time_iso = [ts[1:8]*"T"*ts[9:end]*timezone for ts in timestr_file] 
    time_unix = datetime2unix.(DateTime.(time_iso, "yyyymmddTHHMMSS$(timezone)"))
    MetaData = merge(MetaData, (timestamp = time_unix, time_iso = time_iso))

    # assign data to the right channel
    data_ch1 = [df[!,  MetaData.Ch1.Channel] for df in data if any(map(x-> occursin(MetaData.Ch1.Channel,x), names(df)))]

    # make sure than vertical unit is always in Volts!
    if uparse(MetaData.Ch1.VerticalUnits) != u"V"
        data_ch1 = [ustrip.(uconvert.(u"V", wvf .* uparse(uparse(MetaData.Ch1.VerticalUnits)))) for wvf in data_ch1]
    end 

    if nChannels == 2
        data_ch2 = [df[!, MetaData.Ch2.Channel] for df in data if any(map(x-> occursin(MetaData.Ch2.Channel,x), names(df)))]
        if uparse(MetaData.Ch2.VerticalUnits) != u"V"
            data_ch2 = [ustrip.(uconvert.(u"V", wvf .* uparse(uparse(MetaData.Ch2.VerticalUnits)))) for wvf in data_ch2]
        end 
    end 

    times = fill(0u"µs":timestep:(length(data_ch1[1]) - 1)*timestep, length(data_ch1[1]))
    wvfs_ch1 = ArrayOfRDWaveforms([RDWaveform(time, wvfs) for (time, wvfs) in zip(times, data_ch1)])
    @info "Reading $(length(wvfs_ch1)) waveforms from csv files"
    
    if nChannels == 1
        return wvfs_ch1,  MetaData
    elseif nChannels ==2 
        wvfs_ch2 = ArrayOfRDWaveforms([RDWaveform(time, wvfs) for (time, wvfs) in zip(times, data_ch2)])
        return wvfs_ch1, wvfs_ch2, MetaData
    end 
end

"""
csv_to_lh5(data::LegendData, period::DataPeriod, run::DataRun, category::DataCategoryLike, channel::ChannelId, csv_folder::String; heading::Int = 17, nwvfmax::Union{Int, Float64, Vector{Int64}} = NaN, nChannels::Int = 2, 
        ti::ClosedInterval{<:Quantity} = 0.0u"µs".. 550.0u"µs")
- converts csv files (e.g. from oscilloscope) to lh5 files
- format of csv file matches the one from an oscilloscope...might be different for other systems
- saves the lh5 files in the raw tier defined in "LEGEND_DATA_CONFIG"
INPUTS:
- `data::LegendData` LegendData object. You need `"LEGEND_DATA_CONFIG"` to construct this. this will define later where `.lh5` files are saved
- `period::DataPeriod` data period that you want to assign your data to 
- `run::DataRun` data run that you want to assign your data to 
- `category::DataCategoryLike` data category that you want to assign your data to
- `channel::ChannelId` channel id that you want to assign your data to
- `csv_folder::String` folder where the csv files are located
- `csv_heading::Int` number of lines to skip in csv file
- `nwvfmax::Union{Int, Float64, Vector{Int64}}` number of waveforms to read OR vector of waveforms indices to read
- `nChannels::Int` number of channels in csv files to read (supports only 1 or 2)
- `ti::ClosedInterval{<:Quantity}` time interval to truncate the waveforms to
- `wpf::Int` waveforms per files --> number of waveforms to write per `.lh5` file
"""
function csv_to_lh5(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId, csv_folder::String; csv_heading::Int = 17, nwvfmax::Union{Int, Float64, Vector{Int64}} = NaN, nChannels::Int = 2, 
        ti::ClosedInterval{<:Quantity} = 0.0u"µs".. 550.0u"µs", wpf::Int = 1000)
    files = readdir(csv_folder)
    filter!(x -> occursin(".csv", x), files)
    nwvf_total = length(files)
    nwvf_total = ifelse(nwvfmax < nwvf_total, nwvfmax, nwvf_total)
    nfiles = ceil(Int, nwvf_total/wpf)
    idx_start = [wpf*(i-1)+1 for i = 1:nfiles]
    idx_stop = [wpf*i for i = 1:nfiles]
    idx_stop[end] = nwvf_total
    h5folder = data.tier[DataTier(:raw), category, period, run] * "/"
    eventnumber = collect(1:nwvf_total)

    # create hdf5 file and write
    Threads.@threads for i = 1:nfiles
        # read csv files in csv_folder 
        if nChannels == 1
            wvfs_ch1, MetaData = read_folder_csv_oscilloscope(csv_folder; heading = csv_heading, nChannels = nChannels, nwvfmax = collect(idx_start[i]:idx_stop[i])) # make sure you can read all waveforms 
        elseif nChannels == 2
            wvfs_ch1, wvfs_ch2, MetaData = read_folder_csv_oscilloscope(csv_folder; heading = csv_heading, nChannels = nChannels, nwvfmax = collect(idx_start[i]:idx_stop[i])) # make sure you can read all waveforms 
        end 

        # truncate waveforms 
        if rightendpoint(ti) < wvfs_ch1[1].time[end]
            uflt_trunc = TruncateFilter(ti)
            wvfs_ch1 = ArrayOfRDWaveforms(uflt_trunc.(wvfs_ch1))
            if nChannels == 2
                wvfs_ch2 = ArrayOfRDWaveforms(uflt_trunc.(wvfs_ch2))
            end
        end
        h5name = h5folder * string(FileKey(:l1k65n, period, run, category, Timestamp(MetaData.timestamp[1]))) * "-tier_raw.lh5"
       
        fid = lh5open(h5name, "w")
        fid["$channel/raw/waveform"]  = wvfs_ch1
        fid["$channel/raw/daqenergy"] = maximum.(wvfs_ch1.signal) .- minimum.(wvfs_ch1.signal) #DAQ energy not available in oscilloscope, approx with difference between max and min, needed for compatibility with LEGEND functions
        fid["$channel/raw/eventnumber"]  = eventnumber[idx_start[i]:idx_stop[i]] 
        fid["$channel/raw/baseline"] = fill(NaN, length(wvfs_ch1)) # not available in csv files, but needed for compatibility with LEGEND functions
        for (key, value) in pairs(MetaData)
            if key == :Ch1
                fid["$channel/daq_metadata"] = value
            else
                fid["$channel/raw/$key"] = value
            end
        end
        if nChannels == 2
            fid["$(ChannelId(99))/raw/waveform"]  = wvfs_ch2 # hardcoded for now...
        end
        @info "saved to file: $h5name"
        close(fid)
    end
end