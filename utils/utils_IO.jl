
using HDF5, LegendHDF5IO
using CSV, DataFrames
using LegendDataManagement
using RadiationDetectorSignals
using RadiationDetectorDSP
using Unitful
using IntervalSets
using Dates
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
    if _is_valid_datestr(timestr_file; timezone = timezone)
        time_iso = timestr_file[1:8]*"T"*timestr_file[9:end]*timezone
        time_unix = datetime2unix( DateTime(time_iso, "yyyymmddTHHMMSS$(timezone)"))
    else
        @warn "csv files should have a iso-timekey in the filename. If not, current time is used"
        current_datetime = now()
        time_iso = Dates.format(current_datetime, "yyyymmddTHHMMSS$(timezone)")
        time_unix = round(Int64, Dates.datetime2unix(current_datetime))
    end
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

    @debug "Reading $(length(files)) files from $csv_folder"
    fnames = [joinpath(csv_folder, file) for file in files]
    data = [CSV.read(fname, DataFrame; delim=',', header = heading) for fname in fnames]

    # Get MetaData and timestep
    MetaData, timestep = read_csv_metadata(fnames[1], heading = heading, nChannels = nChannels, timezone = timezone)

   # SOME BASIC DATA CLEANING:
    good_wvf = [any(map(x-> occursin(MetaData.Ch1.Channel,x), names(df))) .&&  # good files need to have column named with a channel 
            !any(ismissing.(df[!,  MetaData.Ch1.Channel])) .&& # no missing data points
            all(isfinite.(df[!, MetaData.Ch1.Channel])) for df in data] # no infinite or nan data points
    if any(.!good_wvf)
        # @info "Moving $(sum(MvIdx)) waveforms to ignore folder because of: missing data points, infinite or nan data points, missing channel"
        # dst_folder = csv_folder * "ignore/"
        # if !ispath(dst_folder)
        #     mkpath(dst_folder)
        # end
        # mv.(fnames[MvIdx], [joinpath(dst_folder, file) for file in files[MvIdx]])
        @info "Ignore $(sum(.!good_wvf)) waveforms because of: missing data points, infinite or nan data points, missing channel"
        fnames = fnames[good_wvf]
        data = data[good_wvf]
    end

    # add timestamps for all files (not just first)
    timestr_file = [split.(split.(fn,".csv")[1],"_")[end][1:end-3] for fn in fnames]
    if _is_valid_datestr(timestr_file[1]; timezone = timezone)
        time_iso = [ts[1:8]*"T"*ts[9:end]*timezone for ts in timestr_file] 
        time_unix = datetime2unix.(DateTime.(time_iso, "yyyymmddTHHMMSS$(timezone)"))
    else 
        @warn "csv files should have a iso-timekey in the filename. If not, current time is used"
        current_datetime = fill(now(), length(fnames)) .+ Second.(1:length(fnames))
        time_iso = Dates.format.(current_datetime, "yyyymmddTHHMMSS$(timezone)")
        time_unix = round.(Int64, Dates.datetime2unix.(current_datetime))
    end
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
   
    if nChannels == 1
        return wvfs_ch1, MetaData, good_wvf
    elseif nChannels ==2 
        wvfs_ch2 = ArrayOfRDWaveforms([RDWaveform(time, wvfs) for (time, wvfs) in zip(times, data_ch2)])
        return wvfs_ch1, wvfs_ch2, MetaData, good_wvf
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
        ti::ClosedInterval{<:Quantity} = 0.0u"µs".. 5000.0u"µs", wpf::Int = 1000)
    files = readdir(csv_folder)
    filter!(x -> occursin(".csv", x), files)
    nwvf_total = length(files)
    nwvf_total = ifelse(nwvfmax < nwvf_total, nwvfmax, nwvf_total)
    nfiles = ceil(Int, nwvf_total/wpf) # number of created hdf5 files, based on selected number of waveforms per hdf5 files
    idx_start = [wpf*(i-1)+1 for i = 1:nfiles]
    idx_stop = [wpf*i for i = 1:nfiles]
    idx_stop[end] = nwvf_total
    eventnumber = collect(1:nwvf_total)
    h5folder = data.tier[DataTier(:raw), category, period, run] * "/"
    if !ispath(h5folder)
        mkpath(h5folder)
        @info "created folder: $h5folder"
    end

    # create hdf5 file and write
    #Threads.@threads 
    for i = 1:nfiles
        # read csv files in csv_folder 
        if nChannels == 1
            wvfs_ch1, MetaData, good_wvf = read_folder_csv_oscilloscope(csv_folder; heading = csv_heading, nChannels = nChannels, nwvfmax = collect(idx_start[i]:idx_stop[i])) # make sure you can read all waveforms 
        elseif nChannels == 2
            wvfs_ch1, wvfs_ch2, MetaData, good_wvf = read_folder_csv_oscilloscope(csv_folder; heading = csv_heading, nChannels = nChannels, nwvfmax = collect(idx_start[i]:idx_stop[i])) # make sure you can read all waveforms 
        end 

        # truncate waveforms 
        if rightendpoint(ti) < wvfs_ch1[1].time[end]
            @info "Truncate waveforms to time interval: $ti (original was: $(wvfs_ch1[1].time[1])..$(wvfs_ch1[1].time[end]))"
            uflt_trunc = TruncateFilter(ti)
            wvfs_ch1 = ArrayOfRDWaveforms(uflt_trunc.(wvfs_ch1))
            if nChannels == 2
                wvfs_ch2 = ArrayOfRDWaveforms(uflt_trunc.(wvfs_ch2))
            end
        end
        filekey = string(FileKey(data.name, period, run, category, Timestamp(MetaData.timestamp[1])))
        h5name = h5folder * filekey * "-tier_raw.lh5"
      
        fid = lh5open(h5name, "w")
        fid["$channel/raw/waveform"]  = wvfs_ch1
        fid["$channel/raw/daqenergy"] = maximum.(wvfs_ch1.signal) .- minimum.(wvfs_ch1.signal) #DAQ energy not available in oscilloscope, approx with difference between max and min, needed for compatibility with LEGEND functions
        fid["$channel/raw/eventnumber"]  = eventnumber[idx_start[i]:idx_stop[i]][good_wvf]
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
        @info "saved $(length(wvfs_ch1)) waveforms in .lh5 files with filekey: $filekey , folder: $h5folder"
        close(fid)
    end
end

"""
    _is_valid_datestr(timestr_file::String; timezone::String = "PT")
check if the time string in the filename is a valid date string
"""
function _is_valid_datestr(timestr_file::AbstractString; timezone::String = "PT")
    if length(timestr_file) < 14 
        return false
    else
        time_iso = timestr_file[1:8]*"T"*timestr_file[9:end]*timezone
        try
            time_unix = datetime2unix( DateTime(time_iso, "yyyymmddTHHMMSS$(timezone)"))
            return true
        catch e
            time_unix = datetime2unix( DateTime(time_iso, "yyyymmddTHHMMSS$(timezone)"))
            return false
        end 
    end
end 