
"""
    process_dsp(dsp_config::DSPConfig, data::LegendData, period::DataPeriod, run::DataRun, category::Symbol, channel::ChannelId; reprocess::Bool = false)
- run the DSP processing for all raw files in the given period, run, category and channel.
- based on "simple_dsp" function
- save the results in the jldsp tier
- if reprocess is false, it will skip the files that are already processed
INPUTS:
    - `data::LegendData` LegendData object. You need `"LEGEND_DATA_CONFIG"` to construct this, e.g. `l200 = LegendData(:l200)`
    - `period::DataPeriod` data period, e.g. `DataPeriod(1)`
    - `run::DataRun` data run, e.g. `DataRun(1)`
    - `category::Symbol` data category, e.g. `DataCategory(:cal)`
    - `channel::ChannelId` channel id, e.g. `ChannelId(1)` (depending on your data!)
    - `reprocess::Bool` reprocess the files or not
OUTPUTS:
    - save the DSP results in the jldsp tier
    - print the progress
    - print the completion message
"""
function process_dsp(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; pars_filter::PropDict = PropDict(), reprocess::Bool = false)
    dsp_config = DSPConfig(data.metadata.config.dsp.dsp_config.default)

    filekeys = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])
    dsp_folder =  data.tier[DataTier(:jldsp), category , period, run] * "/"
    τ_pz = mvalue(data.par.rpars.pz[period, run, channel].τ)

    if reprocess == false
        dsp_files = dsp_folder .* string.(filekeys) .* "-tier_jldsp.lh5"
        if sum(isfile.(dsp_files)) == length(filekeys)
            @info "all raw files are already processed. You're already done!"
            println("to read dsp files: use read_ldata(data, :jldsp, category , period, run, channel)")
            return 
        else
            @info "$(sum(isfile.(dsp_files))) raw files are already processed - skip"
        end
        filekeys = filekeys[.!isfile.(dsp_files)]
    end

    # Threads.@threads 
    for f in eachindex(filekeys) 
        filekey = filekeys[f]
        data_raw = read_ldata(data, DataTier(:raw), filekey, channel)
        dsp_par = simple_dsp(Table(data_raw), dsp_config; τ_pz = τ_pz, pars_filter = pars_filter)
       
        # save dsp results
        if !ispath(dsp_folder)
            mkpath(dsp_folder)
        end
        dsp_file = dsp_folder * string(filekey) * "-tier_jldsp.lh5"
        fdsp = lh5open(dsp_file, "w")
        for par in columnnames(dsp_par)
            fdsp["$channel/jldsp/$par"]  = getproperty(dsp_par, par)
        end
        fdsp["$channel/jldsp/eventnumber"] = data_raw.eventnumber
        fdsp["$channel/jldsp/timestamp"] = data_raw.timestamp

        close(fdsp)
        println("dsp processing done for $filekey")
    end
    println("DONE! To read the generated dsp files: read_ldata(data, :jldsp, category , period, run, channel)")
end