"""
    process_qualitycuts(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; reprocess::Bool = false, qc_config::PropDict = data.metadata.config.qc.qc_config.default)
apply quality cuts based on dsp parameters 
inputs: 
    data: LegendData object
    period: DataPeriod object
    run: DataRun object
    category: Symbol or DataCategory object
    channel: ChannelId object
    reprocess: Bool, default false
    qc_config: PropDict, default data.metadata.config.qc.qc_config.default
"""
function process_qualitycuts(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; 
                            reprocess::Bool = false, qc_config::PropDict = data.metadata.config.qc.qc_config.default)

    # check if quality cut pars already exist
    qc_file = joinpath(mkpath(data_path(data.par.rpars.qc[period])), "$(string(run)).json")
     if isfile(qc_file) && !reprocess
         @info "Quality cuts (qc) file already exist for $category period $period - run $run - channel $channel - you're done!"
         return
     end

    # load dsp parameters 
    dsp_par = Table(read_ldata(data, :jldsp, category, period, run, channel))

    # calculate quality cuts
    _, qc = apply_qc(dsp_par, qc_config)

    # add event number and timestamp
    qc[:timestamp] = dsp_par.timestamp
    qc[:eventnumber] = dsp_par.eventnumber

    # save results 
    result_qc = PropDict(Dict("$channel" => qc))
    writelprops(data.par.rpars.qc[period], run, result_qc)
    @info "Saved pars to disk"

    return qc 
end