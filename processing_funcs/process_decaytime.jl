function process_decayime(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId,
            min_τ::Quantity{T}, max_τ::Quantity{T}, nbins::Int, rel_cut_fit::Real, peak::Symbol,
            bl_window::ClosedInterval{<:Unitful.Time{<:T}}, tail_window::ClosedInterval{<:Unitful.Time{<:T}}; reprocess::Bool = false) where T <: Real
    
    # check if decaytime pars already exist
    pz_file = joinpath(mkpath(data_path(data.par.rpars.pz[period])), "$(string(run)).json")
    if isfile(pz_file) && !reprocess
        @info "Decay time pars (pz) already exist for $category period $period - run $run - channel $channel - you're done!"
        return
    end

    @debug "Create pars db"
    mkpath(joinpath(data_path(data.par.rpars.pz), string(period)))
    
    data_peak  = read_ldata((peak), asic, :jlpeaks, category, period, run, channel)
    wvfs = data_peak.waveform
    decay_times = dsp_decay_times(wvfs, bl_window, tail_window)
    # filter!(x -> 0.0u"µs" < x < 10* max_τ, decay_times)

    cuts_τ = cut_single_peak(decay_times, min_τ, max_τ,; n_bins=nbins, relative_cut=rel_cut_fit)
    result, report = fit_single_trunc_gauss(decay_times, cuts_τ)
    p = plot(report, size = (600, 500), legend = :topright, xlabel = "Decay time (µs)", fillcolor = :deepskyblue2, color = :darkorange, dpi = 200)
    plot!(p, xguidefontsize = 16, yguidefontsize = 16, xtickfontsize = 12, ytickfontsize = 12, legendfontsize = 10,
                legendforegroundcolor = :silver)
    plot!(p, xlabel = " ", xtickfontsize = 1, bottom_margin = -6mm, ylims = (0, ylims()[2]+0.1), subplot = 1)
    title!(p, "Decay Time Distribution $peak ($period, $run, $channel)", subplot=1, titlefontsize = 12)
    filekey = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])[1]
    savelfig(savefig, p, asic, filekey, channel, :decay_time)

    @info "Found decay time at $(round(u"µs", result.µ, digits=2)) for channel $channel"
    result_pz = (τ = result.μ, fit = result)
    writelprops(data.par.rpars.pz[period], run, PropDict("$channel" => result_pz))
    @info "Saved pars to disk"
end
# export process_decayime

function process_decayime(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId,
    pz_config::PropDict, bl_window::ClosedInterval{<:Unitful.Time{<:T}}, tail_window::ClosedInterval{<:Unitful.Time{<:T}}; kwargs...) where T <: Real
   
    process_decayime(data, period, run, category, channel, pz_config.min_tau, pz_config.max_tau, pz_config.nbins, pz_config.rel_cut_fit, Symbol(pz_config.peak), bl_window, tail_window; kwargs...)
end 

function process_decayime(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId,
    pz_config::PropDict, dsp_config::DSPConfig; kwargs...)
   
    process_decayime(data, period, run, category, channel, pz_config.min_tau, pz_config.max_tau, pz_config.nbins, pz_config.rel_cut_fit, Symbol(pz_config.peak), dsp_config.bl_window, dsp_config.tail_window; kwargs...)
end 

function process_decayime(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; 
    pz_config::PropDict = data.metadata.config.dsp.dsp_config.pz.default, 
    dsp_config::DSPConfig = DSPConfig(data.metadata.config.dsp.dsp_config.default), kwargs...)
    # use default values for pz_config and dsp_config from metadata 
    process_decayime(data, period, run, category, channel, pz_config.min_tau, pz_config.max_tau, pz_config.nbins, pz_config.rel_cut_fit, Symbol(pz_config.peak), dsp_config.bl_window, dsp_config.tail_window; kwargs...)
end 

