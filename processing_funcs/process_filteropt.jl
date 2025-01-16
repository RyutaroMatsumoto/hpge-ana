"""
    process_filteropt(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId, dsp_config::DSPConfig, τ_pz::Quantity{T}, peak::Symbol; rt_opt_mode::Symbol = :blnoise, reprocess::Bool = false, filter_types::Vector{Symbol} = [:trap, cusp])
    process_filteropt(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; kwargs...)            

Filter optimization for filter_types
- load waveforms from peakfile, shift baseline and pole-zero 
- optimize rise-time for minimum baseline noise after filtering
- optimize flat-top time for FWHM of peak
- save results to disk
- sanity plots for rise-time and flat-top time optimization
Inputs: 
- `data::LegendData`: LegendData object
- `period::DataPeriod`: data period
- `run::DataRun`: data run
- `category::Union{Symbol, DataCategory}`: data category, e.g. :cal
- `channel::ChannelId`: channel id
Optional:
- `dsp_config::DSPConfig`: DSP configuration object. If not specified will take default from metadata
- `τ_pz::Quantity{T}`: pole-zero decay time.  If not specified will take default from metadata
- `peak::Symbol`: peak to optimize for (needs existing peakfile!). If not specified will take default from metadata
- `rt_opt_mode::Symbol`: mode for rise-time optimization (:blnoise or :pickoff) --> two different strategies for optimization
- `reprocess::Bool`: reprocess the files or not
- `filter_types::Vector{Symbol}`: filter types to optimize for
"""
function process_filteropt(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId, dsp_config::DSPConfig, τ_pz::Quantity{T}, peak::Symbol; 
                rt_opt_mode::Symbol = :bl_noise, reprocess::Bool = false, filter_types::Vector{Symbol} = [:trap, :cusp, :zac]) where T<:Real 
    det = _channel2detector(data, channel)
    @info "Optimize filter for period $period, run $run, channel $channel /det $det - $filter_types"

    # check if decaytime pars already exist
    fltopt_file = joinpath(mkpath(data_path(data.par.rpars.fltopt[period])), "$(string(run)).json")
    if isfile(fltopt_file) && !reprocess
        @info "Filter optimization file already exist for $category period $period - run $run - channel $channel - you're done!"
        return
    end

    # prepare results dict
    mkpath(joinpath(data_path(data.par.rpars.fltopt), string(period)))
    result_filteropt_dict = Dict{Symbol, NamedTuple}()
    @debug "Created path for filter optimization results"

    # load waveforms from peakfile
    data_peak  = read_ldata((peak), data, :jlpeaks, category, period, run, channel)
    wvfs = data_peak.waveform
    @debug "Loaded waveforms for peak $peak"

    function process_filteropt_fltr(filter_type::Symbol)
        @info "Optimize filter $filter_type"
        _, def_ft = get_fltpars(PropDict(), filter_type, dsp_config)
        filekey = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])[1]
        # STEP 1: rise-time optimization --> min. baseline noise after filtering 
        if rt_opt_mode == :bl_noise 
            result_rt, report_rt = filteropt_rt_optimization_blnoise(filter_type, wvfs, dsp_config, τ_pz; ft = def_ft)
            Plots_theme()
            p = scatter(ustrip.(collect(report_rt.enc_grid_rt)), report_rt.enc, size = (600, 400), dpi = 150, 
                        ms = 4, color = :deepskyblue2, markerstrokecolor = :deepskyblue2,
                        xlabel = "Rise time ($(unit(report_rt.rt)))", ylabel = "Noise (a.u.)", 
                        title = "Noise sweep ($filter_type); fixed ft = $(def_ft) \n $period, $run, $channel, $peak peak",
                        label = "Data")
            plot!([ustrip(report_rt.rt), ustrip(report_rt.rt)], [ylims()[1], report_rt.min_enc], color = :red2, linewidth = 2, linestyle = :dot,
                            label = @sprintf("Optimal rise-time = %.1f %s", ustrip(report_rt.rt), unit(report_rt.rt)))
            plot!([xlims()[1], ustrip(report_rt.rt)], [report_rt.min_enc, report_rt.min_enc], color = :red2, linewidth = 2, linestyle = :dot, label = false)
            savelfig(savefig,p , data, filekey, channel,  Symbol("noise_sweep_$(filter_type)_blnoise"))
            @info "Save sanity plot to $(LegendDataManagement.LDMUtils.get_pltfilename(data, filekey, channel, :fltopt_rt))"
        elseif rt_opt_mode == :pickoff
            # for now: do 2 versions or rt...the same using LEGEND DSP function 
            # pro: + compatibility with Juleana; proper fitting of enc instead of rms
            # con: - doesnt work well for smaller number of waveforms. result seems to be more noisy, slower?
            enc_grid = getfield(Main, Symbol("dsp_$(filter_type)_rt_optimization"))(wvfs, dsp_config,  τ_pz; ft=def_ft)
            enc_min, enc_max = _quantile_truncfit(enc_grid; qmin = 0.02, qmax = 0.98)
            e_grid_rt   = getproperty(dsp_config, Symbol("e_grid_rt_$(filter_type)"))
            result_rt, report_rt = fit_enc_sigmas(enc_grid, e_grid_rt, enc_min, enc_max, round(Int,size(enc_grid)[2]/5), 0.1)
            @info "Found optimal rise-time: $(result_rt.rt) at fixed ft = $def_ft"
            p = plot(report_rt)
            title!(p, get_plottitle(filekey, det, "Noise Sweep"; additiional_type=string(filter_type)))
            savelfig(savefig, p, data, filekey, channel, Symbol("noise_sweep_$(filter_type)"))
        end 

        # 2. flat top time optimixation 
        e_grid_ft   = getproperty(dsp_config, Symbol("e_grid_ft_$(filter_type)"))
        e_grid = getfield(Main, Symbol("dsp_$(filter_type)_ft_optimization"))(wvfs, dsp_config, τ_pz, mvalue(result_rt.rt))
        e_min, e_max = _quantile_truncfit(e_grid; qmin = 0.02, qmax = 0.98)
        result_ft, report_ft = fit_fwhm_ft(e_grid, e_grid_ft, result_rt.rt,  e_min, e_max, 0.1; peak = data_peak.gamma_line[1])
        @info "Found optimal flattop-time: $(result_ft.ft) with FWHM $(round(u"keV", result_ft.min_fwhm, digits=2))"

        p = plot(report_ft, legendfontsize = 13, xlabel = "Flat-top time (µs)", ylabel = "FWHM $peak (keV)")
        ymin = ifelse(mvalue(ustrip(result_ft.min_fwhm)) - 0.2 < 0, 0, mvalue(ustrip(result_ft.min_fwhm)) - 1)
        ymax = maximum(mvalue.(ustrip.(report_ft.fwhm))) + 0.2
        ylims!(ymin, ymax)
        # title!("$filter_type" * @sprintf("filter: flattop-time optimization; rt = %.2f %s", ustrip(result_rt.rt), unit(result_rt.rt)) * "\n $period, $run, $channel, $peak peak")
        title!(get_plottitle(filekey, det, "$peak FT Scan"; additiional_type=string(filter_type)))
        savelfig(savefig, p, data, filekey, channel, Symbol("fwhm_ft_scan_$(filter_type)"))

        # return result for this filter type
        merge(result_rt, result_ft)
    end 

    # filter optimization: rise-time, and flat-top times
    for filter_type in filter_types
        result_filteropt_dict[filter_type] =  process_filteropt_fltr(filter_type)
    end
    result = PropDict(Dict("$channel" => result_filteropt_dict))
    writelprops(data.par.rpars.fltopt[period], run, result)
    @info "Saved pars to disk"
end

function process_filteropt(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; kwargs...)
    # load meta data and configs 
    dsp_config = DSPConfig(data.metadata.config.dsp.dsp_config.default)
    peak =  Symbol(data.metadata.config.dsp.dsp_config.pz.default.peak)

    # load decay time for pole-zero correction 
    τ_pz = mvalue(get_values(data.par.rpars.pz[period, run, channel]).τ)
    @debug "Loaded decay time for pole-zero correction: $τ_pz"

    process_filteropt(data, period, run, category, channel, dsp_config, τ_pz, peak; kwargs...) 
end