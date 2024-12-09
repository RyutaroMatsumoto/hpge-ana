# function process_filteropt(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId)
    dsp_config = DSPConfig(data.metadata.config.dsp.dsp_config.default)
    pz_config = data.metadata.config.dsp.dsp_config.pz.default
    energy_config = data.metadata.config.energy.energy_config.default
    peak = Symbol(pz_config.peak)
    τ_pz = mvalue(data.par.rpars.pz[period, run, channel].τ)

    @debug "Create pars db"
    mkpath(joinpath(data_path(data.par.rpars.fltopt), string(period)))

    # load waveforms from peak 
    data_peak  = read_ldata((peak), asic, :jlpeaks, category, period, run, channel)
    wvfs = data_peak.waveform

    # waveforms: shift baseline and deconvolute (pole-zero)
    bl_window  = dsp_config.bl_window
    bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))
    wvfs = shift_waveform.(wvfs, -bl_stats.mean)
    deconv_flt = InvCRFilter(τ_pz)
    wvfs = deconv_flt.(wvfs)

    # filter optimization: rise-time, and flat-top times
    filter_types = [:trap, :cusp] 
    filter_type = :trap
    def_rt, def_ft = get_fltpars(PropDict(), filter_type, dsp_config)
    e_grid_rt   = getproperty(dsp_config, Symbol("e_grid_rt_$(filter_type)"))
    e_grid_ft   = getproperty(dsp_config, Symbol("e_grid_ft_$(filter_type)"))
    wvfs_timestep = diff(wvfs[1].time)[1]
    # if filter_type == :cusp
        flt_length_cusp             = dsp_config.flt_length_cusp
        cusp_scale                  = ustrip(NoUnits, flt_length_cusp/wvfs_timestep)
        τ_cusp = 10000000.0u"µs"
    # end

    # STEP 1: rise-time optimization --> min. baseline noise after filtering 
    result = filteropt_rt_blnoise(wvfs, filter_type, dsp_config; cusp_scale = cusp_scale, τ_cusp = τ_cusp)
    plot(ustrip.(collect(result.grid_rt)), result.noise )
    

    # STEP 2: flat-top optimization 
    signal_estimator = SignalEstimator(PolynomialDNI(dsp_config.kwargs_pars.sig_interpolation_order, dsp_config.kwargs_pars.sig_interpolation_length))
    t50 = get_threshold(wvfs, maximum.(wvfs.signal) .* 0.5; mintot=dsp_config.kwargs_pars.tx_mintot)
    fwhm = zeros(length(collect(e_grid_ft)), 1).*u"keV"
    for (i, ft) in enumerate(e_grid_ft)
    (i, ft) = (1, 1.5u"µs")
        if filter_type == :trap
            uflt_trap_rtft = TrapezoidalChargeFilter(rt_opt, ft)
            e_ftr = signal_estimator.(uflt_trap_rtft.(wvfs), t50 .+ (rt_opt + ft/2))
        elseif filter_type == :cusp
            uflt_cusp_rtft = CUSPChargeFilter(rt_opt, ft, τ_cusp, flt_length_cusp, cusp_scale)
            e_ftr = signal_estimator.(uflt_cusp_rtft.(wvfs), t50 .+ (flt_length_cusp /2))    
        end 
        
        e_simple = filter(x-> x > (data_peak.gamma_line[1] - 100u"keV"), e_ftr .* data_peak.gamma_line[1]/median(filter(isfinite,e_ftr)))
        window_sizes = [(energy_config.co60_left_window_sizes[1], energy_config.co60_right_window_sizes[1])]
        peakhists, peakstats, h_calsimple, bin_widths = get_peakhists_th228(e_simple, [data_peak.gamma_line[1]], window_sizes; binning_peak_window = 2.5u"keV", e_unit=u"keV");
        result_fit, report_fit = LegendSpecFits.fit_peaks(peakhists, peakstats, [data_peak.gamma_line[1]]; calib_type = :th228, fit_func = [:gamma_def], uncertainty = false)
        fwhm[i] = mvalue([result_fit[k].fwhm for k in keys(result_fit)][1]) .*u"keV"
    end
    scatter(ustrip.(collect(e_grid_ft)), fwhm)



    # result_pz = (τ = result.μ, fit = result)
    # writelprops(data.par.rpars.pz[period], run, PropDict("$channel" => result_pz))
    # @info "Saved pars to disk"
# end