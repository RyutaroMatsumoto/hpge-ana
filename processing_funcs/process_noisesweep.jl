"""
    process_noisesweep(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId, dsp_config::DSPConfig, peak::Symbol;reprocess::Bool = true, filter_types::Vector{Symbol} = [:trap, cusp])
Noise sweep for filter_types
- load waveforms, shift baseline
- optimize rise-time for minimum baseline noise after filtering
- save results to disk
- sanity plots for noise noise_sweep_
Inputs: 
- `data::LegendData`: LegendData object
- `period::DataPeriod`: data period
- `run::DataRun`: data run
- `category::Union{Symbol, DataCategory}`: data category, e.g. :cal
- `channel::ChannelId`: channel id
Optional:
- `dsp_config::DSPConfig`: DSP configuration object. If not specified will take default from metadata
- `peak::Symbol`: peak to optimize for (needs existing peakfile!). If not specified will take default from metadata
- `reprocess::Bool`: reprocess the files or not
- `filter_types::Vector{Symbol}`: filter types to optimize for
"""

function process_noisesweep(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId, dsp_config::DSPConfig, peak::Symbol; 
                reprocess::Bool = false, filter_types::Vector{Symbol} = [:trap, :cusp, :zac])  
    det = _channel2detector(data, channel)
    @info "Noise sweep for period $period, run $run, channel $channel /det $det - $filter_types"

    # check if decaytime pars already exist
    noisesweep_file = joinpath(mkpath(data_path(data.par.rpars.noisesweep[period])), "$(string(run)).json")
    if isfile(noisesweep_file) && !reprocess
        @info "Noise sweep file already exist for $category period $period - run $run - channel $channel - you're done!"
        return
    end

    # prepare results dict
    mkpath(joinpath(data_path(data.par.rpars.noisesweep),string(period)))
    result_noisesweep_dict = Dict{Symbol, NamedTuple}()
    @debug "Created path for noisesweep results"

    # load waveforms from peak
    filekeys = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])
    if peak == :all 
        data_peak = read_ldata(data, DataTier(:raw), filekeys, channel)
        #data_peak = merge(data_peak, (gamma_line = [1170*u"keV"],))
    else
        data_peak  = read_ldata((peak), data, :jlpeaks, category, period, run, channel)
    end 
    wvfs = data_peak.waveform
    @debug "Loaded waveforms for peak $peak"

    function process_noisesweep_fltr(filter_type::Symbol)
        @info "Noise sweep for filter $filter_type"
        _, def_ft = get_fltpars(PropDict(), filter_type, dsp_config)
        filekey = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])[1]
        plt_folder = LegendDataManagement.LDMUtils.get_pltfolder(data, filekey, :noisesweep) * "/"
        if !isdir(plt_folder)
            mkdir(plt_folder)
        end
        # STEP 1: rise-time optimization --> min. baseline noise after filtering 
        result_rt, report_rt = enc_noisesweep(filter_type, wvfs, dsp_config; ft = def_ft)
        
        # Plots_theme()
        rt_inter = range(ustrip.(report_rt.enc_grid_rt[1]), stop = ustrip(maximum(report_rt.enc_grid_rt[findall(isfinite.(report_rt.enc))])), step = 0.05); 
        p = Figure()
        ax = Axis(p[1, 1], 
            xlabel = "Rise time ($(unit(report_rt.rt)))", ylabel = "Noise (a.u.)",
            limits = ((ustrip.(extrema(report_rt.enc_grid_rt))[1] - 0.2, ustrip.(extrema(report_rt.enc_grid_rt))[2] + 0.2), (nothing, nothing)),
            title = "Noise sweep ($filter_type), $period-$run-$channel, $peak peak \n" * @sprintf("fixed ft = %.2f %s, optimal rt = %.1f %s", ustrip(def_ft), unit(def_ft), ustrip(report_rt.rt), unit(report_rt.rt)), )
        lines!(ax, rt_inter, report_rt.f_interp.(rt_inter), color = :deepskyblue2, linewidth = 3, linestyle = :solid, label = "Interpolation")
        Makie.scatter!(ax, ustrip.(collect(report_rt.enc_grid_rt)), report_rt.enc,  color = :black, label = "Data")
        axislegend()
        pname = plt_folder * split(LegendDataManagement.LDMUtils.get_pltfilename(data, filekeys[1], channel, Symbol("noise_sweep_$(filter_type)")),"/")[end]
        d = LegendDataManagement.LDMUtils.get_pltfolder(data, filekeys[1], Symbol("noise_sweep_$(filter_type)"))
        ifelse(isempty(readdir(d)), rm(d), NamedTuple())

        save(pname, p)
        @info "Save sanity plot to $pname"

        # return result for this filter type
        #merge(result_rt, result_ft)
    end 

    # filter optimization: rise-time, and flat-top times
    for filter_type in filter_types
        result_noisesweep_dict[filter_type] =  process_noisesweep_fltr(filter_type)
    end
    result = PropDict("$channel" => result_noisesweep_dict)
    writelprops(data.par.rpars.noisesweep[period], run, result)
    @info "Saved pars to disk"
end

# function process_filteropt(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; kwargs...)
    
#     filekeys = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])
#     dsp_config = DSPConfig(dataprod_config(data).dsp(filekeys[1]).default)
#     pz_config = dataprod_config(data).dsp(filekeys[1]).pz.default

#     peak =  Symbol(pz_config.peak)
#     τ_pz = mvalue(get_values(data.par.rpars.pz[period, run, channel]).τ)
#     @debug "Loaded decay time for pole-zero correction: $τ_pz"

#     process_filteropt(data, period, run, category, channel, dsp_config, τ_pz, peak; kwargs...) 
# end