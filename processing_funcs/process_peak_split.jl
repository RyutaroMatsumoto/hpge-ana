using Plots 
using LegendSpecFits
"""
    process_peak_split(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId, ecal_config::PropDict, dsp_config::DSPConfig, qc_config::PropDict; reprocess::Bool = false)
    process_peak_split(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; kwargs...)
Create peak files containting only waveforms from peaks in the calibration spectrum.
1. Read raw data
2. Find peaks in the calibration spectrum using rough energy estimate `data_ch.daqenergy`
3. Do simple DSP for peaks only 
4. apply quality cuts based on simple DSP
5. Save waveforms after QC to peak files
## inputs:
- data: LegendData object
- period: DataPeriod object
- run: DataRun object
- category: Symbol or DataCategory object, e.g. :cal for calibration 
- channel: ChannelId for germanium detector 
- ecal_config: PropDict, energy calibration configuration
- dsp_config: DSPConfig, DSP configuration
- qc_config: PropDict, quality cut configuration

## keyword arguments:
reprocess: Bool, default is false
## Output:
- h_uncals histograms of peaksearch
- peakpos 
Also, peak files containting only waveforms from peaks in the calibration spectrum are saved to jlpeaks folder.
"""
function process_peak_split(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId, ecal_config::PropDict, dsp_config::DSPConfig, qc_config::PropDict; reprocess::Bool = false, plotHist::Bool = true)
                        
    filekeys = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])
    peak_folder =  asic.tier[DataTier(:jlpeaks), category , period, run] * "/"
    if !ispath(peak_folder)
        mkpath(peak_folder)
        @info "create path: $peak_folder"
    end
    peak_files = peak_folder .* string.(filekeys) .* "-tier_jlpeaks.lh5"

    if reprocess == false
        if sum(isfile.(peak_files)) == length(filekeys)
            @info "all peak files are already processed. You're already done!"
            println("to read peak files: use read_ldata(data, :jlpeaks, category , period, run, channel)")
            return 
        else
            @info "$(sum(isfile.(peak_files))) peak files are already processed - skip"
        end
        filekeys = filekeys[.!isfile.(peak_files)]
    end

    if Symbol(ecal_config.source) == :co60
        gamma_lines =  ecal_config.co60_lines
        gamma_names =  ecal_config.co60_names
        left_window_sizes = ecal_config.co60_left_window_sizes 
        right_window_sizes = ecal_config.co60_right_window_sizes 
    elseif Symbol(ecal_config.source) == :th228
        gamma_lines =  ecal_config.th228_lines
        gamma_names =  ecal_config.th228_names
        left_window_sizes = ecal_config.left_window_sizes 
        right_window_sizes = ecal_config.right_window_sizes 
    end

    h_uncals = Vector{Histogram}(undef, length(filekeys))
    peakpos = Vector{Vector{<:Real}}(undef, length(filekeys))
    for f in eachindex(filekeys) 
        filekey = filekeys[f]
        peak_file = peak_files[f]
        data_ch = read_ldata(data, DataTier(:raw), filekey, channel)
        wvfs = data_ch.waveform
        e_uncal = filter(x -> x >= qc_config.e_trap.min , data_ch.daqenergy)
        if isempty(e_uncal)
            @warn "No energy values >= $(qc_config.e_trap.min) found for $filekey - skip"
            continue
        end
        eventnumber = data_ch.eventnumber
        timestamp = data_ch.timestamp
       
          # binning and peak search windows and histogram settings 
        bin_min = quantile(e_uncal, ecal_config.left_bin_quantile)
        bin_max = quantile(e_uncal, ecal_config.right_bin_quantile)
        peak_min = quantile(e_uncal, ecal_config.left_peak_quantile)
        peak_max = quantile(e_uncal, ecal_config.right_peak_quantile)
        bin_width = LegendSpecFits.get_friedman_diaconis_bin_width(filter(in(bin_min..bin_max), e_uncal))
        if (peak_max-peak_min)/bin_width < ecal_config.nbins_min
            bin_width = (peak_max-peak_min)/ ecal_config.nbins_min
        end

        # peak search
        h_uncals[f] = fit(Histogram, e_uncal, 0:bin_width:maximum(e_uncal)) # histogram over full energy range; stored for plot 
        try
            h_peaksearch = fit(Histogram, e_uncal, 0:bin_width:peak_max) # histogram for peak search
            _, peakpos[f] = RadiationSpectra.peakfinder(h_peaksearch, σ= ecal_config.peakfinder_σ, backgroundRemove=true, threshold = ecal_config.peakfinder_threshold)
        catch e 
            @warn "peakfinder failed for $filekey - use larger window. err $e"
            h_peaksearch = fit(Histogram, e_uncal, 0:bin_width:(peak_max*1.5)) # histogram for peak search
            _, peakpos[f] = RadiationSpectra.peakfinder(h_peaksearch, σ= ecal_config.peakfinder_σ, backgroundRemove=true, threshold = ecal_config.peakfinder_threshold)
        end 
        if length(peakpos[f]) !== length(gamma_lines)
            error("Number of peaks found $(length(peakpos[f])); expected gamma lines $(length(gamma_lines)) \n you could try to modify peakfinder_threshold and/or peakfinder_σ (file $f)")
        else 
            @info "Found $(length(peakpos[f])) peaks for $filekey"
        end 
        cal_simple = mean(gamma_lines./sort(peakpos[f]))
        e_simplecal = e_uncal .* cal_simple
       
        # save
        if !ispath(peak_folder)
            mkpath(peak_folder)
            @info "create path: $peak_folder"
        end

        fid = lh5open(peak_file, "w") 
        for i in eachindex(gamma_lines)
            peakIdx = findall((gamma_lines[i] - left_window_sizes[i]) .<= e_simplecal .< (gamma_lines[i] + right_window_sizes[i]))
            # do simple dsp. only for peaks. 
            dsp_par = simple_dsp_qc(Table(waveform = wvfs[peakIdx]), dsp_config)
            qc_cuts = apply_qc(dsp_par, qc_config)
            qc_flag = qc_cuts.wvf_keep.all
            qc_idx = findall(qc_flag)
            qc_surv = qc_cuts.qc_surv.all
            fid["$channel/jlpeaks/$(gamma_names[i])/waveform"] = wvfs[peakIdx][qc_idx]
            fid["$channel/jlpeaks/$(gamma_names[i])/daqenergy"] = e_uncal[peakIdx][qc_idx]
            fid["$channel/jlpeaks/$(gamma_names[i])/e_simplecal"] = e_simplecal[peakIdx][qc_idx]
            fid["$channel/jlpeaks/$(gamma_names[i])/eventnumber"] = eventnumber[peakIdx][qc_idx]
            fid["$channel/jlpeaks/$(gamma_names[i])/timestamp"] = timestamp[peakIdx][qc_idx]
            fid["$channel/jlpeaks/$(gamma_names[i])/gamma_line"] = fill(gamma_lines[i], length(qc_idx))
            fid["$channel/jlpeaks/$(gamma_names[i])/qc_flag"] = qc_flag[qc_idx]
            fid["$channel/jlpeaks/$(gamma_names[i])/qc_surv"] = fill(qc_surv, length(qc_idx))
            for par in columnnames(dsp_par)
                fid["$channel/jlpeaks/$(gamma_names[i])/$par"]  = getproperty(dsp_par, par)
            end
        end
        println("peak file processing done for $filekey")
        close(fid)
    end
    if plotHist == true 
        @info "plot histograms"
      
        plts = Vector{Plots.Plot}(undef, length(filekeys))
        Plots_theme()
        for f in eachindex(filekeys) 
            if isassigned(h_uncals, f)
                plts[f] = Plots.plot(h_uncals[f], size = (600,400), xlabel = "Energy (a.u.)", linecolor = :dodgerblue, color = :dodgerblue, label = false)
                title!(get_plottitle(filekeys[f], _channel2detector(data, channel), "peak split"))
                vline!([peakpos[f]], color = :red2, label = "peakfinder result", linewidth = 2, linestyle = :dash)
                plt_folder = LegendDataManagement.LDMUtils.get_pltfolder(data, filekeys[f], :peak_split) * "/"
                pname = plt_folder * "peak_split_$(filekeys[f]).png"
                savefig(plts[f], pname)
            end 
        end 
        return plts 
    end 
   
end

function process_peak_split(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; kwargs...)
    filekeys = search_disk(FileKey, data.tier[DataTier(:raw), category , period, run])
    qc_config = dataprod_config(data).qc(filekeys[1]).default
    ecal_config = dataprod_config(data).energy(filekeys[1]).default
    dsp_config = DSPConfig(dataprod_config(data).dsp(filekeys[1]).default)
    @info "use default configs: qc, ecal, dsp"
    process_peak_split(data, period, run, category, channel, ecal_config, dsp_config, qc_config; kwargs...)
end

