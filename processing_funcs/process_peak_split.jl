"""
    process_peak_split(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; reprocess::Bool = false)
Create peak files containting only waveforms from peaks in the calibration spectrum.
1. Read raw data
2. Find peaks in the calibration spectrum using rough energy estimate `data_ch.daqenergy`
3. Do simple DSP for peaks only 
4. apply quality cuts based on simple DSP
5. Save waveforms after QC to peak files
Input:
- data: LegendData object
- period: DataPeriod object
- run: DataRun object
- category: Symbol or DataCategory object, e.g. :cal for calibration 
channel: ChannelId for germanium detector 
reprocess: Bool, default is false
Output:
- peak files containting only waveforms from peaks in the calibration spectrum
"""
function process_peak_split(data::LegendData, period::DataPeriod, run::DataRun, category::Union{Symbol, DataCategory}, channel::ChannelId; 
                        source::Symbol = :co60, reprocess::Bool = false, 
                        ecal_config = data.metadata.config.energy.energy_config.default,
                        dsp_config = DSPConfig(data.metadata.config.dsp.dsp_config.default))
   
    qc_config = data.metadata.config.qc.qc_config.default
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

    h_uncals = Vector{Histogram}(undef, length(filekeys))
    peakpos = Vector{Vector{<:Real}}(undef, length(filekeys))
    for f in eachindex(filekeys) 
        filekey = filekeys[f]
        peak_file = peak_files[f]

        data_ch = read_ldata(data, DataTier(:raw), filekey, channel)
        wvfs = data_ch.waveform
        e_uncal = data_ch.daqenergy
        eventnumber = data_ch.eventnumber
        timestamp = data_ch.timestamp
        if source == :co60
            gamma_lines =  ecal_config.co60_lines
            gamma_names =  ecal_config.co60_names
            left_window_sizes = ecal_config.co60_left_window_sizes 
            right_window_sizes = ecal_config.co60_right_window_sizes 
        elseif source == :th228
            gamma_lines =  ecal_config.th228_lines
            gamma_names =  ecal_config.th228_names
            left_window_sizes = ecal_config.left_window_sizes 
            right_window_sizes = ecal_config.right_window_sizes 
        end
        # nbins_def =  length(e_uncal)/10 > 100 ?  round(Int, length(e_uncal)/10) : 100
        bin_width = get_friedman_diaconis_bin_width(filter(in(quantile(e_uncal, 0.1)..quantile(e_uncal, 0.9)), e_uncal))
        nbins = ceil(Int, (maximum(e_uncal)-minimum(e_uncal))/bin_width)
        h_uncals[f] = fit(Histogram, e_uncal, nbins=nbins)
        _, peakpos[f] = RadiationSpectra.peakfinder(h_uncals[f], σ=ecal_config.peakfinder_σ, backgroundRemove=true, threshold=ecal_config.peakfinder_threshold)
        if length(peakpos[f]) !== length(gamma_lines)
            error("Number of peaks found $(length(peakpos[f])); expected gamma lines $(length(gamma_lines)) \n you could try to modify peakfinder_threshold and/or peakfinder_σ")
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
            qc_cuts, _ = apply_qc(dsp_par, qc_config)
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
    return h_uncals, peakpos 
end
