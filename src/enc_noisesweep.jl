function rms(x::Vector{T}) where {T <: Real} 
    sqrt(mean(x.^2))
end

"""
"""
function enc_noisesweep(filter_type::Symbol, wvfs::ArrayOfRDWaveforms, dsp_config::DSPConfig; 
             ft::Quantity{T}= 0.0u"µs", τ_cusp::Quantity{<:AbstractFloat} = 10000000.0u"µs", τ_zac::Quantity{<:AbstractFloat} = 10000000.0u"µs" ) where T<:Real

    # gather config parameters 
    grid_rt   = getproperty(dsp_config, Symbol("e_grid_rt_$(filter_type)"))
    cap = getproperty(energy_config, Symbol("C_Pulser"))
    temp = getproperty(energy_config, Symbol("Temp"))
    if ft == 0.0u"µs"
        (_, ft) = get_fltpars(PropDict(), :trap, dsp_config)
    end
    bl_window = dsp_config.bl_window

    # cusp-specific filter parameters 
    flt_length_cusp = getproperty(dsp_config, Symbol("flt_length_cusp"))
    cusp_scale = ustrip(NoUnits, flt_length_cusp/step(wvfs[1].time))

    # zac-specific filter parameters 
    flt_length_zac = getproperty(dsp_config, Symbol("flt_length_zac"))
    zac_scale = ustrip(NoUnits, flt_length_zac/step(wvfs[1].time))

    # waveforms: shift baseline and deconvolute (pole-zero)
    bl_stats = signalstats.(wvfs, leftendpoint(bl_window), rightendpoint(bl_window))
    wvfs_shift = shift_waveform.(wvfs, -bl_stats.mean)
 
    # calculate baseline rms after filtering 
    # filtered waveforms are already truncated to valid windows. no additional truncation needed.
    noise_au = zeros(length(grid_rt))
    noise = zeros(length(grid_rt))
    for (i, rt) in enumerate(collect(grid_rt))
        if filter_type == :trap
            wvfs_flt =  TrapezoidalChargeFilter(rt, ft).(wvfs_shift)   # filtered waveforms 
        elseif filter_type == :cusp
            wvfs_flt = CUSPChargeFilter(rt, ft, τ_cusp, flt_length_cusp, cusp_scale).(wvfs_shift) 
        elseif filter_type == :zac
            wvfs_flt = ZACChargeFilter(rt, ft, τ_zac, flt_length_zac, zac_scale).(wvfs_shift)
        end
        valid_bins = findall(leftendpoint(bl_window) .<=  wvfs_flt[1].time .<=rightendpoint(bl_window)) # bins within baseline of filtered waveform 
        bl_trap = filter.(isfinite, map(x-> x.signal[valid_bins], wvfs_flt)) # baselines - cut bins in beginning and end 
        noise_au[i] = rms(vcat(bl_trap...))
        noise[i] = pulser_ADC_to_keV(noise_au[i],cap,temp)
    end

    # find optimal shaping time: 
    if  all(.!isfinite.(noise))
        @error "Error: no finite noise values found for filter type $filter_type - try adjusting grid_rt."
    end
    result, report = let enc = noise[findall(isfinite.(noise))], rt = ustrip.(collect(grid_rt))[findall(isfinite.(noise))]
        if length(enc) >= 4
            f_interp = BSplineKit.interpolate(rt, enc, BSplineOrder(4))
        else
            f_interp = LinearInterpolation(rt, enc)
        end
        result = optimize(f_interp, minimum(rt), maximum(rt)) 
        rt_opt = Optim.minimizer(result)*u"µs" # optimial rise time for trap filter
        min_enc = Optim.minimum(result)
        (rt = rt_opt, min_enc = min_enc), (rt = rt_opt, min_enc = min_enc, enc_grid_rt = grid_rt, enc = noise, f_interp = f_interp)
    end
    return result, report 
end

