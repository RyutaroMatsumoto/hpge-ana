function rms(x::Vector{T}) where {T <: Real} 
    sqrt(mean(x.^2))
end

# filter optimization: rise-time, and flat-top times
function filteropt_rt_blnoise(wvfs::ArrayOfRDWaveforms, filter_type::Symbol, ft::Quantity{T}, 
     grid_rt::StepRangeLen{Quantity{<:Float64}}, bl_window::ClosedInterval{<:Unitful.Time{T}}, flt_length_cusp::Quantity{T}; 
     cusp_scale::T = 3750.0, τ_cusp::Quantity{<:AbstractFloat} = 10000000.0u"µs") where T<:Real

    wvfs_timestep = diff(wvfs[1].time)[1]
    # STEP 1: rise-time optimization --> min. baseline noise after filtering 
    bl_start = findfirst(wvfs[1].time .>= leftendpoint(bl_window))
    bl_stop = findfirst(wvfs[1].time .>= rightendpoint(bl_window))
    noise = zeros(length(grid_rt))
    for (i, rt) in enumerate(collect(grid_rt))
        binscut = ceil(Int, (rt + ft /2)/wvfs_timestep)
        bin_start = bl_start + binscut
        bin_stop = bl_stop - binscut
        if filter_type == :trap
            wvfs_flt =  TrapezoidalChargeFilter(rt, ft).(wvfs)   # filtered waveforms 
        elseif filter_type == :cusp
            wvfs_flt = CUSPChargeFilter(rt, ft, τ_cusp, flt_length_cusp, cusp_scale).(wvfs) 
        end
        bl_trap = filter.(isfinite, map(x-> x.signal[bin_start:1:bin_stop], wvfs_flt)) # baselines - cut bins in beginning and end 
        noise[i] = rms(vcat(bl_trap...))
    end

    # find optimal shaping time: 
    rt_opt_rough = grid_rt[findfirst(noise .== minimum(noise))] # rt = 0.935 µs for ft = 1.5 µs; rt = 0.8 µs ft = 2 µs 
    f_interp = interpolate(ustrip.(collect(grid_rt)), noise, BSplineOrder(4))
    result = optimize(f_interp, 0.5, 5.0)
    rt_opt = Optim.minimizer(result)*u"µs" # optimial rise time for trap filter
    noise_min = Optim.minimum(result)
    
    result = (grid_rt = grid_rt, noise = noise, ft = ft, rt_opt = rt_opt, noise_min = noise_min, )
    return result
end

function filteropt_rt_blnoise(wvfs::ArrayOfRDWaveforms, filter_type::Symbol, dsp_config::DSPConfig; kwargs...)
    (_, ft)        = get_fltpars(PropDict(), filter_type, dsp_config)
    grid_rt   = getproperty(dsp_config, Symbol("e_grid_rt_$(filter_type)"))
    bl_window = dsp_config.bl_window
    flt_length_cusp = getproperty(dsp_config, Symbol("flt_length_cusp"))
    filteropt_rt_blnoise(wvfs, filter_type, ft, grid_rt, bl_window, flt_length_cusp; kwargs...)
end 

function filteropt_rt_blnoise(wvfs::ArrayOfRDWaveforms, filter_type::Symbol,
    ft::Quantity{T}, grid_rt::StepRangeLen{Quantity{<:Float64}}, bl_window::ClosedInterval{<:Unitful.Time{T}}; kwargs...) where T<:Real
    if filter_type !== :cusp
        filteropt_rt_blnoise(wvfs, filter_type, ft, grid_rt, bl_window, 0.0u"µs"; kwargs...)
    else 
        @info "Cusp filter requires additional argument: flt_length_cusp"
    end
end