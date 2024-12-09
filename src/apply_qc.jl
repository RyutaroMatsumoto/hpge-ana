
"""
    quality cuts based on dsp parameters
"""
function apply_qc(dsp_par::Union{PropDict, Table}, config_qc::PropDict)
    # get qc efficiencies
    qc_surv = Dict() 
    wvf_keep = (all = ones(Bool, length(dsp_par.t0)), )
    wvf_keep_all = ones(Bool, length(dsp_par.t0))

    for par in keys(config_qc)
        par_val = ustrip.(getproperty(dsp_par, par))
        par_val[par_val .== nothing] .= 0.0 * unit(par_val[1])
        local wvf_keep_tmp = ones(Bool, length(dsp_par.t0))
        # try 
        if haskey(getproperty(config_qc, par), :min)
            qc_surv_min = sum(par_val .> getproperty(config_qc, par).min)/length(dsp_par)  
            wvf_keep_tmp = wvf_keep_tmp .& (ustrip.(par_val) .> getproperty(config_qc, par).min)
        else
            qc_surv_min = 1 
        end
        if haskey(getproperty(config_qc, par), :max)
            qc_surv_max = sum(par_val .< getproperty(config_qc, par).max)/length(dsp_par)
            wvf_keep_tmp = wvf_keep_tmp  .& (ustrip.(par_val) .< getproperty(config_qc, par).max)
        else 
            qc_surv_max = 1
        end
        if haskey(getproperty(config_qc, par), :absmin)
            qc_surv_absmin = sum(abs.(par_val) .> getproperty(config_qc, par).absmin)/length(dsp_par)
            wvf_keep_tmp = wvf_keep_tmp .& (abs.(ustrip.(par_val)) .> getproperty(config_qc, par).absmin)
        else
            qc_surv_absmin = 1 
        end
        if haskey(getproperty(config_qc, par), :absmax)
            qc_surv_absmax = sum(abs.(par_val) .< getproperty(config_qc, par).absmax)/length(dsp_par)
            wvf_keep_tmp = wvf_keep_tmp .& (abs.(ustrip.(par_val)) .< getproperty(config_qc, par).absmax)
        else 
            qc_surv_absmax = 1
        end
        
        wvf_keep_all = wvf_keep_all .& wvf_keep_tmp
        wvf_keep = merge(wvf_keep, NamedTuple{(Symbol(par),)}((wvf_keep_tmp, )))
        qc_surv[par] = Dict(:min => qc_surv_min, :max => qc_surv_max, :absmax => qc_surv_absmax, :absmin => qc_surv_absmin)
    end
    qc_surv["all"] = sum(wvf_keep_all)/length(wvf_keep_all)
    wvf_keep = merge(wvf_keep, (all = wvf_keep_all,) )
    cuts = PropDict(Dict(:qc_surv => qc_surv, :wvf_keep => wvf_keep))
    return cuts
end 