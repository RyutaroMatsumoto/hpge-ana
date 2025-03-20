using LegendDataManagement
using LegendDataManagement: readlprops
using LegendDataManagement.LDMUtils
using LegendHDF5IO
using LegendSpecFits
using LegendDSP
using LegendDSP: get_fltpars
using RadiationDetectorSignals
using RadiationDetectorDSP
using Measurements: value as mvalue
using PropDicts
using StatsBase, IntervalSets
using Unitful
using TypedTables
using Makie, LegendMakie, CairoMakie
using Measures
using Optim
using BSplineKit
using Printf



result_rt, report_rt = filteropt_rt_optimization_blnoise(filter_type, wvfs, dsp_config, τ_pz; ft = def_ft)
           
            # Plots_theme()
            rt_inter = range(ustrip.(report_rt.enc_grid_rt[1]), stop = ustrip(maximum(report_rt.enc_grid_rt[findall(isfinite.(report_rt.enc))])), step = 0.05); 
            p = Figure()
            ax = Axis(p[1, 1], 
                xlabel = "Rise time ($(unit(report_rt.rt)))", ylabel = "Noise (a.u.)",
                limits = ((ustrip.(extrema(report_rt.enc_grid_rt))[1] - 0.2, ustrip.(extrema(report_rt.enc_grid_rt))[2] + 0.2), (nothing, nothing)),
                title = "Noise sweep ($filter_type), $period-$run-$channel, $peak peak \n" * @sprintf("fixed ft = %.2f %s, optimal rt = %.1f %s", ustrip(def_ft), unit(def_ft), ustrip(report_rt.rt), unit(report_rt.rt)), )
            lines!(ax, rt_inter, report_rt.f_interp.(rt_inter), color = :deepskyblue2, linewidth = 3, linestyle = :solid, label = "Interpolation")
            Makie.scatter!(ax, ustrip.(collect(report_rt.enc_grid_rt)), report_rt.enc,  color = :black, label = "Data")
            # p = Plots.plot(rt_inter, report_rt.f_interp.(rt_inter), color = :deepskyblue2, linewidth = 3, linestyle = :solid, label = "Interpolation")
        #    scatter!(ax, ustrip.(collect(report_rt.enc_grid_rt)), report_rt.enc, 
        #         size = (600, 400), dpi = 150, 
        #         ms = 4, color = :black, markerstrokecolor = :black,
        #         xlabel = "Rise time ($(unit(report_rt.rt)))", ylabel = "Noise (a.u.)", 
        #         title = "Noise sweep ($filter_type), $period-$run-$channel, $peak peak \n" * @sprintf("fixed ft = %.2f %s, optimal rt = %.1f %s", ustrip(def_ft), unit(def_ft), ustrip(report_rt.rt), unit(report_rt.rt)),
        #         label = "Data",
        #         legend = :top)
            # Plots.xlims!(ustrip.(extrema(report_rt.enc_grid_rt))[1] - 0.2, ustrip.(extrema(report_rt.enc_grid_rt))[2] + 0.2)
            axislegend()
            pname = plt_folder * split(LegendDataManagement.LDMUtils.get_pltfilename(data, filekeys[1], channel, Symbol("noise_sweep_$(filter_type)_blnoise")),"/")[end]
            d = LegendDataManagement.LDMUtils.get_pltfolder(data, filekeys[1], Symbol("noise_sweep_$(filter_type)_blnoise"))
            ifelse(isempty(readdir(d)), rm(d), nothing )

            result_energy = (
                m_cal_simple = m_cal_simple,
                fwhm = result_fwhm,
                cal = result_calib,
                fit  = result_fit,
        )
        end
    
for e_type in e_types
    result_dict[e_type] = _energy_calibration(e_type);
end 

result = PropDict(Dict("$channel" => result_dict))
writelprops(asic.par.rpars.ecal[period], run, result)
@info "Saved pars to disk"

plt_folder = LegendDataManagement.LDMUtils.get_pltfolder(asic, filekey, :noisesweep) * "/"
if !isdir(plt_folder)
    mkdir(plt_folder)
end
plt_name = plt_folder * _get_pltfilename(asic, filekey, channel, Symbol("noisesweep_$(e_type)"))
save(plt_name, fig)
