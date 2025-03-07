
function Makie_theme(; fs = 20, xgridvisible = true, ygridvisible = true)
    PlotTheme = Theme(
        size = (600, 420),
        fontsize = fs,
        dpi = 300,
        margin = 3mm,
        Fonts = (
            regular="Helvetica", math="Helvetica"),
        Axis = (
            titlefont = :regular,
            xgridvisible = xgridvisible,
            ygridvisible = ygridvisible,
            xlabelsize= fs + 6,
            ylabelsize= fs + 6,
            xtickalign=1,
            ytickalign=1,
        ),
        Legend = (
            framecolor = :silver,
            labelsize = fs + 2,
            position = :lt,
        )
    )
    set_theme!(PlotTheme)
end

function Plots_theme(; fs::Int = 16, grid::Symbol = :off)
    default(
        size = (650, 400),
        frame= :semi, 
        dpi = 300,
        margin = 3mm,
        grid = grid,
        xguidefontsize = fs, xtickfontsize = fs - 4 ,
        yguidefontsize = fs, ytickfontsize = fs - 4,
        legendfontsize = fs - 4 , titlefontsize = fs - 6,
        legend = :best, 
        legendforegroundcolor = :silver
    )
end

function LegendMakie.lplot!(report::NamedTuple{(:h_calsimple, :h_uncal, :c, :peak_guess, :peakhists, :peakstats)}; kwargs...)
    report_rn = (h_calsimple = report.h_calsimple, 
                        h_uncal = report.h_uncal,
                        c = report.c,
                        fep_guess = report.peak_guess,
                        peakhists = report.peakhists,
                        peakstats = report.peakstats)
    LegendMakie.lplot!(report_rn; kwargs...) 
end