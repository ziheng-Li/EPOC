using Plots, Interpolations

# include("../ooeoae_expts.jl")

######################################
# key parameters that control plot appearance
include("expt_plot_endP.jl")
##########################################
isdir("figures") || mkdir("figures")
dropbox_output_dir = "figures"

output_figures_dir = joinpath(dropbox_output_dir, "Data_compil")
isdir(output_figures_dir) || mkdir(output_figures_dir)
##########################################

(plot_PTB_C, plot_PTB_U, plot_PTB_TEMP, df_C, df_U, df_TEMP) = plot_data_C_U_O("PTB-C-U-O.xlsx")

p_d13C, new_Age, new_data_mean, new_data_min, new_data_max = box_plot(df_C, :d13Ccarb, 0.1e6; start_age=-254e6, end_age=-244e6, color=:skyblue, marker=:circle, ylims=(-9,9), ylabel="δ¹³Ccarb (‰)", data_alpha=0.05, plot_lines=true)

p_d238U, new_Age, new_data_mean, new_data_min, new_data_max = box_plot(df_U, :d238U, 0.4e6; start_age=-254e6, end_age=-244e6, color=:orange, marker=:diamond, ylims=(-1.2,0.3), ylabel="δ²³⁸Ucarb (‰)", data_alpha=0.6)

p_d18O, new_Age, new_data_mean, new_data_min, new_data_max = box_plot(df_TEMP, :d18O, 0.4e6; start_age=-254e6, end_age=-244e6, color=:green, marker=:circle, ylims=(17,23), ylabel="δ¹⁸Oapatite (‰)", data_alpha=0.6, yflip=true)

# sum plots
l = @layout[
            a{0.3h};
            b;
            c
    ]

p_sum = Plots.plot(p_d13C, p_d238U, p_d18O,  
                        layout = l, left_margin = 5Plots.mm, right_margin = 5Plots.mm, size=(400, 500))
# savefig(p_sum, joinpath(output_figures_dir, "POAU_data_compil.svg"))  
