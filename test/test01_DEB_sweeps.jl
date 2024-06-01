TAG = splitpath(@__FILE__)[end] |> x -> split(x, ".")[1]
begin 
    p = DEBParamCollection()
    out = DynamicHydra.simulator(p)

    plt = @df out plot(
        plot(:t, :S, ylabel = "S"), 
        plot(:t, :H, ylabel = "H"),
        plot(:t, :R, ylabel = "R"), 
        plot(:t, [diffvec(:I), diffvec(:I_p), diffvec(:I_emb), diffvec(:A)], ylabel = "Idot"), 
        plot(:t, :X_p, ylabel = "X_p"), 
        plot(:t, :X_emb, ylabel = "X_emb"),
        xlabel = "t", 
        size = (800,500), 
        layout = (2,4), 
        leftmargin = 5mm
    )
    savefig(plt, "plots/$(TAG)_basesim")
    display(plt)
end

begin # effect of food input
    norm(x) = x ./ sum(x)
    # prepare the plot
    plt = plot(
        layout = grid(1,3, widths = norm([2, 1, 1])),
        leg = false, 
        title = ["Growth" "Reproduction" "Food density"], 
        leftmargin = 5mm, bottommargin = 6mm, 
        size = (1200,350), 
        xlabel = "Time (d)"
        )

    # iterate over nutrient input concentrations
    let Xdot_in = 4800.
        for _ in 1:5
            Xdot_in /= 2
            # generate the predidction
            out = DynamicHydra.simulator(
                DEBParamCollection(
                    glb = GlobalParams(Xdot_in = Xdot_in, t_max = 56.), 
                    spc = SpeciesParams(K_X = 12e3))
                )

            # plot the trajectories
            @df out plot!(plt, :t, :S, ylabel = "S", subplot = 1, leg = :outertopleft, label = "Xdot_in = $(Xdot_in)") 
            @df out plot!(plt, :t, :R, ylabel = "R", subplot = 2)
            @df out plot!(plt, :t, :X_p ./ GlobalParams().V_patch, ylabel = "[X_p]", subplot = 3, 
                yscale = :log10
                )
        end
        hline!(plt, [DynamicHydra.calc_S_max(SpeciesParams())], linestyle = :dash, color = "gray", subplot = 1, label = "S_max")
        display(plt)
        savefig(plt, "plots/$(TAG)_Xdot_in_sweep.png")
    end
end


