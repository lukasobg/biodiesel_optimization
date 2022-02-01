function get_results(m,file)

    # Total biomass bought (sum of x)

    x_sb = value.(m[:x])

    PFAD_biomass = 0
    AF_biomass = 0
    CO_biomass = 0
    for i in 1:length(x_sb[:,1])
        PFAD_biomass += x_sb[i,3]
        AF_biomass += x_sb[i,2]
        CO_biomass += x_sb[i,1]
    end
    total_biomass = PFAD_biomass + AF_biomass + CO_biomass

    println("Bought PFAD biomass: $(PFAD_biomass)")
    println("Bought AF biomass: $(AF_biomass)")
    println("Bought CO biomass: $(CO_biomass)")
    println("Total bought biomass: $(total_biomass)")

    # Total biodiesel / demand

    total_biodiesel = value.(m[:z])
    println("Total biodiesel: $(total_biodiesel)")

    demand = 1500
    println("Total demand: $(demand)")

    # Quality

    quality = value.(m[:rho])
    println("Quality: $(quality)")

    # Total cost

    total_cost = objective_value(m)
    println("Total cost: $(total_cost)")

    data = [(PFAD_biomass, AF_biomass, CO_biomass, total_biomass,
             total_biodiesel, demand,
             quality, total_cost)]

    open(file; write=true) do f
         write(f, "PFAD;AF;CO;total biom;prod biod;demand;quality;total cost\n")
         writedlm(f, data, ';')
    end

end
