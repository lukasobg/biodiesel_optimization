function get_results(m,c,v)


    if result_count(m) == 0
        return [c,v,"infeasible"]
    end

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

    prod_cost = value.(m[:prod_cost])
    trans_cost = value.(m[:trans_cost])
    process_cost = value.(m[:process_cost])
    hydro_cost = value.(m[:hydro_cost])

    println("Production cost: $(prod_cost)")
    println("Transportation cost: $(trans_cost)")
    println("Processing cost: $(process_cost)")
    println("Hydrotreatment cost: $(hydro_cost)")

    total_cost = objective_value(m)
    println("Total cost: $(total_cost)")
    println("should equal: $(prod_cost + trans_cost + process_cost + hydro_cost)")

    data = [c,v,
            PFAD_biomass, AF_biomass, CO_biomass, total_biomass,
            total_biodiesel, demand,
            quality,
            prod_cost, trans_cost, process_cost, hydro_cost, total_cost]

    return data

end
