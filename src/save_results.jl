function get_results(m,s,c,v)


    if result_count(m) == 0
        return [s,c,v,"infeasible"]
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

    println("----------------------- RESULTS ------------------------")
    @printf "Bought PFAD biomass: %.3f\n" PFAD_biomass
    @printf "Bought AF biomass: %.3f\n" AF_biomass
    @printf "Bought CO biomass: %.3f\n" CO_biomass
    @printf "Total bought biomass: %.3f\n" total_biomass

    # Total biodiesel / demand

    total_biodiesel = value.(m[:z])
    @printf "Total biodiesel: %.3f\n" total_biodiesel

    demand = 1500
    @printf "Total demand: %.3f\n" demand

    # Quality

    quality = value.(m[:rho])
    @printf "Quality: %.3f\n" quality

    # Total cost

    prod_cost = value.(m[:prod_cost])
    trans_cost = value.(m[:trans_cost])
    process_cost = value.(m[:process_cost])
    hydro_cost = value.(m[:hydro_cost])

    @printf "Production cost: %.3f\n" prod_cost
    @printf "Transportation cost: %.3f\n" trans_cost
    @printf "Processing cost: %.3f\n" process_cost
    @printf "Hydrotreatment cost: %.3f\n" hydro_cost

    total_cost = objective_value(m)
    @printf "Total cost: %.3f\n" total_cost
    @printf "should equal: %.3f\n" prod_cost + trans_cost + process_cost + hydro_cost
    println("----------------------- RESULTS ------------------------")
    println()

    data = [s,c,v,
            PFAD_biomass, AF_biomass, CO_biomass, total_biomass,
            total_biodiesel, demand,
            quality,
            prod_cost, trans_cost, process_cost, hydro_cost, total_cost]

    return data

end
