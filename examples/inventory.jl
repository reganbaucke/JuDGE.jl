using Random
using JuMP
using JuDGE

if !isdefined(@__MODULE__, :JuDGE_MP_Solver)
    # Replace this with another file in `/solvers` as appropriate.
    include("solvers/setup_gurobi.jl")
end

function inventory(;
    depth = 4,
    degree = 4,
    price_array = nothing,
    visualise = false,
    risk = RiskNeutral(),
)
    mytree = narytree(depth, degree)

    if price_array == nothing
        Random.seed!(1000)
        price_array = rand(length(collect(mytree)))
    end

    price = Dict(zip(collect(mytree), price_array))

    function sub_problems(n)
        sp = JuMP.Model(JuDGE_SP_Solver)

        @state(
            sp,
            -50 <= Δstock <= 50,
            state_name = stock,
            lb = 0,
            ub = 200,
            initial = 0
        )

        @ongoingcosts(sp, 0.01 * stock)

        @variable(sp, buy >= 0)
        @variable(sp, sell >= 0)
        #@variable(sp, change, Int)
        #@constraint(sp, Δstock==10*change)
        @constraint(sp, trade, Δstock == buy - sell)

        @objective(sp, Min, (buy * 1.01) * price[n] - (sell / 1.01) * price[n])

        return sp
    end

    model = JuDGEModel(
        mytree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_MP_Solver,
        check = true,
        risk = risk,
    )

    JuDGE.solve(model, verbose = 1)
    JuDGE.print_expansions(model, onlynonzero = true, inttol = 10^-5)
    println("\nRe-solved Objective: " * string(resolve_subproblems(model)))

    if visualise
        solution = JuDGE.solution_to_dictionary(model)
        JuDGE.add_to_dictionary!(solution, price, :prices)

        custom_plots = Dict{Symbol,Tuple{String,String,String}}()
        custom_plots[:graph1] = (
            "<div id=\"plotly\"></div>",
            "plotly_graph",
            joinpath(@__DIR__, "plot_functions", "inventory.js"),
        )

        for leaf in JuDGE.get_leafnodes(mytree)
            solution[leaf][:custom_data] = Dict{Symbol,Any}()
            history = JuDGE.history(leaf)

            graphs = []

            graph = Dict{Symbol,Any}()
            graph[:x] = []
            graph[:y] = []
            for i in 1:length(history)
                node = history[i]
                push!(graph[:x], JuDGE.depth(node) + 1)
                push!(graph[:y], solution[node][:stock_master])
            end
            push!(graph[:x], 0.0)
            push!(graph[:y], 0.0)
            graph[:type] = "scatter"
            graph[:name] = "Stock"
            push!(graphs, graph)

            graph = Dict{Symbol,Any}()
            graph[:x] = []
            graph[:y] = []
            for i in 1:length(history)
                node = history[i]
                push!(graph[:x], JuDGE.depth(node) + 1)
                push!(graph[:y], solution[node][:prices])
                push!(graph[:x], JuDGE.depth(node))
                push!(graph[:y], solution[node][:prices])
            end
            graph[:yaxis] = "y2"
            graph[:type] = "scatter"
            graph[:name] = "Prices"
            push!(graphs, graph)

            solution[leaf][:custom_data][:graph1] = graphs
        end

        JuDGE.visualize_tree(mytree, solution, custom = custom_plots)
    end
    return JuDGE.get_objval(model)
end

if !isdefined(@__MODULE__, :running_tests) || !running_tests
    inventory(visualise = true, risk = RiskNeutral())
    inventory(visualise = true, risk = JuDGE.Risk(0.1, bound = 0.0))
end
