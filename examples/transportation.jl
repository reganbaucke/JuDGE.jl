using DelimitedFiles
using JuMP
using JuDGE
using Test

if !isdefined(@__MODULE__, :JuDGE_MP_Solver)
    # Replace this with another file in `/solvers` as appropriate.
    include("solvers/setup_gurobi.jl")
end

function transportation(; visualize = false)
    mytree = narytree(5, 2)

    function invest_supply_cost(node)
        if node.parent == nothing
            return Dict(zip(supply_nodes, [1.0, 2.0]))
        else
            p = node.parent
            for i in 1:length(p.children)
                if p.children[i] == node
                    temp = deepcopy(invest_supply_cost(p))
                    for key in keys(temp)
                        temp[key] *= (0.2 * i + 0.8)
                    end
                    return temp
                end
            end
        end
    end

    function invest_arc_cost(node)
        if node.parent == nothing
            temp = []
            for i in supply_nodes
                for j in demand_nodes
                    push!(temp, (i, j))
                end
            end
            return Dict(zip(temp, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]))
        else
            p = node.parent
            for i in 1:length(p.children)
                if p.children[i] == node
                    temp = deepcopy(invest_arc_cost(p))
                    for key in keys(temp)
                        temp[key] *= (0.2 * i + 0.8)
                    end
                    return temp
                end
            end
        end
    end

    function demand(node)
        if node.parent == nothing
            return d_dict
        else
            p = node.parent
            for i in 1:length(p.children)
                if p.children[i] == node
                    temp = deepcopy(demand(p))
                    for key in keys(temp)
                        temp[key] += i * i
                    end
                    return temp
                end
            end
        end
    end

    function supply(node)
        return s_dict
    end

    data_file = "transportation.csv"
    data = readdlm(joinpath(@__DIR__, data_file), ',')

    supply_nodes = data[3:end, 2]
    s = data[3:end, 1]

    demand_nodes = collect(data[2, 3:end])
    d = collect(data[1, 3:end])

    c = data[3:end, 3:end]

    # Converting arrays to dictionaries
    s_dict = Dict(zip(supply_nodes, s * 2))
    d_dict = Dict(zip(demand_nodes, d))

    c_dict = Dict()
    for i in 1:length(supply_nodes)
        for j in 1:length(demand_nodes)
            c_dict[supply_nodes[i], demand_nodes[j]] = c[i, j]
        end
    end

    ### with judge
    function sub_problems(node)
        model = Model(JuDGE_SP_Solver)

        @expansion(model, new_supply[supply_nodes], Bin) #invest in more supply
        @expansion(
            model,
            0 <= new_capacity[supply_nodes, demand_nodes] <= 2,
            Int
        ) #invest in more arc capacity

        @capitalcosts(
            model,
            sum(
                invest_supply_cost(node)[i] * new_supply[i] for
                i in supply_nodes
            ) + sum(
                invest_arc_cost(node)[i, j] * new_capacity[i, j] for
                i in supply_nodes for j in demand_nodes
            )
        )

        @variable(model, x[supply_nodes, demand_nodes] >= 0)

        @objective(
            model,
            Min,
            sum(
                c_dict[i, j] * x[i, j] for i in supply_nodes, j in demand_nodes
            )
        )
        @constraint(
            model,
            SupplyIncrease[i in supply_nodes],
            sum(x[i, j] for j in demand_nodes) <=
            supply(node)[i] + s_dict[i] * new_supply[i]
        )

        @constraint(
            model,
            CapacityIncrease[i in supply_nodes, j in demand_nodes],
            x[i, j] <= c_dict[i, j] + c_dict[i, j] * new_capacity[i, j]
        )

        @constraint(
            model,
            demandCon[j in demand_nodes],
            sum(x[i, j] for i in supply_nodes) == demand(node)[j]
        )

        return model
    end

    function format_output(s::Symbol, values)
        if s == :new_capacity
            output = Dict{Tuple,Float64}()
            for i in supply_nodes
                for j in demand_nodes
                    output[i, j] = values[i, j] * c_dict[i, j]
                end
            end
            return output
        elseif s == :new_supply
            output = Dict{String,Float64}()
            for i in supply_nodes
                output[i] = values[i] * s_dict[i]
            end
            return output
        end
        return nothing
    end

    judy = JuDGEModel(
        mytree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_MP_Solver,
        discount_factor = 0.9,
    )
    JuDGE.solve(judy, termination = Termination(inttol = 10^-7))

    println(
        "\nObjective: " * string(objective_value(judy.master_problem)) * "\n",
    )
    JuDGE.print_expansions(judy, format = format_output)

    println("\nRe-solved Objective: " * string(resolve_subproblems(judy)))
    solution = JuDGE.solution_to_dictionary(judy)

    if visualize
        JuDGE.visualize_tree(mytree, solution)
        JuDGE.write_solution_to_file(
            judy,
            joinpath(@__DIR__, "transport_solution_decomp.csv"),
        )
    end

    deteq = DetEqModel(
        mytree,
        ConditionallyUniformProbabilities,
        sub_problems,
        JuDGE_DE_Solver,
        discount_factor = 0.9,
    )
    JuDGE.solve(deteq)
    println(
        "Deterministic Equivalent Objective: " *
        string(objective_value(deteq.problem)),
    )
    JuDGE.write_solution_to_file(
        deteq,
        joinpath(@__DIR__, "transport_solution_deteq.csv"),
    )
    return objective_value(judy.master_problem)
end

@test transportation(visualize = false) ≈ 1924.35 atol = 1e-2
