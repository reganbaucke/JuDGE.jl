module JuDGE

using Printf
using Distributed
using JuMP

include("tree.jl")
include("macros.jl")
include("convergence.jl")
include("risk.jl")
include("deteq.jl")

mutable struct Bounds
    UB::Float64
    LB::Float64
end

struct Column
    node::AbstractTree
    coeffs::Dict{Symbol,Any}
    obj::Float64
    var::VariableRef
    solution::Array{Float64,1}
end

struct JuDGEModel
    tree::AbstractTree
    master_problem::JuMP.Model
    sub_problems::Dict{AbstractTree,JuMP.Model}
    bounds::Bounds
    discount_factor::Float64
    master_solver::Any#::Union{DataType,MOI.OptimizerWithAttributes,}
    probabilities::Dict{AbstractTree,Float64}
    risk::Union{Risk,Array{Risk,1}}
    sideconstraints::Union{Nothing,Function}
    log::Array{ConvergenceState,1}
    ext::Dict{Symbol,Any}
end

include("utilities.jl")

function JuDGEModel(
    tree::T where {T<:AbstractTree},
    probabilities::Dict{AbstractTree,Float64},
    master_problem::JuMP.Model,
    sub_problems::Dict{AbstractTree,JuMP.Model},
    solver,
    bounds::Bounds,
    discount_factor::Float64,
    risk::Union{Risk,Array{Risk,1}},
    sideconstraints::Union{Function,Nothing},
    ext::Dict{Symbol,Any},
)
    new_ext = Dict{Symbol,Any}()
    new_ext[:branches] = copy(ext[:branches])
    new_ext[:optimizer_settings] = deepcopy(ext[:optimizer_settings])
    return JuDGEModel(
        tree,
        master_problem,
        sub_problems,
        bounds,
        discount_factor,
        solver,
        probabilities,
        risk,
        sideconstraints,
        Array{ConvergenceState,1}(),
        new_ext,
    )
end

"""
	JuDGEModel(tree::AbstractTree,
               probabilities,
               sub_problem_builder::Function,
               solver;
               discount_factor=1.0,
               risk=RiskNeutral(),
               sideconstraints=nothing,
               check=true,
               perfect_foresight=false)

Define a JuDGE model.

### Required arguments
`tree` is a reference to a scenario tree

`probabilities` is either a function, which returns a dictionary of the probabilities
of all nodes in a tree, or simply the dictionary itself

`sub_problem_builder` is a function mapping a node to a JuMP model for each subproblems

`solver` is a reference to the optimizer used for the master problem (with appropriate settings);
 this can also be a tuple containing two optimizers (one for solving the relaxation, and one for
 solving the binary model)

### Optional arguments
`discount_factor` is a number between 0 and 1 defining a constant discount factor along each arc
in the scenario tree

`risk` can be either a `Risk` object, or a vector of such objects.

`sideconstraints` is a function which specifies side constraints in the master problem, see
[Tutorial 9: Side-constraints](@ref) for further details

`check` is a boolean, which can be set to `false` to disable the validation of the JuDGE model.

`perfect_foresight` is a boolean; this is an experimental feature, which creates an array of
JuDGE models, one for each leaf node. This will enable users to easily compute the EVPI for
the stochastic program. Also can be used for regret-based risk implementations.

### Examples
	judge = JuDGEModel(tree, ConditionallyUniformProbabilities, sub_problems,
                                    Gurobi.Optimizer)
	judge = JuDGEModel(tree, probabilities, sub_problems, CPLEX.Optimizer,
                                    discount_factor=0.9, risk=Risk(0.5,0.1)))
"""
function JuDGEModel(
    tree::T where {T<:AbstractTree},
    probabilities,
    sub_problem_builder::Function,
    solver;
    discount_factor::Float64 = 1.0,
    risk::Union{Risk,Array{Risk,1}} = RiskNeutral(),
    sideconstraints::Union{Function,Nothing} = nothing,
    check::Bool = true,
    perfect_foresight::Bool = false,
)
    println("")
    if !perfect_foresight
        println("Establishing JuDGE model for tree: " * string(tree))
    else
        println(
            "Establishing perfect foresight models for tree: " * string(tree),
        )
    end

    if typeof(probabilities) <: Function
        probabilities = probabilities(tree)
    end
    if typeof(probabilities) != Dict{AbstractTree,Float64}
        error(
            "\'probabilities\' needs to be a dictionary mapping AbstractTree to Float64\nor a function that generates such a dictionary",
        )
    end

    nodes = collect(tree)
    sub_problems = Dict(i => sub_problem_builder(i) for i in nodes)

    if check
        print("Checking sub-problem format...")
        check_specification_is_legal(sub_problems)
        println("Passed")
    else
        println("Skipping checks of sub-problem format")
    end
    scale_objectives(tree, sub_problems, discount_factor)

    if !perfect_foresight
        print("Building master problem...")
        master_problem = build_master(
            sub_problems,
            tree,
            probabilities,
            solver,
            discount_factor,
            risk,
            sideconstraints,
        )
        println("Complete")
        ext = Dict{Symbol,Any}()
        ext[:branches] = Branch[]
        ext[:optimizer_settings] = Dict{Symbol,Any}()
        return JuDGEModel(
            tree,
            master_problem,
            sub_problems,
            Bounds(Inf, -Inf),
            discount_factor,
            solver,
            probabilities,
            risk,
            sideconstraints,
            Array{ConvergenceState,1}(),
            ext,
        )
    else
        scenarios = Dict{AbstractTree,JuDGEModel}()
        pr = Dict{AbstractTree,Float64}()
        print("Building master problems...")
        scen_trees = get_scenarios(tree)
        for t in scen_trees
            leaf = getID(get_leafnodes(t)[1])
            sps = Dict(i => sub_problems[getID(i)] for i in collect(t))
            probs = Dict(i => 1.0 for i in collect(t))
            master_problem = build_master(
                sps,
                t,
                probs,
                solver,
                discount_factor,
                risk,
                sideconstraints,
            )
            ext = Dict{Symbol,Any}()
            ext[:branches] = Branch[]
            ext[:optimizer_settings] = Dict{Symbol,Any}()
            scenarios[leaf] = JuDGEModel(
                t,
                master_problem,
                sps,
                Bounds(Inf, -Inf),
                discount_factor,
                solver,
                probabilities,
                risk,
                sideconstraints,
                Array{ConvergenceState,1}(),
                ext,
            )
            pr[leaf] = probabilities[leaf]
        end
        println("Complete")
        return scenarios, pr
    end
end

include("branchandprice.jl")
include("master.jl")

function scale_objectives(
    tree::T where {T<:AbstractTree},
    sub_problems,
    discount_factor::Float64,
)
    if discount_factor <= 0.0 || discount_factor > 1.0
        error(
            "discount_factor must be greater than 0.0 and less than or equal to 1.0",
        )
    end

    for (node, sp) in sub_problems
        if !haskey(sp.ext, :capitalcosts)
            sp.ext[:capitalcosts] = Dict()
            sp.ext[:capitalcosts][:constant] = 0
        end
        if !haskey(sp.ext, :ongoingcosts)
            sp.ext[:ongoingcosts] = Dict()
            sp.ext[:ongoingcosts][:constant] = 0
        end

        sp.ext[:objective] = @variable(sp, obj)
        sp.ext[:objective_con] =
            @constraint(sp, sp.ext[:objective] - objective_function(sp) == 0)

        @objective(sp, Min, 0.0)
        set_normalized_coefficient(
            sp.ext[:objective_con],
            sp.ext[:objective],
            1.0 / (discount_factor^depth(node)),
        )
    end
end

function Base.map(f, hello::JuMP.Containers.DenseAxisArray)
    return JuMP.Containers.DenseAxisArray(
        map(f, hello.data),
        deepcopy(hello.axes),
        deepcopy(hello.lookup),
    )
end

function Base.map(f, hello::JuMP.Containers.SparseAxisArray)
    return JuMP.Containers.SparseAxisArray(
        Dict(i => f(hello[i]) for i in keys(hello.data)),
    )
end

function add_variable_as_column(master, column)
    for constr in master.ext[:convexcombination][column.node]
        set_normalized_coefficient(constr, column.var, 1.0)
    end

    for (var, val) in column.coeffs
        if typeof(val) == Float64
            set_normalized_coefficient(
                master.ext[:coverconstraint][column.node][var],
                column.var,
                val,
            )
        else
            for (i, val2) in val
                set_normalized_coefficient(
                    master.ext[:coverconstraint][column.node][var][i],
                    column.var,
                    val2,
                )
            end
        end
    end

    for node in collect(column.node)
        if typeof(node) == Leaf
            set_normalized_coefficient(
                master.ext[:scenprofit_con][node],
                column.var,
                column.obj,
            )
        end
    end
end

# constraint that's used to recreate a discrete subproblem variable in the master problem
function add_mixed_cover(master, sp, column)
    if !(column.node in keys(master.ext[:discrete_con]))
        master.ext[:discrete_var][column.node] = Dict{Int,VariableRef}()
        master.ext[:discrete_con][column.node] = Dict{Int,ConstraintRef}()
        for i in 1:length(sp.ext[:discrete])
            master.ext[:discrete_var][column.node][i] = @variable(master)
            master.ext[:discrete_con][column.node][i] = @constraint(
                master,
                0.0 == master.ext[:discrete_var][column.node][i]
            )
        end
    end

    for i in 1:length(sp.ext[:discrete])
        set_normalized_coefficient(
            master.ext[:discrete_con][column.node][i],
            column.var,
            column.solution[i],
        )
    end
end

function add_column(
    master::JuMP.Model,
    sub_problem::JuMP.Model,
    node::AbstractTree;
    branches = nothing,
)
    coeffs = Dict{Symbol,Any}()

    sol = [0.0]

    if sub_problem.ext[:form] == :mixed
        sol = JuMP.value.(sub_problem.ext[:discrete])
    end

    for (name, variable) in sub_problem.ext[:expansions]
        if typeof(variable) <: AbstractArray
            coeffs[name] = Dict{Any,Float64}()
            vals = JuMP.value.(variable)
            for i in eachindex(variable)
                if sub_problem.ext[:options][name][4] != :Con
                    rval = round(vals[i])
                    if rval != 0.0
                        coeffs[name][i] = rval
                    end
                else
                    coeffs[name][i] = vals[i]
                end
            end
        else
            if sub_problem.ext[:options][name][4] != :Con
                rval = round(JuMP.value(variable))
                if rval != 0.0
                    coeffs[name] = rval
                end
            else
                coeffs[name] = JuMP.value(variable)
            end
        end
    end

    UB = haskey(node.ext, :col_max) ? node.ext[:col_max] : 1.0

    column = Column(
        node,
        coeffs,
        JuMP.value(sub_problem.ext[:objective]),
        JuMP.add_variable(
            master,
            JuMP.build_variable(error, UnitIntervalInformation(UB = UB)),
        ),
        sol,
    )

    if branches != nothing
        for b in branches
            if b.filter != nothing && b.filter(column) == :ban
                set_upper_bound(column.var, 0.0)
                println("New column banned on current branch")
                break
            end
        end
    end

    add_variable_as_column(master, column)
    push!(master.ext[:columns][node], column)
    if sub_problem.ext[:form] == :mixed
        add_mixed_cover(master, sub_problem, column)
    end
    return column
end

"""
	solve(judge::JuDGEModel;
	      termination::Termination=Termination(),
	      max_no_int::Int=typemax(Int),
	      blocks::Union{Nothing,Array{Array{AbstractTree,1},1}}=nothing,
	      warm_starts::Bool=false,
	      optimizer_attributes::Union{Nothing,Function}=nothing,
	      mp_callback::Union{Nothing,Function}=nothing,
	      prune::Float64=Inf,
	      heuristic::Union{Nothing,Function}=nothing,
	      verbose::Int=2)

Solve a JuDGEModel `judge` without branch-and-price.

### Required Arguments
`judge` is the JuDGE model that we wish to solve.

### Optional Arguments

`termination` is a `Termination` object containing all the stopping conditions.

`max_no_int` is the maximum number of iterations yielding a fractional solution before a MIP solve is
performed on the master. By default, the MIP solves will not occur until the relaxed bound gap is less
than the `relgap` / `absgap` stopping conditions. To override this, set `max_no_int` to the negative
of the number desired value.

`blocks` specifies the groups of nodes to solve in each iteration (these groups can be generated using
`JuDGE.get_groups()`, or created manually), after all nodes have been solved, a full pricing iteration
is used to compute an updated lower bound. See `advanced.jl` for more details.

`warm_starts` boolean specifing whether to use warm starts for subproblems and binary solves of master
problem.

`optimizer_attributes` can be set to a specific function that dynamically changes optimizer attributes
for the subproblems; this should only be used by people who have examined the `advanced.jl` example.

`mp_callback` is a user-defined function that specifies termination conditions for MIP solves of the
master problem. See examples/advanced.jl.

`prune` is used to stop the algorithm before convergence, if a known upper bound for the problem is
specified.

`heuristic` is a user-defined function that typically would perform an improvement heuristic when a
new incumbent is found.

`verbose` if 0, all output from solve will be suppressed, if 1, the subproblem solve process will be
suppressed. Default is 2.

### Examples
    JuDGE.solve(jmodel, termination=Termination(rlx_abstol=10^-6))
	JuDGE.solve(jmodel, termination=Termination(rlx_abstol=10^-6), max_no_int=-5)
"""
function solve(
    judge::JuDGEModel;
    termination::Termination = Termination(),
    max_no_int::Int = typemax(Int),
    blocks::Union{Nothing,Array{Array{AbstractTree,1},1}} = nothing,
    warm_starts::Bool = false,
    optimizer_attributes::Union{Nothing,Function} = nothing,
    mp_callback::Union{Nothing,Function} = nothing,
    prune::Float64 = Inf,
    heuristic::Union{Nothing,Function} = nothing,
    verbose::Int = 2,
)
    current = InitialConvergenceState()
    empty!(judge.log)
    push!(judge.log, current)
    if verbose > 0
        print("")
        print("")
        println(
            "Relaxed ObjVal  |   Upper Bound   Lower Bound  |  Absolute Diff   Relative Diff  |  Fractional  |      Time     Iter",
        )
    end
    # set up times for use in convergence
    initial_time = time()
    obj = Inf
    heuristicobj = Inf
    nodes = collect(judge.tree)
    objduals = Dict{AbstractTree,Float64}()
    redcosts = Dict{AbstractTree,Float64}()
    for node in nodes
        objduals[node] = -1
        redcosts[node] = -1
    end
    no_int_count = 0
    if optimizer_attributes != nothing
        optimizer_attributes(judge, false, true)
    end

    max_char = length(nodes[end].name) + length(string(length(nodes)))
    function get_whitespace(name::String, number::Int)
        blank = "  "
        spaces = max_char - length(name) - length(string(number))
        for i in 1:spaces
            blank *= " "
        end
        return blank
    end
    if blocks == nothing
        blocks = [nodes]
    end

    block = -1

    if judge.master_problem.ext[:mip]
        remove_binary(judge)
    end
    set_banned_variables!(judge)

    remove_branch_constraints!(judge)
    b_con = add_branch_constraints!(judge)

    optimize!(judge.master_problem)
    status = termination_status(judge.master_problem)
    if status == MOI.NUMERICAL_ERROR
        println("\nMaster problem returned a MOI.NUMERICAL_ERROR")
        return
    end

    while true
        if block <= 0
            nodes2 = nodes
        else
            nodes2 = blocks[block]
        end
        # perform the main iterations
        for i in 1:length(nodes2)
            node = nodes2[i]
            sp = judge.sub_problems[node]
            updateduals(judge.master_problem, sp, node, status, current.iter)
            if verbose == 2
                overprint(
                    "Solving subproblem for node " *
                    node.name *
                    get_whitespace(node.name, i) *
                    string(i) *
                    "/" *
                    string(length(nodes2)),
                )
            end

            optimize!(sp)

            if termination_status(sp) != MOI.OPTIMAL &&
               termination_status(sp) != MOI.INTERRUPTED &&
               termination_status(sp) != MOI.TIME_LIMIT
                @warn(
                    "Solve for subproblem " *
                    node.name *
                    " exited with status " *
                    string(termination_status(sp))
                )
                println("Subproblem is infeasible or unbounded")
                for (sp2, con) in b_con
                    delete(sp2, con)
                end
                return :sp_infeasible
            end
            if warm_starts
                vars = all_variables(sp)
                set_start_value.(vars, JuMP.value.(vars))
            end
        end
        if verbose == 2
            overprint("")
        end
        frac = 0
        if status != MOI.INFEASIBLE_OR_UNBOUNDED &&
           status != MOI.INFEASIBLE &&
           status != MOI.DUAL_INFEASIBLE
            for node in nodes2
                objduals[node] = objective_bound(judge.sub_problems[node])
                redcosts[node] = objective_value(judge.sub_problems[node])
            end
            if status != MOI.OPTIMAL
                @warn(
                    "Master problem did not solve to optimality: " *
                    string(status)
                )
            elseif block == 0
                getlowerbound(judge, objduals)
            end
        elseif judge.bounds.LB > -Inf
            println("\nMaster problem is infeasible or unbounded")
            for (sp, con) in b_con
                delete(sp, con)
            end
            judge.bounds.LB = Inf
            judge.bounds.UB = Inf
            return
        end

        num_var = num_variables(judge.master_problem)
        for node in nodes2
            if redcosts[node] < -10^-10 ||
               status == MOI.INFEASIBLE_OR_UNBOUNDED ||
               status == MOI.INFEASIBLE ||
               status == MOI.DUAL_INFEASIBLE
                column = add_column(
                    judge.master_problem,
                    judge.sub_problems[node],
                    node,
                )
                if warm_starts
                    set_start_value(column.var, 0.0)
                end
            end
        end

        optimize!(judge.master_problem)
        status = termination_status(judge.master_problem)
        if status == MOI.NUMERICAL_ERROR
            println("\nMaster problem returned a MOI.NUMERICAL_ERROR")
            for (sp, con) in b_con
                delete(sp, con)
            end
            return
        elseif status != MOI.INFEASIBLE_OR_UNBOUNDED &&
               status != MOI.INFEASIBLE &&
               status != MOI.DUAL_INFEASIBLE
            frac = fractionalcount(judge, termination.inttol)
            obj = objective_value(judge.master_problem)
            if frac == 0
                if heuristic != nothing && heuristicobj > obj
                    if heuristic(judge) < 0.0 &&
                       termination.allow_frac ∉
                       [:first_fractional, :no_binary_solve]
                        solve_master_binary(
                            judge,
                            initial_time,
                            termination,
                            warm_starts,
                            nothing,
                            verbose,
                        )
                    else
                        optimize!(judge.master_problem)
                    end
                    obj = objective_value(judge.master_problem)
                    heuristicobj = obj
                end
                judge.bounds.UB = obj

                no_int_count = 0
            else
                no_int_count += 1
            end
        end
        current = ConvergenceState(
            obj,
            judge.bounds.UB,
            judge.bounds.LB,
            time() - initial_time,
            current.iter + 1,
            frac,
        )
        if verbose > 0
            display(current)
        end
        push!(judge.log, current)
        if prune < judge.bounds.LB
            if verbose > 0
                println("\nDominated by incumbent.")
            end
            for (sp, con) in b_con
                delete(sp, con)
            end
            return
        elseif has_converged(termination, current)
            if frac > 0
                solve_master_binary(
                    judge,
                    initial_time,
                    termination,
                    warm_starts,
                    nothing,
                    verbose,
                )
            end
            if heuristic != nothing
                if heuristic(judge) < 0.0 &&
                   termination.allow_frac ∉
                   [:first_fractional, :no_binary_solve]
                    solve_master_binary(
                        judge,
                        initial_time,
                        termination,
                        warm_starts,
                        nothing,
                        verbose,
                    )
                    heuristicobj = judge.bounds.UB
                else
                    optimize!(judge.master_problem)
                end
            end
            if verbose > 0
                println("\nConvergence criteria met.")
            end
            for (sp, con) in b_con
                delete(sp, con)
            end
            return
        elseif termination.allow_frac == :first_fractional && frac > 0
            if verbose > 0
                println("\nFractional solution found.")
            end
            for (sp, con) in b_con
                delete(sp, con)
            end
            return
        elseif (
            (
                max_no_int > 0 &&
                no_int_count >= max_no_int &&
                current.num_frac > 0
            ) && (
                current.rlx_abs < termination.abstol ||
                current.rlx_rel < termination.reltol
            )
        ) || (
            max_no_int < 0 &&
            no_int_count >= -max_no_int &&
            current.num_frac > 0
        )
            current = solve_master_binary(
                judge,
                initial_time,
                termination,
                warm_starts,
                mp_callback,
                verbose,
            )
            if heuristic != nothing && heuristicobj > judge.bounds.UB
                if heuristic(judge) < 0.0
                    current = solve_master_binary(
                        judge,
                        initial_time,
                        termination,
                        warm_starts,
                        mp_callback,
                        verbose,
                    )
                else
                    optimize!(judge.master_problem)
                end
                heuristicobj = judge.bounds.UB
            end

            if has_converged(termination, current)
                if verbose > 0
                    println("\nConvergence criteria met.")
                end
                for (sp, con) in b_con
                    delete(sp, con)
                end
                return
            elseif termination.allow_frac == :binary_solve
                remove_binary(judge)
                optimize!(judge.master_problem)
            end
            no_int_count = 0
        end

        if optimizer_attributes == nothing
            if block == 0 && num_var == num_variables(judge.master_problem)
                solve_master_binary(
                    judge,
                    initial_time,
                    termination,
                    warm_starts,
                    nothing,
                    verbose,
                )
                if verbose > 0
                    println("\nStalled: exiting.")
                end
                for (sp, con) in b_con
                    delete(sp, con)
                end
                return
            end
        elseif optimizer_attributes(
            judge,
            num_var == num_variables(judge.master_problem),
            false,
        )
            if verbose > 0
                println("\nStalled: exiting.")
            end
            for (sp, con) in b_con
                delete(sp, con)
            end
            return
        end
        if length(blocks) == 1
            block = 0
        else
            block = (block + 1) % (length(blocks) + 1)
        end
    end
    for (sp, con) in b_con
        delete(sp, con)
    end
end

function getlowerbound(judge::JuDGEModel, objduals::Dict{AbstractTree,Float64})
    lb = objective_value(judge.master_problem)
    for (i, objdual) in objduals
        if :sum_max in keys(i.ext)
            lb += objdual * i.ext[:sum_max]
        else
            lb += objdual
        end
    end
    if lb > judge.bounds.LB
        judge.bounds.LB = lb
    end
end

function fractionalcount(jmodel::JuDGEModel, inttol::Float64)
    count = 0

    for (node, sp) in jmodel.sub_problems
        if sp.ext[:form] == :binarycolumns
            for col in jmodel.master_problem.ext[:columns][node]
                val = JuMP.value(col.var)
                count += min(val - floor(val), ceil(val) - val) > inttol ? 1 : 0
            end
            return count
        end
    end

    for node in collect(jmodel.tree)
        for x in keys(jmodel.master_problem.ext[:expansions][node])
            if jmodel.sub_problems[jmodel.tree].ext[:options][x][4] != :Con
                var = jmodel.master_problem.ext[:expansions][node][x]
                if typeof(var) <: AbstractArray
                    for key in eachindex(var)
                        val = JuMP.value(var[key])
                        count +=
                            min(val - floor(val), ceil(val) - val) > inttol ?
                            1 : 0
                    end
                else
                    val = JuMP.value(var)
                    count +=
                        min(val - floor(val), ceil(val) - val) > inttol ? 1 : 0
                end
            end
        end
        if jmodel.sub_problems[node].ext[:form] == :mixed
            for (i, var) in jmodel.master_problem.ext[:discrete_var][node]
                val = JuMP.value(var)
                count += min(val - floor(val), ceil(val) - val) > inttol ? 1 : 0
            end
        end
    end

    return count
end

function updateduals(master, sub_problem, node, status, iter)
    if status != MOI.INFEASIBLE_OR_UNBOUNDED &&
       status != MOI.INFEASIBLE &&
       status != MOI.DUAL_INFEASIBLE
        for (name, var) in sub_problem.ext[:expansions]
            if typeof(var) <: AbstractArray
                for i in eachindex(var)
                    set_objective_coefficient(
                        sub_problem,
                        var[i],
                        -dual(master.ext[:coverconstraint][node][name][i]),
                    )
                end
            else
                set_objective_coefficient(
                    sub_problem,
                    var,
                    -dual(master.ext[:coverconstraint][node][name]),
                )
            end
        end
        nodes = collect(node)
        total = 0.0
        for n in nodes
            if typeof(n) == Leaf
                total -= dual(master.ext[:scenprofit_con][n])
            end
        end
        set_objective_coefficient(
            sub_problem,
            sub_problem.ext[:objective],
            total,
        )

        cc_sum = 0.0
        for constr in master.ext[:convexcombination][node]
            set = typeof(constraint_object(constr).set)
            if set <: MOI.LessThan
                cc_sum += dual(constr)
            elseif set <: MOI.GreaterThan
                cc_sum -= dual(constr)
            elseif set <: MOI.EqualTo
                cc_sum += dual(constr)
            end
        end

        set_objective_function(
            sub_problem,
            objective_function(sub_problem) -
            objective_function(sub_problem).constant - cc_sum,
        )
    else
        if iter % 2 == 0
            oc = -10.0^10
        else
            oc = 10.0^10
        end
        for (name, var) in sub_problem.ext[:expansions]
            if typeof(var) <: AbstractArray
                for i in eachindex(var)
                    set_objective_coefficient(sub_problem, var[i], oc)
                end
            else
                set_objective_coefficient(sub_problem, var, oc)
            end
        end
        set_objective_coefficient(sub_problem, sub_problem.ext[:objective], 1.0)
    end
end

function get_objval(jmodel::JuDGEModel; risk = jmodel.risk)
    scenario_objs = Dict{Leaf,Float64}()

    for (leaf, var) in jmodel.master_problem.ext[:scenprofit_var]
        scenario_objs[leaf] = JuMP.value(var)
    end

    return compute_objval(scenario_objs, jmodel.probabilities, risk)
end

"""
	resolve_subproblems(judge::JuDGEModel)

Once a JuDGE model has converged, it is necessary to re-solve the subproblems to find the optimal decisions within each node.

### Required Arguments
`jmodel` is the JuDGE model that we wish to solve.

### Examples
    resolve_subproblems(judge)
"""
function resolve_subproblems(jmodel::JuDGEModel)
    fix_expansions(jmodel)
    return resolve_fixed(jmodel)
end

function fix_expansions(jmodel::JuDGEModel)
    if termination_status(jmodel.master_problem) != MOI.OPTIMAL &&
       termination_status(jmodel.master_problem) != MOI.INTERRUPTED &&
       termination_status(jmodel.master_problem) != MOI.LOCALLY_SOLVED &&
       termination_status(jmodel.master_problem) != MOI.INTERRUPTED
        error("You need to first solve the decomposed model.")
    end

    for node in collect(jmodel.tree)
        sp = jmodel.sub_problems[node]
        set_objective_function(
            sp,
            objective_function(sp) - objective_function(sp).constant,
        )
        set_objective_coefficient(sp, sp.ext[:objective], 1.0)
        for (name, var) in jmodel.master_problem.ext[:expansions][node]
            var2 = sp.ext[:expansions][name]
            if typeof(var) <: AbstractArray
                for i in eachindex(var)
                    value = 0.0
                    if sp.ext[:options][name][1] == :state
                        prev = node.parent
                        if prev == nothing
                            var3 =
                                jmodel.master_problem.ext[:expansions][node][name][i]
                            value = JuMP.value(var3) - sp.ext[:options][name][7]
                        else
                            var3 =
                                jmodel.master_problem.ext[:expansions][node][name][i]
                            var4 =
                                jmodel.master_problem.ext[:expansions][prev][name][i]
                            value = JuMP.value(var3) - JuMP.value(var4)
                        end
                    else
                        con_obj = constraint_object(
                            jmodel.master_problem.ext[:coverconstraint][node][name][i],
                        )
                        for prev in history(node)
                            var3 =
                                jmodel.master_problem.ext[:expansions][prev][name][i]
                            if var3 in keys(con_obj.func.terms)
                                value +=
                                    JuMP.value(var3) * -con_obj.func.terms[var3]
                            end
                        end
                    end
                    if sp.ext[:options][name][1] == :shutdown
                        JuMP.set_lower_bound(var2[i], value)
                    elseif sp.ext[:options][name][1] == :expansion
                        JuMP.set_upper_bound(var2[i], value)
                    elseif sp.ext[:options][name][1] == :enforced ||
                           sp.ext[:options][name][1] == :state
                        JuMP.fix(var2[i], value, force = true)
                    end
                    set_objective_coefficient(sp, var2[i], 0.0)
                end
            elseif isa(var, VariableRef)
                value = 0.0
                if sp.ext[:options][name][1] == :state
                    prev = node.parent
                    if prev == nothing
                        var3 =
                            jmodel.master_problem.ext[:expansions][node][name]
                        value = JuMP.value(var3) - sp.ext[:options][name][7]
                    else
                        var3 =
                            jmodel.master_problem.ext[:expansions][node][name]
                        var4 =
                            jmodel.master_problem.ext[:expansions][prev][name]
                        value = JuMP.value(var3) - JuMP.value(var4)
                    end
                else
                    con_obj = constraint_object(
                        jmodel.master_problem.ext[:coverconstraint][node][name],
                    )
                    for prev in history(node)
                        var3 =
                            jmodel.master_problem.ext[:expansions][prev][name]
                        if var3 in keys(con_obj.func.terms)
                            value +=
                                JuMP.value(var3) * -con_obj.func.terms[var3]
                        end
                    end
                end
                if sp.ext[:options][name][1] == :shutdown
                    JuMP.set_lower_bound(var2, value)
                elseif sp.ext[:options][name][1] == :expansion
                    JuMP.set_upper_bound(var2, value)
                elseif sp.ext[:options][name][1] == :enforced ||
                       sp.ext[:options][name][1] == :state
                    JuMP.fix(var2, value, force = true)
                end
                set_objective_coefficient(sp, var2, 0.0)
            end
        end
    end
end

function resolve_fixed(jmodel::JuDGEModel)
    for n in collect(jmodel.tree)
        JuMP.optimize!(jmodel.sub_problems[n])
    end

    scenario_objs = Dict{Leaf,Float64}()
    for (leaf, con) in jmodel.master_problem.ext[:scenprofit_con]
        obj = 0.0
        for node in history(leaf)
            obj += objective_value(jmodel.sub_problems[node])
            for key in keys(jmodel.master_problem.ext[:expansions][node])
                var = jmodel.master_problem.ext[:expansions][node][key]
                if isa(var, VariableRef)
                    obj += JuMP.value(var) * normalized_coefficient(con, var)
                elseif typeof(var) <: AbstractArray
                    for v in eachindex(var)
                        obj +=
                            JuMP.value(var[v]) *
                            normalized_coefficient(con, var[v])
                    end
                end
            end
        end
        scenario_objs[leaf] = obj
    end

    return compute_objval(scenario_objs, jmodel.probabilities, jmodel.risk)
end

include("model_verification.jl")

# pretty printing
function Base.show(io::IO, ::MIME"text/plain", judge::JuDGEModel)
    print(io, "JuDGE Model with:\n")
    println(io, "  Tree: ", judge.tree)
    print(io, "  Expansion variables: ")
    keys = judge.sub_problems[judge.tree].ext[:expansions]
    for (i, ii) in keys
        print(io, "$(i) ")
    end
end

function Base.show(io::IO, judge::JuDGEModel)
    print(io, "JuDGE Model with:\n")
    println(io, "  Tree: ", judge.tree)
    print(io, "  Expansion variables: ")
    keys = judge.sub_problems[judge.tree].ext[:expansions]
    for (i, ii) in keys
        print(io, "$(i) ")
    end
end

include("output.jl")

export @expansion,
    @shutdown,
    @enforced,
    @state,
    @expansionconstraint,
    @capitalcosts,
    @ongoingcosts,
    JuDGEModel,
    Risk,
    RiskNeutral,
    Leaf,
    Tree,
    AbstractTree,
    narytree,
    ConditionallyUniformProbabilities,
    UniformLeafProbabilities,
    get_node,
    tree_from_leaves,
    tree_from_nodes,
    tree_from_file,
    DetEqModel,
    resolve_subproblems,
    Termination

end
