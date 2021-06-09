function build_master(
    sub_problems::T where {T<:Dict},
    tree::T where {T<:AbstractTree},
    probabilities::Dict{AbstractTree,Float64},
    solver,
    discount_factor::Float64,
    risk::Union{Risk,Array{Risk,1}},
    sideconstraints,
)
    if typeof(solver) <: Tuple
        model = Model(solver[1])
    else
        model = Model(solver)
    end
    @objective(model, Min, 0)

    model.ext[:columns] = Dict{AbstractTree,Array{Column,1}}()
    for node in collect(tree)
        model.ext[:columns][node] = Array{Column,1}()
    end
    leafs = Array{Leaf,1}()

    model.ext[:expansions] = Dict{AbstractTree,Dict{Symbol,Any}}()

    model.ext[:options] = Dict{Symbol,Tuple}()
    for (name, options) in sub_problems[tree].ext[:options]
        model.ext[:options][name] = options
    end

    form = nothing
    for (name, options) in model.ext[:options]
        if options[4] == :Bin
            if form == :continuous
                form = :mixed
                break
            else
                form = :binary
            end
        elseif options[4] == :Int
            form = :mixed
            break
        else
            if form == :binary
                form = :mixed
                break
            else
                form = :continuous
            end
        end
    end

    model.ext[:discrete_con] = Dict{AbstractTree,Dict{Int,ConstraintRef}}()
    model.ext[:discrete_var] = Dict{AbstractTree,Dict{Int,VariableRef}}()

    for (node, sp) in sub_problems
        if haskey(node.ext, :sum_max) ||
           haskey(node.ext, :sum_min) ||
           haskey(node.ext, :sum_set)
            sp.ext[:form] = :binarycolumns
        else
            sp.ext[:form] = form
        end
        if sp.ext[:form] ∉ [:binary, :binarycolumns]
            vars = all_variables(sp)
            if form == :continuous
                for var in vars
                    if is_binary(var) || is_integer(var)
                        sp.ext[:form] = :mixed
                        break
                    end
                end
            end

            if sp.ext[:form] == :mixed
                sp.ext[:discrete] = VariableRef[]
                sp.ext[:discrete_branch] = VariableRef[]
                ignore = VariableRef[]

                for (name, variable) in sp.ext[:expansions]
                    if typeof(variable) <: AbstractArray
                        for index in eachindex(variable)
                            push!(ignore, variable[index])
                        end
                    else
                        push!(ignore, variable)
                    end
                end

                for var in vars
                    if (is_binary(var) || is_integer(var)) && var ∉ ignore
                        push!(sp.ext[:discrete], var)
                        push!(sp.ext[:discrete_branch], var)
                    end
                end
            end
        end
    end

    for (node, sp) in sub_problems
        if typeof(node) == Leaf
            push!(leafs, node)
        end
        model.ext[:expansions][node] = Dict{Symbol,Any}()

        for (name, variable) in sp.ext[:expansions]
            if sp.ext[:options][name][4] == :Bin
                model.ext[:expansions][node][name] =
                    copy_variable!(model, variable, relaxbinary)
            elseif sp.ext[:options][name][4] == :Int
                model.ext[:expansions][node][name] =
                    copy_variable!(model, variable, relaxinteger)
            else
                model.ext[:expansions][node][name] =
                    copy_variable!(model, variable)
            end
            if model.ext[:options][name][5] != nothing
                if typeof(variable) <: AbstractArray
                    for i in eachindex(variable)
                        set_lower_bound(
                            model.ext[:expansions][node][name][i],
                            model.ext[:options][name][5],
                        )
                    end
                else
                    set_lower_bound(
                        model.ext[:expansions][node][name],
                        model.ext[:options][name][5],
                    )
                end
            end
            if model.ext[:options][name][6] != nothing
                if typeof(variable) <: AbstractArray
                    for i in eachindex(variable)
                        set_upper_bound(
                            model.ext[:expansions][node][name][i],
                            model.ext[:options][name][6],
                        )
                    end
                else
                    set_upper_bound(
                        model.ext[:expansions][node][name],
                        model.ext[:options][name][6],
                    )
                end
            end
        end
    end

    model.ext[:scenprofit_var] = Dict{Leaf,VariableRef}()
    model.ext[:scenprofit_con] = Dict{Leaf,ConstraintRef}()

    for leaf in leafs
        model.ext[:scenprofit_var][leaf] = @variable(model)
        model.ext[:scenprofit_con][leaf] =
            @constraint(model, 0 == model.ext[:scenprofit_var][leaf])

        set_objective_coefficient(
            model,
            model.ext[:scenprofit_var][leaf],
            probabilities[leaf],
        )
    end

    risk_objectives = AffExpr[]
    remain = 1.0
    if typeof(risk) == Risk
        if (risk.α == 1.0 || risk.λ == 0.0) && risk.bound == nothing
            risk = []
        else
            risk = [risk]
        end
    end
    for i in 1:length(risk)
        risk_objective = AffExpr(0.0)
        eta = @variable(model)
        for leaf in leafs
            v = @variable(model)
            w = @variable(model)
            set_lower_bound(v, 0.0)
            set_lower_bound(w, 0.0)
            if risk[i].offset == nothing
                @constraint(model, v >= eta - model.ext[:scenprofit_var][leaf])
                @constraint(model, w >= model.ext[:scenprofit_var][leaf] - eta)
                add_to_expression!(
                    risk_objective,
                    (model.ext[:scenprofit_var][leaf]) * probabilities[leaf],
                )
            else
                @constraint(
                    model,
                    v >=
                    eta - model.ext[:scenprofit_var][leaf] +
                    risk[i].offset[leaf]
                )
                @constraint(
                    model,
                    w >=
                    model.ext[:scenprofit_var][leaf] - risk[i].offset[leaf] -
                    eta
                )
                add_to_expression!(
                    risk_objective,
                    (model.ext[:scenprofit_var][leaf] - risk[i].offset[leaf]) *
                    probabilities[leaf],
                )
            end
            add_to_expression!(
                risk_objective,
                probabilities[leaf] * (v + w / risk[i].α * (1 - risk[i].α)),
            )
        end
        remain -= risk[i].λ
        push!(risk_objectives, risk_objective)
    end

    for leaf in leafs
        nodes = history(leaf)
        for n in eachindex(nodes)
            node = nodes[n]
            sp = sub_problems[node]
            df = discount_factor^depth(node)
            for (name, variable) in sp.ext[:expansions]
                interval =
                    max(
                        1,
                        n - sp.ext[:options][name][3] -
                        sp.ext[:options][name][2] + 1,
                    ):n-sp.ext[:options][name][2]
                disc = Dict{Int,Float64}()
                for i in interval
                    disc[i] = df / discount_factor^(i - 1)
                end
                if typeof(variable) <: AbstractArray
                    for i in eachindex(variable)
                        cost_coef =
                            df * coef(sp.ext[:capitalcosts], variable[i])
                        for j in interval
                            cost_coef +=
                                disc[j] *
                                coef(sp.ext[:ongoingcosts], variable[i])
                        end
                        set_normalized_coefficient(
                            model.ext[:scenprofit_con][leaf],
                            model.ext[:expansions][node][name][i],
                            cost_coef,
                        )
                    end
                else
                    cost_coef = df * coef(sp.ext[:capitalcosts], variable)
                    for j in interval
                        cost_coef +=
                            disc[j] * coef(sp.ext[:ongoingcosts], variable)
                    end
                    set_normalized_coefficient(
                        model.ext[:scenprofit_con][leaf],
                        model.ext[:expansions][node][name],
                        cost_coef,
                    )
                end
            end
        end
    end

    # create the cover constraints
    model.ext[:coverconstraint] = Dict{AbstractTree,Dict{Symbol,Any}}()
    for (node, sp) in sub_problems
        model.ext[:coverconstraint][node] = Dict{Symbol,Any}()
        past = history(node)
        for (name, variable) in sp.ext[:expansions]
            interval =
                sp.ext[:options][name][2]+1:min(
                    sp.ext[:options][name][2] + sp.ext[:options][name][3],
                    length(past),
                )
            if typeof(variable) <: AbstractArray
                model.ext[:coverconstraint][node][name] = Dict()
                for i in eachindex(variable)
                    if sp.ext[:options][name][1] == :shutdown
                        model.ext[:coverconstraint][node][name][i] =
                            @constraint(
                                model,
                                0 >= sum(
                                    model.ext[:expansions][past[index]][name][i]
                                    for index in interval
                                )
                            )
                    elseif sp.ext[:options][name][1] == :expansion
                        model.ext[:coverconstraint][node][name][i] =
                            @constraint(
                                model,
                                0 <= sum(
                                    model.ext[:expansions][past[index]][name][i]
                                    for index in interval
                                )
                            )
                    elseif sp.ext[:options][name][1] == :enforced
                        model.ext[:coverconstraint][node][name][i] =
                            @constraint(
                                model,
                                0 == sum(
                                    model.ext[:expansions][past[index]][name][i]
                                    for index in interval
                                )
                            )
                    elseif sp.ext[:options][name][1] == :state
                        if node.parent == nothing
                            model.ext[:coverconstraint][node][name][i] =
                                @constraint(
                                    model,
                                    0 ==
                                    -sp.ext[:options][name][7] +
                                    model.ext[:expansions][node][name][i]
                                )
                        else
                            model.ext[:coverconstraint][node][name][i] =
                                @constraint(
                                    model,
                                    0 ==
                                    model.ext[:expansions][node][name][i] -
                                    model.ext[:expansions][node.parent][name][i]
                                )
                        end
                    end
                end
                # if typeof(node)==Leaf && !sp.ext[:options][name][4]
                # 	for i in eachindex(variable)
                # 		@constraint(model, sum(model.ext[:expansions][past][name][i] for past in history_function(node))<=1)
                # 	end
                # end
            else
                if sp.ext[:options][name][1] == :shutdown
                    model.ext[:coverconstraint][node][name] = @constraint(
                        model,
                        0 >= sum(
                            model.ext[:expansions][past[index]][name] for
                            index in interval
                        )
                    )
                elseif sp.ext[:options][name][1] == :expansion
                    model.ext[:coverconstraint][node][name] = @constraint(
                        model,
                        0 <= sum(
                            model.ext[:expansions][past[index]][name] for
                            index in interval
                        )
                    )
                elseif sp.ext[:options][name][1] == :enforced
                    model.ext[:coverconstraint][node][name] = @constraint(
                        model,
                        0 == sum(
                            model.ext[:expansions][past[index]][name] for
                            index in interval
                        )
                    )
                elseif sp.ext[:options][name][1] == :state
                    if node.parent == nothing
                        model.ext[:coverconstraint][node][name] = @constraint(
                            model,
                            0 ==
                            -sp.ext[:options][name][7] +
                            model.ext[:expansions][node][name]
                        )
                    else
                        model.ext[:coverconstraint][node][name] = @constraint(
                            model,
                            0 ==
                            model.ext[:expansions][node][name] -
                            model.ext[:expansions][node.parent][name]
                        )
                    end
                end
                # if typeof(node)==Leaf && !sp.ext[:options][name][4]
                # 	@constraint(model, sum(model.ext[:expansions][past][name] for past in history_function(node))<=1)
                # end
            end
        end
    end

    model.ext[:convexcombination] = Dict{AbstractTree,Array{ConstraintRef,1}}()
    model.ext[:mip] = false
    model.ext[:branch_cons] = ConstraintRef[]

    for node in keys(sub_problems)
        model.ext[:convexcombination][node] = ConstraintRef[]
        if haskey(node.ext, :sum_max)
            push!(
                model.ext[:convexcombination][node],
                @constraint(model, 0 <= node.ext[:sum_max])
            )
        end
        if haskey(node.ext, :sum_min)
            push!(
                model.ext[:convexcombination][node],
                @constraint(model, 0 >= node.ext[:sum_min])
            )
        end
        if haskey(node.ext, :sum_set)
            if length(model.ext[:convexcombination][node]) > 0
                error(":eq and :le/:ge specifed for node " * node.name)
            end
            push!(
                model.ext[:convexcombination][node],
                @constraint(model, 0 == node.ext[:sum_set])
            )
        end
        if length(model.ext[:convexcombination][node]) == 0
            push!(
                model.ext[:convexcombination][node],
                @constraint(model, 0 == 1)
            )
        end
    end

    if typeof(sideconstraints) <: Function
        map(Main.eval, unpack_expansions(model.ext[:expansions])) #bring expansion variables into global scope
        sideconstraints(model, tree)
        map(Main.eval, clear_expansions(model.ext[:expansions]))
    end

    objective_fn = objective_function(model) * remain

    for i in 1:length(risk_objectives)
        objective_fn += risk_objectives[i] * risk[i].λ
        if risk[i].bound != nothing
            surplus = @variable(model)
            set_lower_bound(surplus, 0)
            @constraint(model, risk_objectives[i] <= risk[i].bound + surplus)
            objective_fn += risk[i].penalty * surplus
        end
    end

    set_objective_function(model, objective_fn)

    return model
end

function solve_master_binary(
    judge::JuDGEModel,
    initial_time::Float64,
    termination::Termination,
    warm_starts::Bool,
    mp_callback,
    verbose,
)
    set = 0
    current = nothing
    if (
        termination.allow_frac == :binary_solve ||
        termination.allow_frac == :binary_solve_return_relaxation
    ) #&& judge.log[end].int>0
        if verbose == 2
            print("Solving master problem as MIP")
        end
        set = 1
        if solve_binary(
            judge,
            termination.abstol,
            termination.reltol,
            warm_starts,
            mp_callback,
        )
            current = ConvergenceState(
                judge.log[end].obj,
                judge.bounds.UB,
                judge.bounds.LB,
                time() - initial_time,
                judge.log[end].iter + 1,
                0,
            )
            if verbose > 0
                if verbose == 2
                    overprint("")
                end
                display(current, relaxation = false)
            end
            push!(judge.log, current)
        else
            if verbose == 2
                overprint("")
            end
        end
    end
    if termination.allow_frac == :binary_solve_return_relaxation && set == 1
        remove_binary(judge)
        optimize!(judge.master_problem)
    end
    return current
end

function solve_binary(
    judge::JuDGEModel,
    abstol::Float64,
    reltol::Float64,
    warm_starts::Bool,
    mp_callback,
)
    if typeof(judge.master_solver) <: Tuple
        JuMP.set_optimizer(judge.master_problem, judge.master_solver[2])
    end

    judge.master_problem.ext[:mip] = true

    for node in keys(judge.master_problem.ext[:expansions])
        if judge.sub_problems[node].ext[:form] == :binarycolumns
            for col in judge.master_problem.ext[:columns][node]
                if !haskey(col.node.ext, :col_max)
                    set_binary(col.var)
                else
                    set_integer(col.var)
                end
            end
        end

        for x in keys(judge.master_problem.ext[:expansions][node])
            var = judge.master_problem.ext[:expansions][node][x]
            if judge.sub_problems[judge.tree].ext[:options][x][4] == :Bin
                if typeof(var) <: JuMP.Containers.SparseAxisArray
                    for key in keys(var.data)
                        set_binary(var.data[key])
                    end
                elseif typeof(var) <: AbstractArray
                    for key in keys(var)
                        set_binary(var[key])
                    end
                else
                    set_binary(var)
                end
            elseif judge.sub_problems[judge.tree].ext[:options][x][4] == :Int
                if typeof(var) <: JuMP.Containers.SparseAxisArray
                    for key in keys(var.data)
                        set_integer(var.data[key])
                    end
                elseif typeof(var) <: AbstractArray
                    for key in keys(var)
                        set_integer(var[key])
                    end
                else
                    set_integer(var)
                end
            end
        end
        if judge.sub_problems[node].ext[:form] == :mixed
            for (i, var) in judge.master_problem.ext[:discrete_var][node]
                set_integer(var)
            end
        end
    end

    if typeof(mp_callback) <: Function
        mp_callback(judge, abstol, reltol)
        printright("Initialising")
    end
    optimize!(judge.master_problem)
    status = termination_status(judge.master_problem)
    if status != MOI.INFEASIBLE_OR_UNBOUNDED &&
       status != MOI.INFEASIBLE &&
       status != MOI.DUAL_INFEASIBLE
        obj = objective_value(judge.master_problem)
        if obj < judge.bounds.UB
            judge.bounds.UB = obj
            if warm_starts
                vars = all_variables(judge.master_problem)
                set_start_value.(vars, JuMP.value.(vars))
            end
        end
        return true
    else
        return false
    end
end

function remove_binary(judge::JuDGEModel)
    # if judge.master_problem.ext[:binarycolumns]
    # 	for col in judge.master_problem.ext[:columns]
    # 		if !haskey(col.node.ext,:col_max)
    # 			unset_binary(col.var)
    # 		else
    # 			unset_integer(col.var)
    # 		end
    # 	end
    # end
    for node in keys(judge.master_problem.ext[:expansions])
        if judge.sub_problems[node].ext[:form] == :binarycolumns
            for col in judge.master_problem.ext[:columns][node]
                if !haskey(col.node.ext, :col_max)
                    unset_binary(col.var)
                else
                    unset_integer(col.var)
                end
            end
        end
        for x in keys(judge.master_problem.ext[:expansions][node])
            var = judge.master_problem.ext[:expansions][node][x]
            if judge.sub_problems[judge.tree].ext[:options][x][4] == :Bin
                if typeof(var) <: AbstractArray
                    for key in eachindex(var)
                        unset_binary(var[key])
                    end
                else
                    unset_binary(var)
                end
            elseif judge.sub_problems[judge.tree].ext[:options][x][4] == :Int
                if typeof(var) <: AbstractArray
                    for key in eachindex(var)
                        unset_integer(var[key])
                    end
                else
                    unset_integer(var)
                end
            end
        end
        if judge.sub_problems[node].ext[:form] == :mixed
            for (i, var) in judge.master_problem.ext[:discrete_var][node]
                unset_integer(var)
            end
        end
    end

    judge.master_problem.ext[:mip] = false

    if typeof(judge.master_solver) <: Tuple
        JuMP.set_optimizer(judge.master_problem, judge.master_solver[1])
    end
end
