# helper function which fetches the current objective value coef for a given variable
function getcoef(var::JuMP.VariableRef)
    terms=objective_function(var.model, AffExpr).terms
    if var in keys(terms)
        return terms[var]
    else
        return 0.0
    end
end
