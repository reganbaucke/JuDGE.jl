# The structure of the main judge model
mutable struct JuDGEModel
    tree::Tree
    master::JuMP.Model
    deteq::JuMP.Model
    subprob::Dict{Node,JuMP.Model}
    mastervar::Dict{Node,Dict{Symbol,Any}}
    mastercon::Dict{Node,Dict{Symbol,Any}}
    covercon::Dict{Node,Any}
    isbuilt::Bool
    isfixed::Bool
    isbuiltdeteq::Bool
    lb::Float64
    buildexpansionvariables
    buildsubs
    expansioncosts
    function JuDGEModel(tree::Tree)
        this = new()
        this.tree = tree
        this.isbuilt=false
        this.isfixed=false
        this.isbuiltdeteq=false
        this.lb=-Inf
        return this
    end
end
