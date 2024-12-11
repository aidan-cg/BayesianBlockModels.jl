module BayesianBlockModels

# Write your package code here.
export gibbs_sample
export sample_groups
export sample_link_probs
export sample_theta
export simulate_sbm

include("gibbs.jl")
include("util.jl")
include("sample.jl")
end
