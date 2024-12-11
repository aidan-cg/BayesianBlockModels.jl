module BayesianBlockModels

using StatsBase, Random, Distributions, LinearAlgebra
# Write your package code here.
export gibbs_sample
export sample_groups
export sample_link_probs
export sample_theta
export simulate_sbm
export get_link_counts

include("gibbs.jl")
include("util.jl")
include("sample.jl")
end
