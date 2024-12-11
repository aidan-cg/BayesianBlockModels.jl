"""
Sample group membership probabilities from the full conditional distribution of a stochastic block model.

# Arguments
- `groups::Vector{Int}`: A K-length vector where entry i contains count of nodes in group i .
- `α::Float64`: parameter of Dirichlet prior distribution on group membership probabilities 

# Returns 
- a K-length vector of reals sample from the full posterior of theta 
"""
function sample_theta(group_counts::Vector{Int}, α::Float64)
    # Parameters of full conditional follow from Dirichlet-multinomial conjugacy 
    new_param = group_counts .+ α

    # Define full conditional distribution 
    d = Dirichlet(new_param)
    
    # Sample from the full conditional distribution 
    new_theta = rand(d)

    return(new_theta)
end


"""
Sample link probabilities from the full conditional distribution of a stochastic block model.

# Arguments
- `count_mat::AbstractMatrix{<:Real}`: A K x K matrix where entry (i, j) contains the count of observed links between groups i and j.
- `group_counts::Vector{<:Real}`: A vector of length K where entry i contains the count of nodes in group i.
- `prior_params::Vector{<:Real}`: A vector [α, β], the Beta prior parameters for the link probabilities.

# Returns
- A K x K matrix of posterior samples for the link probabilities π_gh.
"""
function sample_link_probs(count_mat::AbstractMatrix{<:Real}, group_counts::Vector{<:Real}, prior_params::Vector{<:Real})
    K = size(count_mat, 1)
    α, β = prior_params

    # Initialize matrix to hold sampled link probabilities
    link_probs = zeros(Float64, K, K)

    # Sample from Beta posteriors for each group pair (g, h)
    for g in 1:K
        for h in g:K  # Only iterate over upper triangle due to symmetry
            links = count_mat[g, h]           
            if g == h
                # calculation for non-links within the same group
                total_possible_links = group_counts[g] * (group_counts[g] - 1) ÷ 2
                non_links = total_possible_links - links
            else
                # calculation for non-links between two different groups
                total_possible_links = group_counts[g] * group_counts[h]
                non_links = total_possible_links - links
            end
            
            # Posterior parameters
            α_post = α + links
            β_post = β + non_links
            # Sample from the posterior
            link_probs[g, h] = rand(Beta(α_post, β_post))
        end
    end

    return Symmetric(link_probs, :U)
end


"""
    sample_groups(theta::Vector{<:Real}, pi_mat::AbstractMatrix{<:Real}, cur_groups::Vector{<:Real}, 
    A::AbstractMatrix{<:Real}, N::Int, K::Int)

Sample new group assignments from the full conditional distribution of a stochastic block model.

# Arguments
- `theta::Vector{<:Real}`: A vector of length `K` containing the current group weights (probabilities of belonging to each group).
- `pi_mat::AbstractMatrix{<:Real}`: A `K x K` matrix of link probabilities between groups.
- `cur_groups::Vector{<:Real}`: A vector of length `N` where entry `i` indicates the current group assignment of node `i`.
- `A::AbstractMatrix{<:Real}`: An `N x N` adjacency matrix representing the observed network, where `A[i, j] = 1` if there is a link between nodes `i` and `j` and `0` otherwise.
- `N::Int`: The total number of nodes in the network.
- `K::Int`: The total number of groups.

# Returns
- A vector of length `N` containing the updated group assignments for all nodes.
"""
function sample_groups(theta::Vector{<:Real}, pi_mat::AbstractMatrix{<:Real}, cur_groups::Vector{<:Real}, 
    A::AbstractMatrix{<:Real}, N::Int, K::Int)
    for i = 1:N
        # Vector to store unnormalized log posterior probabilities
        log_c_i = zeros(K)

        for h = 1:K
            # Link probabilities for node i to member of group h 
            link_probs = pi_mat[h, cur_groups]

            # Compute log-probability of observed links
            log_probs = sum(log.(A[i, :] .* link_probs .+ (1 .- A[i,:]) .* (1 .- link_probs)))

            # Add prior contribution from theta
            log_c_i[h] = log(theta[h]) + log_probs
        end

        # Normalize to create valid probabilities (using log-sum-exp for stability)
        log_c_i_norm = log_c_i .- logsumexp(log_c_i)
        c_i_norm = exp.(log_c_i_norm)

        # Sample new group assignment for node i
        w = Weights(c_i_norm)
        cur_groups[i] = sample(1:K, w)
    end

    return cur_groups
end
