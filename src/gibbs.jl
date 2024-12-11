"""
Gibbs sampler for a stochastic block model with a fixed number of groups.

# Arguments
- `A::AbstractMatrix{<:Real}`: The adjacency matrix of the network (N x N), where `A[i, j]` is 1 if there is a link between nodes `i` and `j`, and 0 otherwise.
- `K::Int`: The number of groups in the stochastic block model.
- `n_samp::Int`: The number of posterior samples to draw (excluding burn-in iterations).
- `n_burnin::Int`: The number of burn-in iterations to discard before collecting samples.
- `α::Float64`: The concentration parameter for the Dirichlet prior on group weights.
- `a::Float64`: The α parameter for the Beta prior on link probabilities.
- `b::Float64`: The β parameter for the Beta prior on link probabilities.

# Returns
A dictionary containing the following keys:
- `:group_samples`: A matrix of size `n_samp x N`, where each row contains the group assignments for all nodes at a given iteration.
- `:theta_samples`: A matrix of size `n_samp x K`, where each row contains the sampled group weights for a given iteration.
- `:pi_samples`: A 3D array of size `K x K x n_samp`, where each slice `pi_samples[:, :, t]` contains the sampled link probabilities for iteration `t`.
"""
function gibbs_sample(A::AbstractMatrix{<:Real}, K::Int, n_samp::Int, n_burnin::Int, α::Float64, a::Float64, b::Float64)
    N = size(A, 1)

    # Initialize parameter values
    theta_init = 1/K .* ones(K)
    c_init = sample(1:K, Weights(theta_init), N)  # Uniform initial group assignments
    pi_mat_init = ones(K, K) .* 0.5  # Initial link probabilities (uniform prior)

    # Preallocate storage for samples
    group_samples = Matrix{Int}(undef, n_samp, N)  # To store group assignments
    theta_samples = Matrix{Float64}(undef, n_samp, K)  # To store group weights
    pi_samples = Array{Float64}(undef, K, K, n_samp)  # To store link probabilities

    # Gibbs sampling
    for iter in 1:(n_burnin + n_samp)
        # Update link counts and group counts
        link_counts_i = get_link_counts(c_init, A, K)
        group_counts_i = [sum(c_init .== k) for k in 1:K]

        # Update theta, pi_mat, and group assignments
        theta_init = sample_theta(group_counts_i, α)
        pi_mat_init = sample_link_probs(link_counts_i, group_counts_i, [a, b])
        c_init = sample_groups(theta_init, pi_mat_init, c_init, A, N, K)

        # Store samples after burn-in
        if iter > n_burnin
            save_idx = iter - n_burnin
            group_samples[save_idx, :] = c_init
            theta_samples[save_idx, :] = theta_init
            pi_samples[:, :, save_idx] = pi_mat_init
        end
    end

    # Return results in a dictionary for easy access
    return Dict(
        :group_samples => group_samples,
        :theta_samples => theta_samples,
        :pi_samples => pi_samples
    )
end