using BayesianBlockModels
using Test

# Test simulate_sbm function
@testset "simulate_sbm tests" begin
    N = 10
    K = 2
    pi = [0.8 0.2; 0.2 0.6]
    A, group_assign = simulate_sbm(N, K, pi)
    
    @test size(A) == (N, N)  # Check adjacency matrix dimensions
    @test length(group_assign) == N  # Check group assignment length
    @test issymmetric(A)  # Check if adjacency matrix is symmetric
    @test all(group_assign .<= K) && all(group_assign .>= 1)  # Check group assignments
end

# Test sample_theta function
@testset "sample_theta tests" begin
    group_counts = [5, 5]
    α = 1.0
    theta = sample_theta(group_counts, α)
    
    @test length(theta) == length(group_counts)  # Check length of theta
    @test isapprox(sum(theta), 1.0; atol=1e-6)  # Check if theta sums to 1
end

# Test sample_link_probs function
@testset "sample_link_probs tests" begin
    count_mat = [3.0 2.0; 2.0 4.0]
    group_counts = [5.0, 4.0]
    prior_params = [1.0, 1.0]
    link_probs = sample_link_probs(count_mat, group_counts, prior_params)
    
    @test size(link_probs) == size(count_mat)  # Check matrix dimensions
    @test issymmetric(link_probs)  # Check if the matrix is symmetric
    @test all(link_probs .>= 0.0) && all(link_probs .<= 1.0)  # Check if probabilities are valid
end

# Test gibbs_sample function
@testset "gibbs_sample tests" begin
    N = 10
    K = 2
    pi = [0.8 0.2; 0.2 0.6]
    A, _ = simulate_sbm(N, K, pi)
    n_samp = 10
    n_burnin = 5
    α = 1.0
    a = 1.0
    b = 1.0
    
    results = gibbs_sample(A, K, n_samp, n_burnin, α, a, b)
    group_samples = results[:group_samples]
    theta_samples = results[:theta_samples]
    pi_samples = results[:pi_samples]
    
    @test size(group_samples) == (n_samp, N)  # Check group samples dimensions
    @test size(theta_samples) == (n_samp, K)  # Check theta samples dimensions
    @test size(pi_samples) == (K, K, n_samp)  # Check pi samples dimensions
end

# Test get_link_counts function
@testset "get_link_counts tests" begin
    A = [0 1 1; 1 0 0; 1 0 0]
    groups = [1, 2, 1]
    K = 2
    count_mat = get_link_counts(groups, A, K)
    
    @test size(count_mat) == (K, K)  # Check matrix dimensions
    @test issymmetric(count_mat)  # Check if the matrix is symmetric
end
