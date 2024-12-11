function get_link_counts(group_counts, A, K)
    count_mat = zeros(K, K)

    # iterate over groups 
    for i = 1:K
        # first check if any nodes belong to group 
        # if so, leave the row empty 
        if i ∉ group_counts
            continue 
        end
        
        groups_i_mask = group_counts .== i 
        # the matrix is symmetric, so only need to calculate upper triangle 
        for j = i:K
            groups_j_mask = group_counts .== j 
            if i == j 
                # within group transition matrix is symmetric, so take lower triangle of it 
                count_mat[i,j] = sum(LowerTriangular(A[groups_i_mask , groups_j_mask]))
            else
                # across group trans mat is not symmetric 
                count_mat[i,j] = sum(A[groups_i_mask , groups_j_mask])
            end
        end

    end

    return(Symmetric(count_mat, :U))
end

function logsumexp(vec)
    max_val = maximum(vec)
    return max_val + log(sum(exp.(vec .- max_val)))
end

function simulate_sbm(N, K, pi)
    if size(pi, 1) != K 
        println("Dimensions of pi and K do not match")
        return -1 
    elseif size(pi, 1) != size(pi, 2)
        println("pi is not a square matrix")
        return -1
    end

    # we assign nodes to groups with uniform probability 
    group_prob = 1/K * ones(K)
    group_values = 1:K 

    w = Weights(group_prob)

    group_assign = sample(group_values, w, N)
    
    # matrix of link probabilities 
    pi_mat = zeros(N, N)
    for i in 1:N
        pi_mat[i,:] = pi[group_assign[i], group_assign]
        pi_mat[i,i] = 0 
    end
    
    # use link probabilities to generate links 
    A = rand(N, N) .< pi_mat 
    A = Symmetric(A .- Diagonal(diag(A)))

    return(A, group_assign)
end 
