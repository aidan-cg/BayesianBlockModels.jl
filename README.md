# BayesianBlockModels.jl

[![Build Status](https://github.com/aidan-cg/BayesianBlockModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/aidan-cg/BayesianBlockModels.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Posterior Sampling for a Standard Stochastic Block Model
This package implements a Gibbs sampler for a standard undirected stochastic block model (SBM). The package works best with small to medium sized adjency matrices (less than 100,000 nodes). The main function `gibbs_sample` allows users to quickly sample from the posterior distribution of an SBM for basic model-based analysis of network data. Gibbs sampling takes advantage of the structure of SBMs, typically increasing the speed of sampling when compared to other MCMC algorithms like Metropolis-Hastings or Hamiltonian Monte Carlo. 

## The Model
The model assumes an undirected network and places the standard Dirichlet prior on group membership parameters and Beta priors on link probability parameters. See [Lee and Wilkinson (2019)](https://appliednetsci.springeropen.com/articles/10.1007/s41109-019-0232-2) for a review of SBMs, including the basics of Gibbs sampling for SBMs. Using this package does not require the user to understand the details of this specific inference method; users only need to understand the parameters of the model and the basics of Bayesian inferenence (e.g. how to set the parameters of a prior distribution, how to intepret samples from the posterior). 

## Basic Usage
Users input the adjency matrix and the hyperparameters into `gibbs_sample` which then outputs samples from the posterior distribution of the model. The following code simulates an adjency matrix from a small SBM (the package includes the `simulate_sbm` function) and then samples from the posterior using `gibbs_sample`

```julia
### Simulate a stochastic block model
# Set number of groups
K = 2

# Set link probabilities 
pi = [0.9 0.1;
      0.1 0.9]

# Set number of nodes
N = 100

# Simulate data from the model
A, true_groups = simulate_sbm(N, K, pi)

### Run the Gibbs sampler
# Number of samples and burn in
n_samp = 1000
n_burnin = 500

# hyperparameter for prior on group probabilities
α = 1.0

# Hyperparameters of Beta prior on link probabilities
a = 2.0
b = 2.0

# Sample from the posterior distribution of the SBM 
results = gibbs_sample(A, K, n_samp, n_burnin, α, a, b)

# Access posterior samples
group_samples = results[:group_samples]
theta_samples = results[:theta_samples]
pi_samples = results[:pi_samples]

```

## Setting the Hyperparameters
The above example uses reasonable rule-of-thumb hyperparameters for the Dirichlet and Beta priors. Users without much experience with SBMs and/or Bayesian inference can use them without much concern as long as the number of observed nodes is reasonably large. 

The parameter K determines the number of groups (sometimes called communities) in the model. Setting this parameter poorly (e.g. setting it to 3 when there are actually 5 groups) will lead to model misspecification, meaning the posterior samples cannot be used for inference. Thus, users must be careful when setting this parameter. I recommend setting K to be larger than you would expect the true value to be (if you expect 3-5 groups, you could set K to 8), and the model will hopefully set the probabilities of the excess groups close to 0 (i.e. only 5 groups will have significant membership probability at the mean of the posterior distribution). Such an approach works best with stronger shrinkage priors than are currently implemented; future work will allow for priors better suited to this approach, but users should be able to find success with the current approach with some experimentation. 

A general note: SBMs are incredibly simple models and are unlikely to accurately represent complex real-world data generating processes. My (opinionated) reccomendation for less familiar users is to use SBMs (and therefore this package) as a data exploration/first analysis tool but to be cautious of drawing strong conclusions from such an SBM. 