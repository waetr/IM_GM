# Influence Allocation

## Main Idea

Given a set of nodes and a diffusion model, we want to allocate each node's contribution to the propagation fairly.

## Motivation

Application: 
* Help advertisers pay influencers a fair price
* Rank seeds in IM / CIM
* Used as a feature for downstream prediction tasks (e.g., user retention, behavior conversion), maybe supported by Tencent.

Flaws of existing solutions:
* Most scoring methods (e.g., degree, pagerank, Chen's Shapley centrality) can't reflect fair contributions, leading to underpayments or overpayments.
* Jing Tang's "influence contribution" is based on the propagation result, which can lead to another issue: the real influential node may never have been selected in the previous propagation event, causing inaccurate prediction. 

## Our Method

In the triggering model, for each live/block edge graph instance, we allocate seeds based on the number of nodes they first influence.

The overall allocation to a seed is its expected value across all live/block edge graph instances.

We can prove that the sum of allocations for individual seeds equals the influence spread, while experiments show that other methods (degree, pagerank, centrality, single-node spread) result in significant bias.

## Remaining Issues

* How can we demonstrate the fairness of our method compared to marginal spread? (Note that the sum of marginal spread also equals the overall spread).
* When multiple seeds have the same shortest distance to a node, how should we distribute the allocation for that node to each seed?

## Future Work