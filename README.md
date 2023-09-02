# Semi-Supervised-Learning-in-Multilayer-Hypergraphs-by-Coordinate-Descent

This repository contains the codes of the paper "Semi Supervised Learning in Multilayer Hypergraphs by Coordinate Descent" by Sara Venturini, Andrea Cristofari, Francesco Rinaldi, Francesco Tudisco.

## Compared Methods

- CCD: Cyclic Coordinate Descent Method
- GCD: Greedy Coordinate Descent Method
- RCD: Random Coordinate Descent Method
- GD: Gradient Descent Method

## Synthetic networks (SBM Stochastic Block Model) 
adjacent\_matrix\_generator: creates a single-layer graph using the Stochastic Block Model 

## Real Datasets
### Multilayer datasets
From https://github.com/melopeo/PM_SSL/tree/master/realworld_datasets \
P. Mercado, F. Tudisco, and M. Hein, Generalized Matrix Means for Semi-Supervised Learning with Multilayer Graphs. In NeurIPS 2019.
### Single layer hypergraphs
In folder Data. \
Incident and adjacency matrices created using code AB_hypergraphs from data in https://www.cs.cornell.edu/~arb/data/ \
P. S. Chodrow, N. Veldt, and A. R. Benson. Generative hypergraph clustering: From blockmodels to modularity. Science Advances, 7(28):eabh1303, 2021.\


## Evaluation of the final partitions
- confusion_matrix: calculates the confusion matrix. It is used in "wrong" function and input of "NMI" function.
- reindex_com: reindexes communities. It is used in "confusion_matrix" function.
- wrong: counts the number of nodes in the wrong community. It is used to calculate the Accuracy of a partition.

## Tests on Synthetic Networks (SBM) 
experiments_ART_2: run the tests on synthetic networks SBM reported in the paper.
- efficiency_plot_ART_2: tests each methods on synthetic network SBM 
- plot_eff_ART: does efficency plots objective function and accuracy in terms of number of flops

## Tests on Real Datasets
### Quadratic regularization (p=2)
experiments_REAL_2: run the tests on real-world networks (multilayer networks and hypergraps) with p=2 reported in the paper.
- efficiency_plot_REAL_2: tests each methods on real-world network 
- plot_eff_REAL: does efficency plots objective function and accuracy in terms of number of flops
### p-Laplacian regularization (p!=2)
experiments_REAL_p: run the tests on real-world networks (multilayer networks and hypergraps) with p!=2 reported in the paper.
- efficiency_plot_REAL_p: tests each methods on real-world network 
- plot_eff_REAL: does efficency plots objective function and accuracy in terms of number of flops

## Functions' definitions
- fun_obj: objective function initialization - just p-Laplacian term
- phi_p: function used to define the p-Laplacian 

## Reference paper
"Semi Supervised Learning in Multilayer Hypergraphs by Coordinate Descent" by Sara Venturini, Andrea Cristofari, Francesco Rinaldi, Francesco Tudisco.

## Authors
- Sara Venturini (e-mail: sara.venturini@math.unipd.it)
- Andrea Cristofari (e-mail: andrea.cristofari@uniroma2.it)
- Francesco Rinaldi (e-mail: rinaldi@math.unipd.it)
- Francesco Tudisco (e-mail: francesco.tudisco@gssi.it)

