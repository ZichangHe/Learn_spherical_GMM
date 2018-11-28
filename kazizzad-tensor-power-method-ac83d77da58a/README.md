#What is this repository about#

Power Method is a robust tensor CP decomposition involves decomposing a tensor into summation of its rank-one component. In this repo we provide simple Matlab code for tensor power method. We show how this code whitens the input tensor by taking the SVD of its diagonal slice and apples iterative algorithm on whitened tensor to compute eigenvalues and eigenvectors. We illustrate how it computes the eigenvalues and eigenvectors of whitened tensor by iterative algorithm. At each iteration it multiplies the input tensor by a vector and normalize the result vector. Then multiply the tensor again by the result vector and normalize it again. It does it on and on until the iterative algorithm converges to a stable vector. We store this final vector and its normalization coefficient. For next step, it deflates the corresponding rank one tensor from the original whitened tensor and update the original whitened tensor. After that it pursues the iteration and deflation until the original whitened tensor vanishes. Then we show how it dewhitens the stored vectors and their corresponding coefficients to come up with eigenvalues and eigenvectors original tensor.
### Contents ###
The PowerIteration.m is main file in this repository and the the variables are described in the code.
The outputs are the eigenvalues,eigenvectors an the figures of power method convergence and the comparison between eigenvalues and eigenvectors with their estimated ones.

###Credit by:###
Furong Huang and Kamyar Azizzadenesheli
###Contact###
kazizzad@uci.edu