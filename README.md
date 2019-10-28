# Extended-target-PMBM-filter
An extended target Poisson multi-Bernoulli mixture filter using Gamma Gaussian inverse Wishart distribution
This repository contains the Matlab implementations of the Extended target Poisson multi-Bernoulli mixture filter proposed in 

Granstrom, Karl, Maryam Fatemi, and Lennart Svensson. "Poisson multi-Bernoulli mixture conjugate prior for multiple extended target filtering." IEEE Transactions on Aerospace and Electronic Systems (2019).

Full text is available at https://arxiv.org/abs/1605.06311

The filters are evaluated using the generalised optimal subpattern-assignment (GOSPA) integrated with the Gaussian Wasserstein distance

A. S. Rahmathullah, A. F. García-Fernández, and L. Svensson, Generalized optimal sub-pattern assignment metric, in 20th International
Conference on Information Fusion, 2017.

Video on GOSPA: https://www.youtube.com/watch?v=M79GTTytvCM

- main.m runs the extended target PMB filter

- data association algorithm: DB-scan + Murty's algorithm
