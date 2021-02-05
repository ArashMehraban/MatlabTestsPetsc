# 3D Elasticty, Hyperelasticity with MATLAB

You can run compressible 3D Linear Elasticity and Hyperelasticity at finite strain using refernce configuration and current configuration.

Current configuration code is work in progress! Its residual evaluation is correct, but the action of Jacobian needs work.

Try:

`runLinElas.m` to run 3D Linear Elasticity problem

`runHyperFS_ref.m` to run Neo-Hookean Hyperelasticity at finite strain using reference config

`runHyperFS_ref_dav.m` to run Neo-Hookean Hyperelasticity at finite strain using Davydov Algorithm 2 

This is done using Dirichlet Boundary conditions only. No Neumann yet!

![](HyperFS_ref_dav.gif)
