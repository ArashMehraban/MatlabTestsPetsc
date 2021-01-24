# 3D Elasticty, Hyperelasticity with MATLAB

You can run compressible 3D Linear Elasticity and Hyperelasticity at finite strain using refernce configuration and current configuration.

Current configuration code is work in progress! Its residual evaluation is correct, but the action of Jacobian needs work.

Try:

`runLinElas.m` to run 3D Linear Elasticity problem

`runHyperFS\_ref.m` to run Neo-Hookean Hyperelasticity at finite strain using reference config

`runHyperFS\_cur.m` to run Neo-Hookean Hyperelasticity at finite strain using current config

This is done using Dirichlet Boundary conditions only. No Neumann yet!

