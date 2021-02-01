%Solving 3D Neo-Hookean Hypereleasticity Finite Strain problem using FEM

clear
clc

format short

%Phyics parameter
phys=struct();
phys.nu = 0.3;
phys.E = 1;

%if nonlinear problem store gradu or tensor C 
appCtx = struct();
appCtx.store = 1;

%NOTE for testing: 
% To test userf implementation, choose:
%  solver.KSP_type = 'test_res';
%  This will ignore userdf implementation and uses MATLAB's fsolve.
%
% % To test userdf implementation, choose:
%  solver.KSP_type = 'test_jac';
%  This assumes userf implementation is correct and assembles a global
%  Jacobian for each step (solver.numSteps) and compares it the 
%  global Jacobian computed from fslove and residual evaluation.
%
solver=struct();
solver.KSP_type = 'gmres'; 
solver.KSP_max_iter = 225;
solver.nonlinear_max_iter=10;
solver.global_res_tol = 1.0e-6;
solver.precond = 'OFF';
solver.numSteps = 1;


%degree of accuracy to solve with
degree = 2;  % 1 for Hex8, 2 for Hex27 
P = degree +1;
[~ , msh] = get_mesh('beam27_1e_l999_r998_6ss.exo','lex');

msh = apply_boundary_on(msh, {'ns998', 'ns999'}, {@bc_zero, @bc_bend});

dof = 3;

%physics implementation for residual evaluation
userf= @HyperFSF_ref_dav;
%physic implementation for Jacobian and action of Jacobian evaluation
userdf=@HyperFSF_dF_ref_dav;
%Forcing function
usrf_force=@HyperFS_force; 

[fem_sol, JACOB__] =  get_fem_sol(msh, dof, P, userf, userdf, usrf_force, solver, phys, appCtx);

