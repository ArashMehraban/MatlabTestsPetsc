%Solving 3D Linear Elasticity problem using FEM
clear
clc
format short

%Phyics parameter
phys=struct();
phys.nu = 0.3;
phys.E = 1;

%if nonlinear problem store gradu or tensor C 
appCtx = struct();
appCtx.store = 0;

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
degree = 1;  % 1 for Hex8, 2 for Hex27 
P = degree +1;
% Use this mesh for MMS
[~ , msh] = get_mesh('cylinder8_110e_us.exo','lex');

dof = 3;

%physics implementation for residual evaluation
userf= @LinElasF;
%physic implementation for Jacobian and action of Jacobian evaluation
userdf=@LinElasF_dF;
%Forcing function
usrf_force=@LinElas_force_MMS; 

msh = apply_boundary_on(msh, {'MMS'}, {@bc_MMS});
%msh = apply_boundary_on(msh, {'ns1', 'ns2'}, {@bc_zero, @bc_bend});
    
%compute exact solution based on given_u 
exactSol = compute_exact_solution_3D_elas(msh);

[fem_sol, JACOB__] =  get_fem_sol(msh, dof, P, userf, userdf, usrf_force, solver, phys, appCtx);

error = norm(exactSol - fem_sol)/norm(fem_sol);
L2ErrMsg = strcat('L2 Error: ', num2str(error));
disp(L2ErrMsg);





    