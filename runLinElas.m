%Solving 3D Linear Elasticity problem using FEM
clear
clc
format short

%physic implementation for residual evaluation
userf= @LinElasF;
%physic implementation for Jacobian and action of Jacobian evaluation
userdf=@LinElasF_dF;
%Forcing function
usrf_force=@LinElas_force; 
%Phyics parameter
phys=struct();
phys.nu = 0.3;
phys.E = 1;

%if nonlinear problem store gradu or tensor C 
store = 0; 

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
[~ , msh] = get_mesh('cylinder8_110e_us','exo','lex');

dof = 3;
%get all Dirichlet boundary node sets
dir_bndry_nodes = get_all_dir_ns(msh);
    
%NOTE: modify userf function according to given_u
%NOTE: Proposed solution to manufacture rhs from in userf function:
given_u{1}=@(x,y,z)exp(2*x).*sin(3*y).*cos(4*z);
given_u{2}=@(x,y,z)exp(3*y).*sin(4*z).*cos(2*x); 
given_u{3}=@(x,y,z)exp(4*z).*sin(2*x).*cos(3*y);
                                      
%get constructed dir_bndry_vals and exac Solutions on remaining nodes 
[dir_bndry_val, exactSol] = get_exact_sol(msh,dir_bndry_nodes, given_u);

[fem_sol, JACOB__] =  get_fem_sol(msh, dof, dir_bndry_nodes, dir_bndry_val,P,userf,userdf,usrf_force, solver, phys, store);

error = norm(exactSol - fem_sol)/norm(fem_sol);
L2ErrMsg = strcat('L2 Error: ', num2str(error));
disp(L2ErrMsg);





    