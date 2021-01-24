%Solving 3D Neo-Hookean Hypereleasticity Finite Strain problem using FEM
%with 2 structured elements
clear
clc
format short

%physic implementation for residual evaluation
userf= @HyperFSF_ref;
%physic implementation for Jacobian and action of Jacobian evaluation
userdf=@HyperFSF_dF_ref;
%Forcing function
usrf_force=@HyperFS_force; 
%Phyics parameter
phys=struct();
phys.nu = 0.3;
phys.E = 1;

%if nonlinear problem store gradu or tensor C 
store = 1; 

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
degree = 1; % 1 for Hex8, 2 for Hex27
P = degree +1;
[~ , msh] = get_mesh('cube8_8','exo','lex');

dof = 3;
%get all Dirichlet boundary node sets
dir_bndry_nodes = get_all_dir_ns(msh);
    
%NOTE: modify userf function according to given_u
given_u{1}=@(x,y,z)0;
given_u{2}=@(x,y,z)0;
given_u{3}=@(x,y,z)0.5*z;
%get vertex coordinates from mesh                                     

[dir_bndry_val, ~] = get_exact_sol(msh,dir_bndry_nodes, given_u);
dir_bndry_val{1} = 0* dir_bndry_val{1};

[fem_sol, JACOB__] =  get_fem_sol(msh, dof, dir_bndry_nodes, dir_bndry_val,P,userf,userdf,usrf_force, solver, phys, store);




    