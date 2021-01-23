%Solving 3D Linear Elasticity problem using FEM
clear
clc
format short

%degree of accuracy to solve with
degree = 1;
P = degree +1;

%physic implementation for residual evaluation
userf= @LinElas;
%physic implementation for Jacobian and action of Jacobian evaluation
userdf=@LinElas_dF;
%Forcing function
usrf_force=@LinElas_force; 
%Phyics parameter
phys=struct();
phys.nu = 0.3;
phys.E = 1;

%if nonlinear problem store gradu or tensor C 
store = 0; 

%NOTE: You can use MATLAB's nonlinear solver fsolve to test the correctness
%      of the residual evaluation function (get_global_res) by setting
%      the KSPS_type = 'newton_override'. With that option, solve time will
%      increase especially for larger meshes and the action of Jacobian
%      function (get_global_Jv) will be ignored.
solver=struct();
solver.KSP_type = 'gmres'; %solver.KSP_type = 'newton_override'; for fsolve
solver.KSP_max_iter = 225;
solver.nonlinear_max_iter=10;
solver.global_res_tol = 1.0e-6;
solver.precond = 'OFF';
solver.numSteps = 1;

% Use this mesh for MMS
[~ , msh] = get_mesh('cylinder8_368e_us','exo','lex');
%[~ , msh] = get_mesh('cube8_8','exo','lex');

dof = 3;
%get all Dirichlet boundary node sets
dir_bndry_nodes = get_all_dir_ns(msh);
    
%NOTE: modify userf function according to given_u
%NOTE: Proposed solution to manufacture rhs from in userf function:
given_u{1}=@(x,y,z)exp(2*x).*sin(3*y).*cos(4*z);
given_u{2}=@(x,y,z)exp(3*y).*sin(4*z).*cos(2*x); 
given_u{3}=@(x,y,z)exp(4*z).*sin(2*x).*cos(3*y);

%get vertex coordinates from mesh                                     
vtx_coords = msh.vtx_coords;                                          
%get constructed dir_bndry_vals and exac Solutions on remaining nodes 
[dir_bndry_val, exactSol] = get_exact_sol(vtx_coords,dir_bndry_nodes, given_u);

[fem_sol, Jac] =  get_fem_sol(msh, dof, dir_bndry_nodes, dir_bndry_val,P,userf,userdf,usrf_force, solver, phys, store);

error = norm(exactSol - fem_sol)/norm(fem_sol);
L2ErrMsg = strcat('L2 Error: ', num2str(error));
disp(L2ErrMsg);





    