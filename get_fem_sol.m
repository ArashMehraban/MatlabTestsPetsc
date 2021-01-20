function u =  get_fem_sol(msh, dof, dir_bndry_nodes, dir_bndry_val,P,userf,userdf,solver, phys, store)
%GET_FEM_SOL returns the solution to the PDE
%
%input :                msh: mesh object
%      :         sz_u_field: size of unknown field (eg. : Poisson 1, Plane Strain 2)
%      :    dir_bndry_nodes: Dirichlet Boundary nodes
%      :      dir_bndry_val: Dirichlet Boundary values
%      :num_quadr_pts_in_1d: number of quadrature points in 1 dimension
%      :              userf: any knowns physics resources as a function input for global residual 
%      :             userdf: any knowns physics resources as a function inout for consistent tangent 
%      :     solver.KSP_max_iter: user specified maximum number of itertions before GMRES stops     
%      :        max_iter_nw: user specified maximum number newton steps before Newton stops
%      :     global_res_tol: user specified tolerance for norm of global_res before solver stops  
%
%output:              u: FEM Solution u
%           gl_res_norm: An array where each element contains the global residual norm after an iteration of GMRES
%               JACOB__: Jacobian of the system returned by the solver
%
% NOTE: return values that end with __ (e.g. JACOB__) are not intended to
% be used in the program per se. Be very careful before using them. 


% % % % % % % % solver = {'gmres', solver.KSP_max_iter,max_iter_nw,tol,global_res_tol};
    
  %get number of nodes in mesh
  num_nodes = msh.num_nodes;  
    
    
  %get the unknown u guess vector and global to local mapping for each u
  [u, global_idx_map] = get_global_map(num_nodes,dir_bndry_nodes,dof);
   
  if(strcmp(solver.KSP_type,'gmres'))      
      if solver.precond == 1   
         Jac = assemble_global_jac(size(u,1),global_idx_map, msh,P,dof);
      end
      
        
     t = linspace(0,1,solver.numSteps+1);
     dt = t(:,2:end);
     
     step_dir_bndry_val = dir_bndry_val;

     for m=1:solver.numSteps
    
      if solver.numSteps>1
        for ii=1:size(dir_bndry_val,1)
           step_dir_bndry_val{ii} = dir_bndry_val{ii}*dt(m);
        end
      end
    
        iter=1;
        gl_res_norm = zeros(1);
        gl_res_norm_iter=1;
        
        if(solver.KSP_max_iter > max(size(u)))
            solver.KSP_max_iter = max(size(u));
            sz_u = max(size(u));
            message = ['Number of solver.KSP_max_iter  was set to the number of unknowns: ', num2str(sz_u)];
            disp(message);            
        end
        
        while(true)
            [stored, global_res] = get_global_res(u, global_idx_map, msh, dir_bndry_val,P,userf,phys, store); 
            
            global_res_norm = norm(global_res,inf);
        
            %store the norm of global_res after each iteration
            gl_res_norm(gl_res_norm_iter) = global_res_norm;
            if(gl_res_norm_iter>1)
               if(gl_res_norm(gl_res_norm_iter) > 20*gl_res_norm(gl_res_norm_iter-1)) 
                  disp('Rtol is diverging');
                  break;
               end
            end
            gl_res_norm_iter = gl_res_norm_iter+1;
                     
            Rtol = abs(global_res_norm/ gl_res_norm(1));
            ralative_tolerance = ['Rtol: ', num2str(Rtol)];
            disp(ralative_tolerance);
            fprintf('\n');
        
            if(global_res_norm < global_res_tol || abs(global_res_norm/ gl_res_norm(1)) < solver.global_res_tol || iter > nonlinear_max_iter )
               break;
            end

            jv_fun = @(dlta_u)get_global_Jv(dlta_u, global_idx_map, msh, P, userdf,AppCtx);
            
            if strcmp(solver.precond ,'OFF')
                tic
                dlta_u= gmres(jv_fun, global_res,[], [],solver.KSP_max_iter,[]);
                toc
                u=u-dlta_u;                
            else
                disp('Jacobian Start...')
                tic
                global_jac = get_global_jac(Jac, global_idx_map, msh,P,dof, userdf,AppCtx);
                toc
                disp('Jacobian End.')

                disp('LU Start...')
                tic
                [L,U] = lu(global_jac);
                toc
                disp('LU End.')
                precond_jac_fun =@(u)(U\(L\u));  
                tic
                dlta_u= gmres(jv_fun, global_res,[], [],solver.KSP_max_iter,precond_jac_fun);
                toc
                u=u-dlta_u;               
            end
            
            iter=iter+1;
        end
    end
  end       
end