function [u,JACOB__] =  get_fem_sol(msh, dof, P, userf, userdf, usrf_force, solver, phys, appCtx)
%GET_FEM_SOL returns the FEM solution to the PDE and gobal Jacobian if
%desired.
%
%input :            msh: mesh object
%      :            dof: size of unknown field (eg. : Poisson 1, Plane Strain 2)
%      :dir_bndry_nodes: Dirichlet Boundary nodes
%      :  dir_bndry_val: Dirichlet Boundary values
%      :              P: number of quadrature points in 1 dimension
%      :          userf: user-implemented physics for residual evaluation 
%      :         userdf: user-implemented physics for global Jacobian evaluation (the derivative of userf)
%
%output:              u: FEM Solution u
%               JACOB__: Jacobian of the system returned by the solver
%
    
  %dofMap is a struct
  dofMap = create_dof_map(msh,dof);
  
  %u: initial unknown vector
  u = zeros(dofMap.u_sz,1);
  
  %For testing userdf in compute_Jacobian_action (action of Jacobian)
  %Allocating Space
  if(strcmp(solver.KSP_type,'test_jac'))
      stepJacs = cell(solver.numSteps,1);
      for m=1:solver.numSteps
          stepJacs{m}= zeros(size(u,1),size(u,1));
      end
      Iu = eye(size(u,1),size(u,1));
      
      JACOB__ = struct();
      JACOB__.stepJacs = stepJacs;
      JACOB__.stepFsolveJacs = stepJacs;
      JACOB__.diff_norms = zeros(solver.numSteps,1);
      JACOB__.diff_norms_computation = 'norm(JACOB__.stepFsolveJacs{step_i} - JACOB__.stepJacs{step_i})/norm(JACOB__.stepJacs{step_i})';
  end
  
  %compute load increments based on number of steps (solver.numSteps)
  pseudo_time = linspace(0,1, solver.numSteps+1);
  load_increments = pseudo_time(:,2:end);
  
  
  if(isfield(appCtx,'vtk_filename'))
      vtk_files = cell(solver.numSteps+1,1);
      for s=1:solver.numSteps+1
          vtk_files{s} = strcat(appCtx.vtk_filename , num2str(s), '.vtk');    
      end
       %Store undeformed geometry in vtk_files{1}
        u_closure = get_closure_u(u,dofMap); 
        %graph at no boundary applied
        processVtk(vtk_files{1},appCtx.origConn, msh,u_closure);
  end
  
 
  
  for m=1:solver.numSteps
    
    if(strcmp(solver.KSP_type,'gmres'))      
      if solver.precond == 1   
         Jac = assemble_global_jac(size(u,1),global_idx_map, msh,P,dof);
      end
    
        iter=1;
        gl_res_norm = zeros(1);
        gl_res_norm_iter=1;
        
        if(solver.KSP_max_iter > max(size(u)))
            solver.KSP_max_iter = max(size(u));
            sz_u = max(size(u));
            message = ['Number of solver.KSP_max_iter was set to the number of unknowns: ', num2str(sz_u)];
            disp(message);            
        end
        
        while(true)
            [global_res, stored] = compute_residual(u, msh, dofMap, P ,userf,usrf_force,phys, load_increments(m), appCtx); 
            JACOB__ = 'Matrix-Free Applied';
            
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
        
            if(global_res_norm < solver.global_res_tol || abs(global_res_norm/ gl_res_norm(1)) < solver.global_res_tol || iter > solver.nonlinear_max_iter )
               break;
            end

            Jacobian_action_fun = @(dlta_u)compute_Jacobian_action(dlta_u, msh, dofMap, P ,userdf, phys, stored);
            
            if strcmp(solver.precond ,'OFF')
                tic
                dlta_u= gmres(Jacobian_action_fun, global_res,[], [],solver.KSP_max_iter,[]);
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
                dlta_u= gmres(Jacobian_action_fun, global_res,[], [],solver.KSP_max_iter,precond_jac_fun);
                toc
                u=u-dlta_u;               
            end
            
            iter=iter+1;
        end
        
        if(isfield(appCtx,'vtk_filename'))
            u_closure = get_closure_u(u,dofMap); 
            u_closure = insert_boundary(u_closure, msh, load_increments(m), appCtx);
            processVtk(vtk_files{m+1},appCtx.origConn, msh,u_closure);
        end
        
    end %end for if gmres solver
   
   %For testing global residual evaluation. This will ignore compute_Jacobian_action
   if(strcmp(solver.KSP_type,'test_res'))
      fun = @(u)compute_residual(u, msh, dofMap, P ,userf,usrf_force,phys, load_increments(m), appCtx);
      [u,~,~,~,JACOB__] = fsolve(fun, u);
   end %end of if(strcmp(solver.KSP_type,'test_res'))
   
   %For testing action of global Jacobian. This compares the global Jacobian
   %produced from compute_Jacobian_action with Jacobian (JACOB__) from fsolve 
   if(strcmp(solver.KSP_type,'test_jac'))
       %compute golbal residual for each step
       [~, stored] = compute_residual(u, msh, dofMap, P ,userf,usrf_force,phys, load_increments(m), appCtx);
       %populate global Jacobian based on the global residual for each step
       for i=1:size(u,1)                          
           JACOB__.stepJacs{m}(i,:)=compute_Jacobian_action(Iu(:,i), msh, dofMap, P ,userdf, phys, stored);
       end
       %Jacobian from fsolve for each step
       fun = @(u)compute_residual(u, msh, dofMap, P ,userf,usrf_force,phys, load_increments(m), appCtx);
       [u,~,~,~,fsolvJ] = fsolve(fun, u);% 
       JACOB__.stepFsolveJacs{m} = fsolvJ; 
       
       %norm of difference between stepFsolveJacs and stepJacs
       JACOB__.diff_norms(m) = norm(JACOB__.stepFsolveJacs{m} -JACOB__.stepJacs{m})/norm(JACOB__.stepJacs{m});
   end %if(strcmp(solver.KSP_type,'test_jac'))   
  end %end for numSteps for-loop
end