function [u,JACOB__] =  get_fem_sol(msh, dof, dir_bndry_nodes, dir_bndry_val,P,userf,userdf,usrf_force,solver, phys, store)
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
    
  %u: initial unknown vector 
  %global_idx_map: dof-map for local-to-global and global-to-local
  [u, global_idx_map] = get_global_map(msh.num_nodes,dir_bndry_nodes,dof);
  
  t = linspace(0,1,solver.numSteps+1);
  dt = t(:,2:end);
     
  step_dir_bndry_val = dir_bndry_val;
  
  %For testing userdf in get_global_Jv (action of Jacobian)
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
      JACOB__.diff_norms_computation = 'norm(JACOB__.stepFsolveJacs{step_i} -JACOB__.stepJacs{step_i})/norm(JACOB__.stepJacs{step_i})';
  end

  for m=1:solver.numSteps
    
    if solver.numSteps>1
      for ii=1:size(dir_bndry_val,1)
         step_dir_bndry_val{ii} = dir_bndry_val{ii}*dt(m);
      end
    end
   
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
            message = ['Number of solver.KSP_max_iter  was set to the number of unknowns: ', num2str(sz_u)];
            disp(message);            
        end
        
        while(true)
            [global_res, stored] = get_global_res(u, global_idx_map, msh, step_dir_bndry_val,P,userf,usrf_force,phys, store); 
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

            jv_fun = @(dlta_u)get_global_Jv(dlta_u, global_idx_map, msh,P, userdf,phys, stored);
            
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
    end %end for if gmres solver
   
   %For testing global residual evaluation. This will ignore get_global_Jv
   if(strcmp(solver.KSP_type,'test_res'))
      fun = @(u)get_global_res(u, global_idx_map, msh, step_dir_bndry_val,P,userf,usrf_force,phys, store);
      %options = optimoptions(@fsolve,'Algorithm','Levenberg-Marquardt');%, 'TolX',tol); %,'Jacobian','on');
      [u,~,~,~,JACOB__] = fsolve(fun, u);% ,options);
   end %end of if(strcmp(solver.KSP_type,'test_res'))
   
   %For testing action of global Jacobian. This assumes global residual
   %evaluation is correctly implemented and compare the global Jacobian
   %produced from get_global_Jv with JACOB__ from fsolve 
   if(strcmp(solver.KSP_type,'test_jac'))
       %compute golbal residual for each step
       [~, stored] = get_global_res(u, global_idx_map, msh, step_dir_bndry_val,P,userf,usrf_force,phys, store);
       %populate global Jacobian based on the global residual for each step
       for i=1:size(u,1)
           JACOB__.stepJacs{m}(i,:)=get_global_Jv(Iu(:,i), global_idx_map, msh,P, userdf,phys, stored);
       end
       %Jacobian from fsolve for each step
       fun = @(u)get_global_res(u, global_idx_map, msh, step_dir_bndry_val,P,userf,usrf_force,phys, store);
       [u,~,~,~,fsolvJ] = fsolve(fun, u);% 
       JACOB__.stepFsolveJacs{m} = fsolvJ; 
       
       %norm of difference between stepFsolveJacs and stepJacs
       JACOB__.diff_norms(m) = norm(JACOB__.stepFsolveJacs{m} -JACOB__.stepJacs{m})/norm(JACOB__.stepJacs{m});
   end %if(strcmp(solver.KSP_type,'test_jac'))   
  end %end for numSteps for-loop
end