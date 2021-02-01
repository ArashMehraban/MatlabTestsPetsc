function [global_res, stored] = compute_residual(u, msh, dofMap, P ,userf,usrf_force,phys, load_increment, appCtx)
% GET_GLOBAL_RES evaluates the global residual and the consistent tangent
%  input:              u: vector of unknowns 
%                  userf: user supplied code for physics
%
% output: global_res: global residual
%       : 
     
     %Allocate space for global residual     
     global_res = zeros(dofMap.u_sz,1);
 
     u_closure =  get_closure_u(u,dofMap);
    
     u_closure = insert_boundary(u_closure, msh, load_increment, appCtx);
        
     %compute Weights, Basis (B) functions and their Derivatives (D0, D1 and D2)
     %D = partial_B/partial_x_i (interlaced)
     [B, D, W] = create_FEbasis_interlaced(P, msh.num_dims);
     
     %allocate space for storage of gradu or tensor C if needed
     if appCtx.store == 1
        stored = zeros(msh.num_elem*size(D,1) , dofMap.dof);
     else
         stored = 0;
     end
     
     
     for i=1:msh.num_elem
         
         elem_conn = msh.conn(i,:);
         %get corresponding vertex coordinates for each element 
         element_vtx_coords = msh.vtx_coords(elem_conn,:);
         
         %get corresponding unknown/solution u for each element
         elem_u = u_closure(elem_conn,:);   
         
         %Compute ForcingVector
         xe= B*element_vtx_coords;
         f0 = usrf_force(xe,phys); %forcing function
         
         
         %Compute Constitutive model
         du= D*elem_u; % element derivative wrt xi eta zeta
         dx= D*element_vtx_coords; %element
         [dets, dXdx] = invJacobianTensor(dx); %dXdx: element inverse Jacobian
         wdetj = W.*dets;
         
         
         for j=1:size(f0,2)
             f0(:,j)=wdetj.*f0(:,j);
         end 
         
         %phys passes the user defined Physics (such as nu, E, etc.)
         [usrfStored, f1] = userf(du, dXdx, wdetj, phys); 
         
         % store compupted gradu from userf to be used in consistent
         % tangent operator
         if appCtx.store == 1
             stored((i-1)*size(du,1)+1:i*size(du,1),:) = usrfStored;
         end
         
         % element residual evaluation
         res_e = B'*f0 + D'*f1;
         
         

         % global residual assembly 
         temp=elem_conn';
         k=1:size(elem_conn,2);
         kk=temp(k);
         in_glb = dofMap.map(kk,:);
         kk =in_glb(in_glb ~= 0);
         %global residual
         global_res(kk) = global_res(kk)+ res_e(in_glb ~= 0);    
     end
end