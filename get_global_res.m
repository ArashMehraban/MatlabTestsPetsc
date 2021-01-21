function [stored , global_res] = get_global_res(u, global_idx_map, msh, dir_bndry_val,P,userf,usrf_force,phys, store)
% GET_GLOBAL_RES evaluates the global residual and the consistent tangent
%  input:              u: vector of unknowns 
%       : global_idx_map: global map of local u's
%       :            msh: mesh object (see get_mesh function)
%       :  dir_bndry_val: Dirchlet boundary values if any
%       :          userf: user supplied code for physics
%
% output: global_res: global residual
%       : 
     
     %Allocate space for global residual for unkowns     
     global_res =zeros(size(u,1),1);
     
     %get all dirichlet boundary node_sets
     dir_bndry_nodes = get_all_dir_ns(msh);
     
     %get closure of u : a vector consisting the unknown and dirichlet boundary values
     u_closure =  get_closure_u(u,dir_bndry_nodes,dir_bndry_val,global_idx_map);     
       
     %get Weights, Basis (B) functions and their Derivatives (D0, D1 and D2)
     %D = partial_B/partial_x_i
     [B, D, W] = get_shape(P, msh.num_dims);
     
     %allocate space for storage of gradu or tensor C if needed
     if store == 1
        stored = zeros(msh.num_elem*size(D,1) , size(global_idx_map,2));
     else
         stored = 0;
     end
     
     for i=1:msh.num_elem
                   
         %get corresponding vertex coordinates for each element 
         element_vtx_coords = msh.vtx_coords(msh.conn(i,:),:);
         
         %get corresponding unknown/solution u for each element
         elem_u = u_closure(msh.conn(i,:),:);   
         
         %Compute ForcingVector
         xe= B*element_vtx_coords;
         f0 = usrf_force(xe,phys); %forcing function
         
         
         %Compute Constitutive model
         du= D*elem_u; % elemt derivative wrt xi eta zeta
         dx= D*element_vtx_coords; %element
         [dets, dXdx] = invJacobianTensor(dx); %dXdx: element inverse Jacobian
         wdetj = W.*dets;
         %phys passes the user defined Physics (such as nu, E, etc.)
         
         for j=1:size(f0,2)
             f0(:,j)=wdetj.*f0(:,j);
         end 
         
         [usrfStored, f1] = userf(du, dXdx, wdetj, phys); 
         
         % store compupted gradu from userf to be used in consistent
         % tangent operator
         if store == 1
             stored((i-1)*size(du,1)+1:i*size(du,1),:) = usrfStored;
         end
         
         % element residual evaluation
         res_e = B'*f0 + D'*f1;

         % global residual assembly 
         temp=msh.conn(i,:)';
         k=1:size(msh.conn(i,:),2);
         kk=temp(k);
         in_glb = global_idx_map(kk,:);
         kk =in_glb(in_glb>=0);
         %global residual
         global_res(kk) = global_res(kk)+ res_e(in_glb>=0);    
     end
end