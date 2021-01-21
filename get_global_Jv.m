function Jv = get_global_Jv(dlta_u, global_idx_map, msh,P, userdf,phys, stored)
% GET_GLOBAL_Jv evaluates the action of Jaccobian (consistent
% tangent) on dlta_u. It produces a vector mfj=J*dlta_u
%  input:              dlta_u: vector of variation of unknowns 
%       :      global_idx_map: global map of local u's
%       :                 msh: mesh object 
%       : num_quadr_pts_in_1d: number of quadrature points in 1D
%       :              userdf: user supplied code for get_userdf function
%
% output: Jv: Jacobian-vector Product: the action of consistant tangent on dlta_u
%               
     
     %Allocate space for global Matric Free Jacobian     
     Jv =zeros(size(dlta_u,1),1);

     %get all dirichlet boundary node_sets
     dir_bndry_nodes = get_all_dir_ns(msh);
      
     delta_u_closure = get_closure_dlta_u(dlta_u,dir_bndry_nodes,global_idx_map); 
     
     %get Weights, Basis (B) functions and their Derivatives (D0, D1 and D2)
     %D = partial_B/partial_x_i
     [B, D, W] = get_shape(P, msh.num_dims);
  
               
     for i=1:msh.num_elem
         
        %get corresponding vertex coordinates for each element 
         element_vtx_coords = msh.vtx_coords(msh.conn(i,:),:);
         
         %get corresponding unknown/solution u and du for each element
         elem_dlta_u = delta_u_closure(msh.conn(i,:),:); 
                  
         grad_dlta_ue = D*elem_dlta_u; % elemt derivative
         dxdX= D*element_vtx_coords; %element
         [dets, dXdx] = invJacobianTensor(dxdX); %dXdx: element inverse Jacobian
         wdetj = W.*dets;
                   
         %stored passes the userf defined Physics (such as nu, E, etc.)
         f = userdf(grad_dlta_ue, stored((i-1)*size(grad_dlta_ue,1)+1:i*size(grad_dlta_ue,1),:) ,dXdx, wdetj, phys);
         
         %Matrix Free Jacobian per element (This is a Vector)
         Jv_e = [B' D']*f;
  
         % global (the action of) jacobian assembly 
         temp=msh.num_elem(i,:)';
         k=1:size(msh.num_elem(i,:),2);
         kk=temp(k);
         in_glb = global_idx_map(kk,:);
         kk =in_glb(in_glb>=0);
         %global action of Jacobian
         Jv(kk) = Jv(kk)+ Jv_e(in_glb>=0); 
     end
end