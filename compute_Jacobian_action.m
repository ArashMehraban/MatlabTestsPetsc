function Jv = compute_Jacobian_action(dlta_u, msh, dofMap, P ,userdf, phys, stored)
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
     %Allocate space for action of Jacobian      
     Jv = zeros(dofMap.u_sz,1);
 
     dlta_u_closure =  get_closure_u(dlta_u,dofMap);
        
     %compute Weights, Basis (B) functions and their Derivatives (D0, D1 and D2)
     %D = partial_B/partial_x_i (interlaced)
     [B, D, W] = create_FEbasis_interlaced(P, msh.num_dims);
      
               
     for i=1:msh.num_elem
         
         elem_conn = msh.conn(i,:);
         %get corresponding vertex coordinates for each element 
         element_vtx_coords = msh.vtx_coords(elem_conn,:);
         
         %get corresponding unknown/solution u for each element
         elem_dlta_u = dlta_u_closure(elem_conn,:);
         i;
                           
         ddu = D*elem_dlta_u; % elemt derivative
         dxdX= D*element_vtx_coords; %element
         [dets, dXdx] = invJacobianTensor(dxdX); %dXdx: element inverse Jacobian
         wdetj = W.*dets;
         
         dlta_ue = B*elem_dlta_u;
         
         if stored == 0
             f = userdf(dlta_ue, ddu, stored ,dXdx, wdetj, phys);
         else
             %stored passes the userf defined Physics (such as nu, E, etc.)
             f = userdf(dlta_ue, ddu, stored((i-1)*size(ddu,1)+1:i*size(ddu,1),:) ,dXdx, wdetj, phys);
         end
         
         %Matrix Free Jacobian per element (This is a Vector)
         Jv_e = [B' D']*f;
           
         % global (the action of) jacobian assembly 
         temp=elem_conn';
         k=1:size(elem_conn,2);
         kk=temp(k);
         in_glb = dofMap.map(kk,:);
         kk =in_glb(in_glb ~= 0);
         %global action of Jacobian
         Jv(kk) = Jv(kk)+ Jv_e(in_glb ~= 0); 
     end
end