function exactSol = compute_exact_solution_3D_elas(msh)
% GET_EXCAT_SOL produces the exact solution based on a given_u 

    given_u{1}=@(x,y,z)exp(2*x).*sin(3*y).*cos(4*z);
    given_u{2}=@(x,y,z)exp(3*y).*sin(4*z).*cos(2*x); 
    given_u{3}=@(x,y,z)exp(4*z).*sin(2*x).*cos(3*y);
   
    x = msh.vtx_coords(:,1);
    y = msh.vtx_coords(:,2);
    z = msh.vtx_coords(:,3);
    
    exactSol = zeros(msh.num_nodes, size(msh.vtx_coords,2));
    
    for i=1:msh.num_dims
       exactSol(:,i) = given_u{i}(x,y,z);
    end
    
    dir_bndry_nodes_cell_array = get_boundary_nodes(msh);
    dir_bndry_nodes_idx = unique(cell2mat(dir_bndry_nodes_cell_array));
    exactSol(dir_bndry_nodes_idx,:) = [];
      
    exactSol = exactSol(:);
end