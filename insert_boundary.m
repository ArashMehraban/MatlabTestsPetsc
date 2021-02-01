function u_closure = insert_boundary(u_closure, msh, load_increment, AppCtx)

  sz_bd = size(msh.apply_boundary_on,2);
  bndry_nodes = get_boundary_nodes(msh);
  for i=1:sz_bd
      boundary_dof = u_closure(bndry_nodes{i},:);
      boundary_vtx = msh.vtx_coords(bndry_nodes{i},:);
      u_closure(bndry_nodes{i},:) = msh.boundary_funtions{i}(boundary_dof, boundary_vtx, load_increment, AppCtx); 
  end
  
end