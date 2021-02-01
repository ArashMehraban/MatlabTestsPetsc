function bndry_nodes = get_boundary_nodes(msh)


  sz_bd = size(msh.apply_boundary_on,2);
  bndry_nodes = cell(sz_bd,1);
  for i=1:sz_bd
      bndry_nodes{i} = msh.(msh.apply_boundary_on{i});      
  end
end