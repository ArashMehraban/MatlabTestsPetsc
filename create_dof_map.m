function dofMap = create_dof_map(msh,dof)
   
   dofMap=struct();
   dof_map = -ones(msh.num_nodes, dof)';
   dir_bndry_nodes_cell_array = get_boundary_nodes(msh);
   dir_bndry_nodes_idx = unique(cell2mat(dir_bndry_nodes_cell_array));
   dof_map(:,dir_bndry_nodes_idx) = 0;
   u_closure_sz = size(dof_map,2)*dof;
   j=1;
   for i=1:u_closure_sz
      if(dof_map(i) ~= 0)
          dof_map(i) = j;
          j = j+1;          
      end
   end
   dofMap.map = dof_map';
   
   uknown_nodes_idx = (1:msh.num_nodes)';
   uknown_nodes_idx(uknown_nodes_idx(dir_bndry_nodes_idx))=[];
   
   dofMap.u_sz = max(max(dof_map)); 
   dofMap.u_closure_sz = u_closure_sz;
   dofMap.boundary_nodes_idx = dir_bndry_nodes_idx;
   dofMap.u_nodes_idx = uknown_nodes_idx;
   dofMap.dof = dof;
end