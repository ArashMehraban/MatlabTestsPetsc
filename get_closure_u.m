function u_closure =  get_closure_u(u,dofMap)

   u_closure = reshape(zeros(dofMap.u_closure_sz,1),[],dofMap.dof);
   u_closure(dofMap.u_nodes_idx,:) = reshape(u',dofMap.dof,[])';
end