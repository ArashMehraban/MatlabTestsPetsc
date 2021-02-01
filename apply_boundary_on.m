function msh = apply_boundary_on(msh, boundaryName, boundaryFunction)
   
   if(strcmp(boundaryName{1},'MMS'))
       msh.apply_boundary_on = msh.ns_names;
       sz_boundary_functions = size(msh.ns_names,2);
       tmp_bd_func = cell(1,sz_boundary_functions);
       for i=1:sz_boundary_functions
           tmp_bd_func{i} = boundaryFunction{1};
       end
       msh.boundary_funtions = tmp_bd_func;
   else
       msh.apply_boundary_on = boundaryName;
       sz_boundary_functions = size(boundaryName,2);
       tmp_bd_func = cell(1,sz_boundary_functions);
       msh.apply_boundary_on = boundaryName;
       for i=1:sz_boundary_functions
           tmp_bd_func{i} = boundaryFunction{i};
       end
       msh.boundary_funtions = tmp_bd_func;
   end
end
   