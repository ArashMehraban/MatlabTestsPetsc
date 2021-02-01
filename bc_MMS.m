function boundary_dof = bc_MMS(boundary_dof, boundary_vtx, load_increment, AppCtx)
  x=boundary_vtx(:,1);
  y=boundary_vtx(:,2);
  z=boundary_vtx(:,3);
  boundary_dof(:,1) = exp(2*x).*sin(3*y).*cos(4*z);
  boundary_dof(:,2) = exp(3*y).*sin(4*z).*cos(2*x);
  boundary_dof(:,3) = exp(4*z).*sin(2*x).*cos(3*y);

end
