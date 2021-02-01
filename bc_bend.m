function boundary_dof = bc_bend(boundary_dof, boundary_vtx, load_increment, AppCtx)
  
  x = boundary_vtx(:,1);

%   boundary_dof(:,1) = 0;
  %boundary_dof(:,3) = -0.2; % 50% the size of the beam in the x-coords bent in -Y direction (beam_8e_l999_r998_6ss.png)  
  boundary_dof(1,3) = 0;
  boundary_dof(2,3) = -0.1;
  boundary_dof(3,3) = 0;
  boundary_dof(4,3) = 0.1;
  boundary_dof(5,3) = 0.1;
  boundary_dof(6,3) = -0.1;
  boundary_dof(7,3) = -0.1;
  boundary_dof(8,3) = 0.1;
  boundary_dof(9,3) = 0;

end