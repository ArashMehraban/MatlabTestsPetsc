function boundary_dof = bc_bend(boundary_dof, boundary_vtx, load_increment, AppCtx)
  
  %boundary_dof(:,1) = -0.3*load_increment;  
  boundary_dof(:,2) = -2*load_increment; % 50% the size of the beam in the x-coords bent in -Y direction (beam_8e_l999_r998_6ss.png)  
  %boundary_dof(:,3) = -0.3*load_increment;

end