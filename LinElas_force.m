function f0 = LinElas_force(xe,phys)

      % RHS:
      g1 = 0*xe(:,1);
      g2 = 0*xe(:,1);
      g3 = 0*xe(:,1);
      
      f0 = [-g1, -g2 , -g3]; 
    
end