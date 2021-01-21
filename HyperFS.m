function [usrfStored, f1] = HyperFS(du, dXdx, wdetj, phys)
%USERF_3d_ELAS provides weak form of the linear 3D Elastisity problem to solve 
%
% IMPORTANT:
%
%              [du1/dx | du2/dx | du3/dx]    <-for quadrature 1                           
%              [.....  | .....  |  .....]                             
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%           ---------------------------------
%              [du1/dy | du2/dy | du3/dy]    <-for quadrature 1
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]                     
%      du =    [.....  | .....  |  .....]                        
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%           ---------------------------------
%              [du1/dz | du2/dz | du3/dz]    <-for quadrature 1
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              
    
   nu = phys.nu;
   E = phys.E;
   TwoMu=E/(1+nu);
   muu = TwoMu /2;
   Kbulk=E/(3*(1-2*nu));
   lambda=(3*Kbulk-TwoMu)/3;

   gradu=  0*du; %Allocate sapce for gradu that gets computed here
   usrfStored = 0*du;
   f1 = 0*du; %This dvdX in libCEED
   J = 0*wdetj; %Allocate space for J = det(F)
   [r,c] = size(du);
   blk = r/c; 
   idx = reshape(reshape(1:r,blk,c)',[],1); 
   permuted_du = du(idx,:);
   permuted_dXdx = dXdx(idx,:);
   
   for i = 1:blk
     du_tmp= permuted_dXdx((i-1)*c+1:i*c,:) * permuted_du((i-1)*c+1:i*c,:);
     gradu(idx((i-1)*c+1:i*c),:) = du_tmp;
     usrfStored((i-1)*c+1:i*c,:) = du_tmp; % store for HyperFS_dF
     F = eye(3) + du_tmp;
     J(i) = det(F);
     C = F' * F;
     C_inv = C\eye(3);
     S = muu*eye(3) + (lambda*log(J(i))-muu)*C_inv;
     P = F * S;
     %f1 is dvdX meaning P * dXdx^T
     f1(idx((i-1)*c+1:i*c),:) = permuted_dXdx((i-1)*c+1:i*c,:)'* P * wdetj(i);
   end
       
end
