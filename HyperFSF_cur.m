function [usrfStored, f1] = HyperFSF_cur(du, dXdx, wdetj, phys)
%HyperFSF_cur implements constitutive model for the compressible Neo-Hookean
%Hyperelasticity at finite strain using the current configuration 
%
% IMPORTANT:
%
%              [du1/dx | du2/dx | du3/dx]    <-for quadrature 1 
%              [du1/dy | du2/dy | du3/dy]    <-for quadrature 1
%              [du1/dz | du2/dz | du3/dz]    <-for quadrature 1
%            ------------------------------
%              [du1/dx | du2/dx | du3/dx]    <-for quadrature 2 
%              [du1/dy | du2/dy | du3/dy]    <-for quadrature 2
%              [du1/dz | du2/dz | du3/dz]    <-for quadrature 2
%            ------------------------------
%              [.....  | .....  |  .....]           .
%              [.....  | .....  |  .....]           .
%              [.....  | .....  |  .....]           .
%            ------------------------------
%              [.....  | .....  |  .....]           .
%              [.....  | .....  |  .....]           .                  
%      du =    [.....  | .....  |  .....]           .
%            ------------------------------
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%            ------------------------------
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%            ------------------------------
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%            ------------------------------
%              [.....  | .....  |  .....]          .
%              [.....  | .....  |  .....]          .
%              [.....  | .....  |  .....]          .
%   
%
%   du and dXdx have the same layout
%             
    
   nu = phys.nu;
   E = phys.E;
   TwoMu=E/(1+nu);
   muu = TwoMu /2;
   Kbulk=E/(3*(1-2*nu));
   lambda=(3*Kbulk-TwoMu)/3;

   gradu=  0*du;      %Allocate sapce for gradu computed here
   usrfStored = 0*du; %Allocate space for usrfStored as output for HyperFSF_dF_cur 
   f1 = 0*du;         %f1 is dvdX = dXdx_Transpose * P * wdetj (libCEED notation)
   J = 0*wdetj;       %Allocate space for J = det(F)
   [r,c] = size(du);
   blk = r/c; 
   
   for i = 1:blk
     du_tmp= dXdx((i-1)*c+1:i*c,:) * du((i-1)*c+1:i*c,:);
     gradu((i-1)*c+1:i*c,:) = du_tmp;
     usrfStored((i-1)*c+1:i*c,:) = du_tmp; 
     F = eye(3) + du_tmp;
     J(i) = det(F);
     b = F * F';
     tau = muu * b - (muu - 2*lambda*log(J(i)))*eye(3);
     f1((i-1)*c+1:i*c,:) = (dXdx((i-1)*c+1:i*c,:)'* tau * wdetj(i))/J(i);
   end
       
end
