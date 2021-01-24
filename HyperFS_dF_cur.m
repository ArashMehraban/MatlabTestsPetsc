function f = HyperFS_dF_cur(dlta_ue,ddu, stored ,dXdx, wdetj, phys)
%HyperFS_dF_ref implements the derivative of constitutive model for the 
%compressible Neo-Hookean Hyperelasticity at finite strain using the
%current configuration
%
%IMPORTANT:
% Variable stored contains du from HyperFS_cur with following layout
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
%   du, ddu and dXdx have the same layout
%

   nu = phys.nu;
   E = phys.E;
   TwoMu=E/(1+nu);
   muu = TwoMu /2;
   Kbulk=E/(3*(1-2*nu));
   lambda=(3*Kbulk-TwoMu)/3;
    
   I3 = eye(3);
   
   graddu=  0*ddu; %Allocate sapce for graddu computed here
   fds = 0*ddu;    %fds is d_dvdX = dXdx_Transpose * dP * wdetj (libCEED notation)
   J = 0*wdetj;    %Allocate space for J = det(F)
   [r,c] = size(ddu);
   blk = r/c; 
   
   for i = 1:blk
     ddu_tmp= dXdx((i-1)*c+1:i*c,:) * ddu((i-1)*c+1:i*c,:);
     graddu((i-1)*c+1:i*c,:) = ddu_tmp;
     F = eye(3) + stored((i-1)*c+1:i*c,:);  % gradu computed in HyperFS_cur
     delta_b = (ddu_tmp * F' + F * ddu_tmp');
     J(i) = det(F);
     b = F * F';
     F_inv = F\I3; 
     tau = muu * b - (muu - 2*lambda*log(J(i)))*eye(3);
     %dtau = mu*delta_b + 2lambda * F^(-T) dF (Note: dF = ddu_tmp)
     Finv_contract_deltaF = sum(sum(F_inv' .* ddu_tmp));
     %dtau = muu *delta_b + 2* lambda * F_inv' * ddu_tmp;
     dtau = muu *delta_b + 2* lambda * Finv_contract_deltaF * I3;
     dstress = -ddu_tmp * F_inv * tau + dtau;
     fds((i-1)*c+1:i*c,:) = (dXdx((i-1)*c+1:i*c,:)'* dstress * wdetj(i))/J(i);
   end
   
    fu1 = 0*dlta_ue(:,1);
    fu2 = 0*dlta_ue(:,2);
    fu3 = 0*dlta_ue(:,3);
   
   f = [  fu1      fu2      fu3; fds];    
end





