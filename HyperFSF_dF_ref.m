function f = HyperFSF_dF_ref(dlta_ue,ddu, stored ,dXdx, wdetj, phys)
%HyperFSF_dF_ref implements the derivative of constitutive model for the 
%compressible Neo-Hookean Hyperelasticity at finite strain using the
%reference configuration
%
%IMPORTANT:
% Variable stored contains du from HyperFS_ref with following layout
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
   fdP = 0*ddu;    %fdP is d_dvdX = dXdx_Transpose * dP * wdetj (libCEED notation)
   J = 0*wdetj;    %Allocate space for J = det(F)
   [r,c] = size(ddu);
   blk = r/c; 

   
   for i = 1:blk
     ddu_tmp= dXdx((i-1)*c+1:i*c,:) * ddu((i-1)*c+1:i*c,:);
     graddu((i-1)*c+1:i*c,:) = ddu_tmp;
     F = eye(3) + stored((i-1)*c+1:i*c,:); % gradu computed in HyperFS_ref
     deltaE = 0.5*(ddu_tmp' * F + F' * ddu_tmp);
     J(i) = det(F);
     C = F' * F;
     C_inv = C\I3;
     S = muu*I3 + (lambda*log(J(i))-muu)*C_inv;
     %dS = dSdE:deltaE
     %   = lambda(Cinv:deltaE)Cinv + 2(mu-lambda*log(J))Cinv*deltaE*Cinv
     %C_inv : deltaE (tensor contraction: First, pointwise multiply C_inv and deltaE, then, add all entries in the resulting matrix)
     Cinv_contract_deltaE = sum(sum(C_inv .* deltaE)); 
     dS = lambda * Cinv_contract_deltaE * C_inv + 2*(muu-lambda*log(J(i)))*C_inv*deltaE*C_inv;
     %dP = deltaF*S + F*dS where deltaF = ddu_tmp
     dP = ddu_tmp * S + F*dS;
     fdP((i-1)*c+1:i*c,:) = dXdx((i-1)*c+1:i*c,:)'* dP * wdetj(i);
   end
   
    fu1 = 0*dlta_ue(:,1);
    fu2 = 0*dlta_ue(:,2);
    fu3 = 0*dlta_ue(:,3);
   
   f = [  fu1      fu2      fu3; fdP];    
end





