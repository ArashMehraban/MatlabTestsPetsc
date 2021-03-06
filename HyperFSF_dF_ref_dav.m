function f = HyperFSF_dF_ref_dav(dlta_ue,ddu, stored ,dXdx, wdetj, phys)
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
   
   Graddu=  0*ddu; %Allocate sapce for graddu computed here
   fdP = 0*ddu;    
   J = 0*wdetj;    %Allocate space for J = det(F)
   [r,c] = size(ddu);
   blk = r/c; 

   
   for i = 1:blk
     ddu_tmp= dXdx((i-1)*c+1:i*c,:) * ddu((i-1)*c+1:i*c,:);  %ref config
     
     %drF
     Graddu((i-1)*c+1:i*c,:) = ddu_tmp;  %ref config
     
     F = stored((i-1)*c+1:i*c,:); % F computed in HyperFS_ref_dav
     F_inv=F\I3;
     J(i) = det(F);
     b = F*F';
     tau = (muu*b + (lambda*log(J(i))-muu)*I3);
     
     % db = drF*F^{T} + F*drF^{T}
     delta_b = (ddu_tmp * F' + F * ddu_tmp'); % current config
     
     delta_F_inv = -F_inv*ddu_tmp*F_inv;
                                   %(b^(-1)  :  delta_b)
     inv_b_contract_delta_b = sum(sum(inv(b) .* delta_b)); 
     
     dtau = (muu*delta_b + 0.5*lambda*inv_b_contract_delta_b * I3);
     dstress = (tau * delta_F_inv') + dtau*F_inv';
     
     fdP((i-1)*c+1:i*c,:) =  dXdx((i-1)*c+1:i*c,:)'* dstress * wdetj(i);
   end
   
   

    fu1 = 0*dlta_ue(:,1);
    fu2 = 0*dlta_ue(:,2);
    fu3 = 0*dlta_ue(:,3);
   
   f = [  fu1      fu2      fu3; fdP];    
end

% % The issue here is we don't want to store F when evaluating the action of the Jacobian. We have
% % 
% % x_current = x_initial + u
% % 
% % We also have a constitutive relation \tau(\bm b, J) from Eq 14, and we're discretizing Eq 7 of Davydov, which can be written
% % 
% % 
% % \int_{B_current} \sigma : symgrad(du) = \int_{B_initial} J \sigma : symgrad(du)
% % 
% % where \tau = J \sigma and symgrad is in the current configuration. Their notation is clumsy, IMO, since we're really doing the integral in the current configuration so I'd just use \sigma (Cauchy stress). 
% % We know that
% % 
% % dF = \partial du/\partial x_initial
% %    = (\partial du/\partial x_current) (\partial x_current/\partial x_initial)
% %    = (\partial du/\partial x_current) F
% %    = grad(du) F
% % 
% % So you need to be able to evaluate db(du), which is
% % 
% %   db = dF F^T + F dF^T
% %      = grad(du) F F^T + F F^T grad(du)^T
% %      = grad(du) b + b grad(du)^T





