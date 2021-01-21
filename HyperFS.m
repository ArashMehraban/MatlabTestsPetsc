function [usrfStored, f0,f1] = HyperFS(grad_ue, dXdx, wdetj, phys)
%USERF_3d_ELAS provides weak form of the linear 3D Elastisity problem to solve 
%
% IMPORTANT:
%
%              [du1/dx | du2/dx | du3/dx]    <-for Node 1                           
%              [.....  | .....  |  .....]                             
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%           ---------------------------------
%              [du1/dy | du2/dy | du3/dy]    <-for Node 1
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]                     
% grad_ue =    [.....  | .....  |  .....]                        
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%           ---------------------------------
%              [du1/dz | du2/dz | du3/dz]    <-for Node 1
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              [.....  | .....  |  .....]
%              
%     NOTE 1: 
%        We need to compute F = I + nabla_u. Notice that each row of
%        nabla_u is placed in every other 8th row in above grad_ue. 
%     NOTE 2:
%        The nabla_u we take in NOTE 1, is really the TRANSPOSE of the nabla_u 
%         we really need to compute F = I + nabla_u. So transpose it.
%      Summary of NOTE 1 and NOTE 2 gives:
%       Write all operations in this for loop:
%          
%          num_row = size(ue,1);       <-- num_row is 8 here 
%          for i = 1:num_row 
%              F =  eye(3) + grad_ue(i:num_row:end,:)';  <--- This truley F = I + nalba_u      
%              ....
%              ...
%          end
    
    % Young's modulus 
    E=1;
    % Poisson ratio 
    nu=0.3;

        
    num_row = size(ue,1);
    
    TwoMu=E/(1+nu);
    muu = TwoMu /2;
    Kbulk=E/(3*(1-2*nu));
    lambda=(3*Kbulk-TwoMu)/3;
    
    
    f1 = 0 * grad_ue; 

    for i = 1:num_row
        F =  eye(3) + grad_ue(i:num_row:end,:)';
        J = det(F);
        C = F'*F;
        C_inv  = C\eye(3); % C_inverse
        S = muu*eye(3) + (lambda*log(J)-muu)*C_inv; %(2nd Piola)
        P = F*S;     % P is in referential config (1st Piola)
        
        f1(i:num_row:end,:) = P';
    end
    
      % RHS:
      g1 = 0*ue(:,1);
      g2 = 0*ue(:,1);
      g3 = 0*ue(:,1);
      
    f0 = [-g1, -g2 , -g3]; 
    
end