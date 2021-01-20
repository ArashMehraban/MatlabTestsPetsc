function f = HyperFS_dF(dlta_ue, xe, grad_dlta_ue, stored ,dXdx, wdetj, phys)
%USERF_3d_ELAS provides weak form of the linear 3D Elastisity problem to solve 
%
%IMPORTANT:
%     
%              [du1/dx | du2/dx | du3/dx]    <-for Node 1                            [d_dlta_u1/dx | d_dlta_u2/dx | d_dlta_u3/dx]    <-for Node 1
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ] 
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%           ---------------------------------                                     ----------------------------------------------------
%              [du1/dy | du2/dy | du3/dy]    <-for Node 1                            [d_dlta_u1/dy | d_dlta_u2/dy | d_dlta_u3/dy]   <-for Node 1
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]                    
% grad_ue =    [.....  | .....  |  .....]                            grad_dlta_ue =  [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                            ( Variation     [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                             of grad_ue )   [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%           ---------------------------------                                     ----------------------------------------------------
%              [du1/dz | du2/dz | du3/dz]    <-for Node 1                            [d_dlta_u1/dz | d_dlta_u2/dz | d_dlta_u3/dz]   <-for Node 1
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%
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
%              F =  eye(3) + grad_ue(i:num_row:end,:)';            <--- This truley F = I + nalba_u  
%              dltaF =  grad_dlta_ue(i:num_row:end,:)';   <--- This truley dltaF = I + nalba_dlta_u
%              ....
%              ...
%          end
             
 
    % Young's modulus 
    E=1;
    % Poisson ratio 
    nu=0.3;
    
        
    num_row = size(dlta_ue,1);
    
    TwoMu=E/(1+nu);
    muu = TwoMu /2;
    Kbulk=E/(3*(1-2*nu));
    lambda=(3*Kbulk-TwoMu)/3;
    
    I3 = eye(3);
    f1 = 0*grad_dlta_ue;

    for k = 1:num_row
        F =  I3 + grad_ue(k:num_row:end,:)';
        J = det(F);
        C = F'*F;
        C_inv  = C\I3; % C_inverse
        S = muu*I3 + (lambda*log(J)-muu)*C_inv; %(2nd Piola)
        F_inv = F\I3;       
        
        szf = size(F,1);        
              
        dPdF = zeros(szf,szf,szf,szf);
        for i = 1:szf
            for I = 1:szf
                for a = 1:szf
                    for A = 1:szf
                        dPdF(i,I,a,A) = I3(a,i)*S(A,I) + lambda*F_inv(A,a)*F_inv(I,i) - (lambda*log(J)-muu)*(F_inv(A,i)*F_inv(I,a) + I3(a,i)*C_inv(A,I));
                    end
                end
            end
        end
        
        dltaF = grad_dlta_ue(k:num_row:end,:)';
        dpdf = reshape(dPdF,9,9);
        
        mat_dp = dpdf*dltaF(:);
        
% %         dltaF = grad_dlta_ue(k:num_row:end,:)';
% %         dP = zeros(3,3);
% %         %Note: A, and a are dummy indecies
% %         for i = 1:szf
% %             for I = 1:szf
% %                 for a = 1:szf
% %                     for A = 1:szf
% %                         dP(i,I) = dP(i,I) + dPdF(i,I,a,A) * dltaF(a,A);
% %                     end
% %                 end
% %             end
% %         end
% %         
% %         f1(k:num_row:end,:) = dP;



       f1(k:num_row:end,:) = [mat_dp(1),mat_dp(2),mat_dp(3);
                              mat_dp(4),mat_dp(5),mat_dp(6);
                              mat_dp(7),mat_dp(8),mat_dp(9)] ;
    end

   
        %0*dlta_ue is a placeholder. Will not affect anything!
    f = [0*dlta_ue; f1];
        
end








        
%         dPdF = reshape(dPdF,szf^2,szf^2);
%         D = [dPdF(1,:);
%              dPdF(4,:);
%              dPdF(7,:);
%              dPdF(2,:);
%              dPdF(5,:);
%              dPdF(8,:);
%              dPdF(3,:);
%              dPdF(6,:);
%              dPdF(9,:)];
%                                                 %
%                                                 %  D computed aboves with DIJ values
%                                                 %                           [D11, D12, D13, D14, D15, D16, D17, D18 ,D19]
%                                                 %                           [D21, D22, D23, D24, D25, D26, D27, D28 ,D29]
%                                                 %                           [D31, D32, D33, D34, D35, D36, D37, D38 ,D39]
%                                                 %                           [D41, D42, D43, D44, D45, D46, D47, D48 ,D49]
%                                                 %                       D = [D51, D52, D53, D54, D55, D56, D57, D58 ,D59]
%                                                 %                           [D61, D62, D63, D64, D65, D66, D67, D68 ,D69]
%                                                 %                           [D71, D72, D73, D74, D75, D76, D77, D78 ,D79]
%                                                 %                           [D81, D82, D83, D84, D85, D86, D87, D88 ,D89]
%                                                 %                           [D91, D92, D93, D94, D95, D96, D97, D98 ,D99]
%                                                 %                                                                        9x9 matrix
%                                                 %
%         
%              
%                                                 %                   [d_dlta_u1/dx | d_dlta_u1/dy | d_dlta_u1/dz ]
%                                                 %  dltaF =          [d_dlta_u2/dx | d_dlta_u2/dy | d_dlta_u2/dz ]
%                                                 %                   [d_dlta_u3/dx | d_dlta_u3/dy | d_dlta_u3/dz ]
%                                                 %                                                               3x3 matrix
%                                                 
%         
%                                     
%         dltaF = grad_dlta_ue(k:num_row:end,:)';
% %         dltaJ = det(dltaF);
% %         dltaC = dltaF'*dltaF;
% %         dltaC_inv  = dltaC\I3; % dltaC_inverse
% %         dltaS = muu*I3 + (lambda*log(dltaJ)-muu)*dltaC_inv; %(2nd Piola)
% %         dltaP = dltaF*dltaS;            
%                                                 
%         f1(i:num_row:end,:) = reshape(D*dltaF(:),3,3); 




