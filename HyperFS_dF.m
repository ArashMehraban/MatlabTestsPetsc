function f = HyperFS_dF(dlta_ue,ddu, stored ,dXdx, wdetj, phys)
%USERF_3d_ELAS provides weak form of the linear 3D Elastisity problem to solve 
%
%IMPORTANT:
%     
%              [du1/dx | du2/dx | du3/dx]    <-for quadrature 1                      [d_dlta_u1/dx | d_dlta_u2/dx | d_dlta_u3/dx]    <-for quadrature 1
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ] 
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%           ---------------------------------                                     ----------------------------------------------------
%              [du1/dy | du2/dy | du3/dy]    <-for quadrature 1                      [d_dlta_u1/dy | d_dlta_u2/dy | d_dlta_u3/dy]   <-for quadrature 1
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]                    
% grad_ue =    [.....  | .....  |  .....]                            grad_dlta_ue =  [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                            ( Variation     [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                             of grad_ue )   [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%           ---------------------------------                                     ----------------------------------------------------
%              [du1/dz | du2/dz | du3/dz]    <-for quadrature 1                      [d_dlta_u1/dz | d_dlta_u2/dz | d_dlta_u3/dz]   <-for quadrature 1
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%              [.....  | .....  |  .....]                                            [   .....     |    .....     |     .....   ]
%
   nu = phys.nu;
   E = phys.E;
   TwoMu=E/(1+nu);
   muu = TwoMu /2;
   Kbulk=E/(3*(1-2*nu));
   lambda=(3*Kbulk-TwoMu)/3;
    
   I3 = eye(3);
   
   graddu=  0*ddu; %Allocate sapce for gradu that gets computed here
   fP = 0*ddu; %This dvdX in libCEED
   J = 0*wdetj; %Allocate space for J = det(F)
   [r,c] = size(ddu);
   blk = r/c; 

   
   for i = 1:blk
     ddu_tmp= dXdx((i-1)*c+1:i*c,:) * ddu((i-1)*c+1:i*c,:);
     graddu((i-1)*c+1:i*c,:) = ddu_tmp;
     F = eye(3) + stored((i-1)*c+1:i*c,:); %stored contains gradu from userf
     deltaE = 0.5*(ddu_tmp' * F + F' * ddu_tmp);
     J(i) = det(F);
     C = F' * F;
     C_inv = C\I3;
     S = muu*I3 + (lambda*log(J(i))-muu)*C_inv;
     %dS = dSdE:deltaE
     %   = lambda(Cinv:deltaE)Cinv + 2(mu-lambda*log(J))Cinv*deltaE*Cinv
     %C_inv : deltaE (tensor contraction which is a pointwise multiplication of C_inv and deltaE and then add all teh elements of the resulting matrix)
     Cinv_contract_deltaE = sum(sum(C_inv .* deltaE)); 
     dS = lambda * Cinv_contract_deltaE * C_inv + 2*(muu-lambda*log(J(i)))*C_inv*deltaE*C_inv;
     %dP = deltaF*S + F*dS where deltaF = ddu_tmp
     dP = ddu_tmp * S + F*dS;
     %fP is dvdX meaning dXdx^T *dP * wdetj
     fP((i-1)*c+1:i*c,:) = dXdx((i-1)*c+1:i*c,:)'* dP * wdetj(i);
   end
   
    fu1 = 0*dlta_ue(:,1);
    fu2 = 0*dlta_ue(:,2);
    fu3 = 0*dlta_ue(:,3);
   
   f = [  fu1      fu2      fu3; fP];    
end





