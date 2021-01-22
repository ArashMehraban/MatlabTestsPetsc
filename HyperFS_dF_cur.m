function f = HyperFS_dF_cur(dlta_ue,ddu, stored ,dXdx, wdetj, phys)
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
   fs = 0*ddu; %This dvdX in libCEED
   J = 0*wdetj; %Allocate space for J = det(F)
   [r,c] = size(ddu);
   blk = r/c; 
   idx = reshape(reshape(1:r,blk,c)',[],1); 
   permuted_ddu = ddu(idx,:);

   permuted_dXdx = dXdx(idx,:);
   
   %[B, D, W] = get_shape(2, 3);
   
   for i = 1:blk
     ddu_tmp= permuted_dXdx((i-1)*c+1:i*c,:) * permuted_ddu((i-1)*c+1:i*c,:);
     graddu(idx((i-1)*c+1:i*c),:) = ddu_tmp;
     F = eye(3) + stored((i-1)*c+1:i*c,:); %stored contains gradu from userf
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
     %fP is dvdX meaning dXdx^T *dP * wdetj
     fs(idx((i-1)*c+1:i*c),:) = (permuted_dXdx((i-1)*c+1:i*c,:)'* dstress * wdetj(i))/J(i);
   end
   
    fu1 = 0*dlta_ue(:,1);
    fu2 = 0*dlta_ue(:,2);
    fu3 = 0*dlta_ue(:,3);
   
   f = [  fu1      fu2      fu3; fs];    
end





