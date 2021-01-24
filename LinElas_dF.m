function f = LinElas_dF(dlta_ue,ddu, stored ,dXdx, wdetj, phys)
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

   %ddu is variation of du
   graddu=  0*ddu; %Allocate sapce for gradu that gets computed here
   gradduT= 0*ddu; %Allocate space for Transpose of gradu
   dsigma = 0*ddu; %Allocate sapce for sigma computation
   f_dsigma = 0*ddu; %This dvdX in libCEED
   [r,c] = size(ddu);
   blk = r/c; 
%    idx = reshape(reshape(1:r,blk,c)',[],1); 
%    permuted_ddu = ddu(idx,:);
%    permuted_dXdx = dXdx(idx,:);
   
   for i = 1:blk
%      tmp= permuted_dXdx((i-1)*c+1:i*c,:) * permuted_ddu((i-1)*c+1:i*c,:);
     tmp= dXdx((i-1)*c+1:i*c,:) * ddu((i-1)*c+1:i*c,:);
%      graddu(idx((i-1)*c+1:i*c),:) = tmp;
     graddu((i-1)*c+1:i*c,:) = tmp;
%      gradduT(idx((i-1)*c+1:i*c),:) = tmp';
     gradduT((i-1)*c+1:i*c,:) = tmp';
%      permuted_dXdxT((i-1)*c+1:i*c,:) = permuted_dXdx((i-1)*c+1:i*c,:)';
     dXdxT((i-1)*c+1:i*c,:) = dXdx((i-1)*c+1:i*c,:)';
   end
   
   e = 0.5 * (graddu+gradduT); %strain
   
   nu = phys.nu;
   E = phys.E;
   ss = E / ((1 + nu)*(1 - 2*nu));
   
   for i=1:blk
%        tmp = e(idx((i-1)*c+1:i*c),:);
       tmp = e((i-1)*c+1:i*c,:);
       dsigma11 = ss*((1 - nu)*tmp(1,1) + nu*tmp(2,2) + nu*tmp(3,3)); 
       dsigma22 = ss*(nu*tmp(1,1)+ (1 - nu)*tmp(2,2) + nu*tmp(3,3));
       dsigma33 = ss*(nu*tmp(1,1) + nu*tmp(2,2) + (1 - nu)*tmp(3,3));
       dsigma23 = ss*(1 - 2*nu)*tmp(2,3)*0.5;
       dsigma13 = ss*(1 - 2*nu)*tmp(1,3)*0.5;
       dsigma12 = ss*(1 - 2*nu)*tmp(1,2)*0.5;
       dsigma3x3 = [dsigma11 dsigma12 dsigma13;
                   dsigma12 dsigma22 dsigma23;
                   dsigma13 dsigma23 dsigma33];               
%        dsigma(idx((i-1)*c+1:i*c),:) = dsigma3x3;
       dsigma((i-1)*c+1:i*c,:) = dsigma3x3;
   end
   
   for i = 1:blk
%      tmp= permuted_dXdxT((i-1)*c+1:i*c,:) * dsigma(idx((i-1)*c+1:i*c),:);
     tmp= dXdxT((i-1)*c+1:i*c,:) * dsigma((i-1)*c+1:i*c,:);
%      f_dsigma(idx((i-1)*c+1:i*c),:) = tmp*wdetj(i);
     f_dsigma((i-1)*c+1:i*c,:) = tmp*wdetj(i);
   end
   
    fu1 = 0*dlta_ue(:,1);
    fu2 = 0*dlta_ue(:,2);
    fu3 = 0*dlta_ue(:,3);
   
   f = [  fu1      fu2      fu3; f_dsigma];

end