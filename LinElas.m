function [usrfStored, f1] = LinElas(du, dXdx, wdetj, phys)
%HyperFS_ref implements constitutive model for the compressible Neo-Hookean
%Hyperelasticity at finite strain using the reference configuration 
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

   usrfStored = 'dummy';
   gradu=  0*du; %Allocate sapce for gradu that gets computed here
   graduT= 0*du; %Allocate space for Transpose of gradu
   sigma = 0*du; %Allocate sapce for sigma computation
   f1 = 0*du; %This dvdX in libCEED
   [r,c] = size(du);
   blk = r/c; 
%    idx = reshape(reshape(1:r,blk,c)',[],1); 
%    permuted_du = du(idx,:);
%    permuted_dXdx = dXdx(idx,:);
   
   for i = 1:blk
%      tmp= permuted_dXdx((i-1)*c+1:i*c,:) * permuted_du((i-1)*c+1:i*c,:);
       tmp= dXdx((i-1)*c+1:i*c,:) * du((i-1)*c+1:i*c,:);
%      gradu(idx((i-1)*c+1:i*c),:) = tmp;
       gradu((i-1)*c+1:i*c,:) = tmp;
%      graduT(idx((i-1)*c+1:i*c),:) = tmp';
       graduT((i-1)*c+1:i*c,:) = tmp';
%      permuted_dXdxT((i-1)*c+1:i*c,:) = permuted_dXdx((i-1)*c+1:i*c,:)';
       dXdxT((i-1)*c+1:i*c,:) = dXdx((i-1)*c+1:i*c,:)';
   end
   
   e = 0.5 * (gradu+graduT); %strain
   
   nu = phys.nu;
   E = phys.E;
   ss = E / ((1 + nu)*(1 - 2*nu));
   
   for i=1:blk
%        tmp = e(idx((i-1)*c+1:i*c),:);
       tmp = e((i-1)*c+1:i*c,:);
       sigma11 = ss*((1 - nu)*tmp(1,1) + nu*tmp(2,2) + nu*tmp(3,3)); 
       sigma22 = ss*(nu*tmp(1,1)+ (1 - nu)*tmp(2,2) + nu*tmp(3,3));
       sigma33 = ss*(nu*tmp(1,1) + nu*tmp(2,2) + (1 - nu)*tmp(3,3));
       sigma23 = ss*(1 - 2*nu)*tmp(2,3)*0.5;
       sigma13 = ss*(1 - 2*nu)*tmp(1,3)*0.5;
       sigma12 = ss*(1 - 2*nu)*tmp(1,2)*0.5;
       sigma3x3 = [sigma11 sigma12 sigma13;
                   sigma12 sigma22 sigma23;
                   sigma13 sigma23 sigma33];               
%        sigma(idx((i-1)*c+1:i*c),:) = sigma3x3;
       sigma((i-1)*c+1:i*c,:) = sigma3x3;
   end
   
   for i = 1:blk
%      tmp= permuted_dXdxT((i-1)*c+1:i*c,:) * sigma(idx((i-1)*c+1:i*c),:);
     tmp= dXdxT((i-1)*c+1:i*c,:) * sigma((i-1)*c+1:i*c,:);
%      f1(idx((i-1)*c+1:i*c),:) = tmp*wdetj(i);
     f1((i-1)*c+1:i*c,:) = tmp*wdetj(i);
   end
   
end
