function [f, usrfStored] = LinElasF_dF(dlta_ue,ddu, stored ,dXdx, wdetj, phys)
%LinElasF_dF implements derivative of the constitutive model for the 
%compressible Linear Elasticity 
%
% IMPORTANT:
%
%              [ddu1/dx | ddu2/dx | ddu3/dx]    <-for quadrature 1 
%              [ddu1/dy | ddu2/dy | ddu3/dy]    <-for quadrature 1
%              [ddu1/dz | ddu2/dz | ddu3/dz]    <-for quadrature 1
%            ------------------------------
%              [ddu1/dx | ddu2/dx | ddu3/dx]    <-for quadrature 2 
%              [ddu1/dy | ddu2/dy | ddu3/dy]    <-for quadrature 2
%              [ddu1/dz | ddu2/dz | ddu3/dz]    <-for quadrature 2
%            ----------------------------------
%              [......  | ......  |  ......]           .
%              [......  | ......  |  ......]           .
%              [......  | ......  |  ......]           .
%            ----------------------------------
%              [......  | ......  |  ......]           .
%              [......  | ......  |  ......]           .                  
%     ddu =    [......  | ......  |  ......]           .
%            -----------------------------------
%              [......  | ......  |  ......]           .
%              [......  | ......  |  ......]           .
%              [......  | ......  |  ......]           .
%            ----------------------------------
%              [......  | ......  |  ......]           .
%              [......  | ......  |  ......]           .
%              [......  | ......  |  ......]           .
%            ----------------------------------
%              [......  | ......  |  ......]           .
%              [......  | ......  |  ......]           .
%              [......  | ......  |  ......]           .
%            ----------------------------------
%              [......  | ......  |  ......]           .
%              [......  | ......  |  ......]           .
%              [......  | ......  |  ......]           .
%   
%
%   ddu and dXdx have the same layout
%

   usrfStored = 'dummy'; %no storage of gradu is needed for linear problem
   dsigma = 0*ddu;       %Allocate sapce for sigma computation
   fdsigma = 0*ddu;      %f1 is dvdX = d_dXdx_Transpose * sigma * wdetj (libCEED notation)
   [r,c] = size(ddu);
   blk = r/c; 
   
   nu = phys.nu;
   E = phys.E;
   ss = E / ((1 + nu)*(1 - 2*nu));
   
   for i=1:blk
       graddu = dXdx((i-1)*c+1:i*c,:) * ddu((i-1)*c+1:i*c,:); %graddu is 3x3
       strain = 0.5 * (graddu+graddu');
       dsigma11 = ss*((1 - nu)*strain(1,1) + nu*strain(2,2) + nu*strain(3,3)); 
       dsigma22 = ss*(nu*strain(1,1) + (1 - nu)*strain(2,2) + nu*strain(3,3));
       dsigma33 = ss*(nu*strain(1,1) + nu*strain(2,2) + (1 - nu)*strain(3,3));
       dsigma23 = ss*(1 - 2*nu)*strain(2,3)*0.5;
       dsigma13 = ss*(1 - 2*nu)*strain(1,3)*0.5;
       dsigma12 = ss*(1 - 2*nu)*strain(1,2)*0.5;
       dsigma((i-1)*c+1:i*c,:) = [dsigma11 dsigma12 dsigma13;
                                 dsigma12 dsigma22 dsigma23;
                                 dsigma13 dsigma23 dsigma33];               
       fdsigma((i-1)*c+1:i*c,:) = dXdx((i-1)*c+1:i*c,:)' * dsigma((i-1)*c+1:i*c,:) * wdetj(i);
   end
   
   fu1 = 0*dlta_ue(:,1);
   fu2 = 0*dlta_ue(:,2);
   fu3 = 0*dlta_ue(:,3);
   
   f = [fu1, fu2,fu3; fdsigma];   
end
