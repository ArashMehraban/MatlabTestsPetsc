function [usrfStored, f1] = LinElasF(du, dXdx, wdetj, phys)
%LinElasF implements constitutive model for the compressible Linear
%Elasticity 
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

   usrfStored = 'dummy'; %no storage of gradu is needed for linear problems
   sigma = 0*du;         %Allocate sapce for sigma computation
   f1 = 0*du;            %f1 is dvdX = dXdx_Transpose * sigma * wdetj (libCEED notation)
   [r,c] = size(du);
   blk = r/c; 
   
   nu = phys.nu;
   E = phys.E;
   ss = E / ((1 + nu)*(1 - 2*nu));
   
   for i=1:blk
       gradu = dXdx((i-1)*c+1:i*c,:) * du((i-1)*c+1:i*c,:); %gradu is 3x3
       strain = 0.5 * (gradu+gradu');
       sigma11 = ss*((1 - nu)*strain(1,1) + nu*strain(2,2) + nu*strain(3,3)); 
       sigma22 = ss*(nu*strain(1,1) + (1 - nu)*strain(2,2) + nu*strain(3,3));
       sigma33 = ss*(nu*strain(1,1) + nu*strain(2,2) + (1 - nu)*strain(3,3));
       sigma23 = ss*(1 - 2*nu)*strain(2,3)*0.5;
       sigma13 = ss*(1 - 2*nu)*strain(1,3)*0.5;
       sigma12 = ss*(1 - 2*nu)*strain(1,2)*0.5;
       sigma((i-1)*c+1:i*c,:) = [sigma11 sigma12 sigma13;
                                 sigma12 sigma22 sigma23;
                                 sigma13 sigma23 sigma33];               
       f1((i-1)*c+1:i*c,:) = dXdx((i-1)*c+1:i*c,:)' * sigma((i-1)*c+1:i*c,:) * wdetj(i);
   end
   
end
