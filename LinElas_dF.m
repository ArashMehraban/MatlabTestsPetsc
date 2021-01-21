function f = LinElas_dF(grad_dlta_ue, stored ,dXdx, wdetj, phys)
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


end