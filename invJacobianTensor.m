function [dets, invJe] = invJacobianTensor(elemJac, P)

 [r,c] = size(elemJac);
 invJe = zeros(r,c);
 blk = P*P*P;
 dets = zeros(blk,1);
 
 idx = reshape(reshape(1:r,blk,c)',[],1);
 
 permutedElemJac = elemJac(idx,:);
 
 %make 2D array to 3D:
%  sub = [1 2 3; 4 5 6; 7 8 9]
%  m = repmat(sub,8,1)
 %permute(reshape(permute(m,[2,1,3]),3,3,[]),[2,1,3])
  
 for i = 1:blk
     tmp=permutedElemJac((i-1)*c+1:i*c,:);
     invJe(idx((i-1)*c+1:i*c),:) =inv(tmp);
     dets(i) = det(tmp);
 end
 
end