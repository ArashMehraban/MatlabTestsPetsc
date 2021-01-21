function [dets, invJe] = invJacobianTensor(elemJac)

 [r,c] = size(elemJac);
 invJe = zeros(r,c);
 blk = r/c;
 dets = zeros(blk,1);
 
 idx = reshape(reshape(1:r,blk,c)',[],1);
 
 permutedElemJac = elemJac(idx,:);
  
 for i = 1:blk
     tmp=permutedElemJac((i-1)*c+1:i*c,:);
     invJe(idx((i-1)*c+1:i*c),:) =inv(tmp);
     dets(i) = det(tmp);
 end
 
end