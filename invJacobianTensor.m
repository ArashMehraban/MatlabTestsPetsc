function [dets, invJe] = invJacobianTensor(dx)

 [r,c] = size(dx);
 invJe = zeros(r,c);
 blk = r/c;
 dets = zeros(blk,1);

 for i = 1:blk
     tmp=dx((i-1)*c+1:i*c,:);
     invJe((i-1)*c+1:i*c,:) =inv(tmp);
     dets(i) = det(tmp);
 end
 
end