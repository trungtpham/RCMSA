function P = homography_fit(X)

x1 = X(1:3,:);
x2 = X(4:6,:);

% ALT
H = vgg_H_from_x_lin(x1,x2);

% % nonlinear fit
% H = vgg_H_from_x_nonlin(Hlin,x1,x2);

% % Denormalise
% H = T2\H*T1;

P = H(:);


end


% function P = homography_fit(X)
% % Homography fitting function using Direct Linear Transformation (DLT) (see Hartley and Zisserman).
% % X is a 5 by 4 matrix.
% % X(1,i) is the x coordinate in the 1st image of the i-th correspondence.
% % X(2,i) is the y coordinate in the 1st image of the i-th correspondence.
% % X(3,i) is the x coordinate in the 2nd image of the i-th correspondence.
% % X(4,i) is the y coordinate in the 2nd image of the i-th correspondence.
% % X(5,i) is max(width_of_image,height_of_image) (for increasing stability via normalisation).
% % P is the resulting homography relation (a matrix) reshaped to a vector.
% 
% A = zeros(8,9);
% 
% for i=1:4
%     x1 = X(1,i)/X(5,i);
%     y1 = X(2,i)/X(5,i);
%     x2 = X(3,i)/X(5,i);
%     y2 = X(4,i)/X(5,i);
%     A((i-1)*2+1,:) = [ 0 0 0 -x2 -y2 -1 y1*x2 y1*y2 y1 ]; 
%     A((i-1)*2+2,:) = [ x2 y2 1 0 0 0 -x1*x2 -x1*y2 -x1 ];
% end    
%     
% [ u s v ] = svd(A);
% 
% H = reshape(v(:,end),3,3)';
% 
% N = [ 1/X(5,1) 0 0 ; 0 1/X(5,1) 0 ; 0 0 1 ];
% 
% H = N\(H*N); % Equivalent to H = inv(N)*H*N;
% 
% P = reshape(H',9,1);
% 
% end
