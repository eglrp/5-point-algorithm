function [R_true,t_true,idx]=TruePose(R,t,F1,F2,K)

%%%%%%%% input %%%%%%%%%%%%%%%%%%%%%%
% F1&F2: 2D features, 3X1 vector
% K: intrinsic matrix

incorrect = 0;
for i=1:size(R,3)
    b1 = F1;
    b2 = F2;
    b1 = b1/norm(b1);
    b2 = b2/norm(b2);
    A = [R(:,:,i)*b1 -b2 t(:,:,i)];
    [U,S,V] = svd(A);
    d=V(:,end)/V(end,end);   
    if d(1)>0 && d(2)>0 
        R_true = R(:,:,i);
        t_true = t(:,:,i);
        idx = i;
    else
        incorrect = incorrect + 1;
    end
end

if incorrect == 4
    R_true = zeros(3);
    t_true = zeros(3,1);
    idx = 0;
end
