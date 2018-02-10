function [R_true,t_true,idx]=TruePose_RANSAC(R,t,F1,F2,K)

%%%%%%%% input %%%%%%%%%%%%%%%%%%%%%%
% F1&F2: 2D features, 3X1 vector
% K: intrinsic matrix

R_fall = [1   0   0
         0   -1   0
         0    0   -1];

incorrect = 0;
for i=1:size(R,3)
    b1 = F1;
    b2 = F2;
    b1 = b1/norm(b1);
    b2 = b2/norm(b2);
    d = -t(:,:,i)\[R(:,:,i)*b1 -b2];
    if d(1)>0 && d(2)>0 
        if R(:,:,i)'~=R_fall
        R_true = R(:,:,i);
        t_true = t(:,:,i);
        idx = i;
        end
    else
        incorrect = incorrect + 1;
    end
end

if incorrect == 4
    R_true = zeros(3);
    t_true = zeros(3,1);
    idx = 0;
end