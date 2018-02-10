clear all
clc

NumOfPoses = 2;
% Random pose generator


% To ensure f.o.v = 60 deg.
fov = 60*pi/180;
Rfeat = 1; Rcam = Rfeat/sin(fov/2);

W_T_C1 = [];
AbsolutePoses_true = zeros(3,4,NumOfPoses);
for k = 1:NumOfPoses
    tt = 2*pi*rand();
    WRC = [-sin(tt) 0 -cos(tt); cos(tt) 0 -sin(tt); 0 -1 0];
    WpC = Rcam*[cos(tt);sin(tt);0];
    if k == 1
        W_T_C1 = [WRC WpC];
        AbsolutePoses_true(:,:,k) = [eye(3) zeros(3,1)];
    else
        C_T_W = [WRC' -WRC'*WpC];
        AbsolutePoses_true(:,:,k) = C_T_W*[W_T_C1;zeros(1,3) 1];
    end
end
R = AbsolutePoses_true(:,1:3,2);
t = AbsolutePoses_true(:,4,2);

    T=[0 -t(3) t(2);  
        t(3) 0 -t(1);  
        -t(2) t(1) 0  
    ];  
      
    E=T*R  
    [U,S,V]=svd(E);  
    disp('S?=?diag(1,1,0)')  
    S  
    W=[0 -1 0;  
        1 0 0;  
        0 0 1  
    ];  
    P1=[U*W*V' U(:,3)]  
    P2=[U*W'*V' U(:,3)] 
    
    R1 = U*W*V';
    R2 = U*W'*V';
    if(det(R1)<0)
        R1 = -R1;
    end
    if(det(R2)<0)
        R2 = -R2;
    end
    disp('check R..')  
    norm(R1-R)  
    norm(R2-R)  
      
    disp('check t..')  
    norm(U(:,3) - t/norm(t))  
    norm(- U(:,3) - t/norm(t))  