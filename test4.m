clear all; close all;
addpath('../');
addpath('../robotics3D/');


NumOfFeatures = 6;
NumOfPoses = 2;

PoseGraphMatrix = ones(NumOfFeatures, NumOfPoses);

maxIter = 100;
% whole_iter = 1;
for whole_iter = 1:maxIter

%% 3D-WORLD creation


% To ensure f.o.v = 60 deg.
fov = 60*pi/180;
Rfeat = 1; Rcam = Rfeat/sin(fov/2);


% Random pose generator
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


% Random feature generator
FeatureBag_true = zeros(4,NumOfFeatures);
%figure(2); hold on;
for k = 1:NumOfFeatures
    tt = 2*pi*rand();
    phi = pi*rand() - pi/2;
    FeatureBag_true(:,k) = [Rfeat*[cos(phi)*cos(tt);cos(phi)*sin(tt);sin(phi)];1]; % <--- w.r.t the world
    FeatureBag_true(1:3,k) = InversePose(W_T_C1)*[FeatureBag_true(1:3,k);1];
end
figure(1);
set(gcf,'visible','off');
plot3(FeatureBag_true(1,:),FeatureBag_true(2,:),FeatureBag_true(3,:),'*');title('3D points');
% saveas(gcf,['/home/jahid/Fei/summer/5Point_new/test/iter',num2str(whole_iter),' 3D points.jpg']);

% Create feature projections as measurement
featureExtracted_true = cell(NumOfPoses, 1);
for k = 1:NumOfPoses    
    numFeatPosek = find(PoseGraphMatrix(:,k));
    featureExtracted_true{k} = zeros(2, length(numFeatPosek));
    for l = 1:length(numFeatPosek)
        % assign correct feature index for pose graph matrix
        PoseGraphMatrix(numFeatPosek(l),k) = l; 
        
        % Pose Transformation
        Ck_T_W = AbsolutePoses_true(:,:,k);
        Cl_T_W = AbsolutePoses_true(:,:,FeatureBag_true(4,numFeatPosek(l)));
        Ck_T_Cl = Ck_T_W*[InversePose(Cl_T_W);zeros(1,3) 1];
        
        % Image Projection
        Ck_p_fl = Ck_T_Cl*[FeatureBag_true(1:3,numFeatPosek(l));1];
        featureExtracted_true{k}(:,l) = Ck_p_fl(1:2)./Ck_p_fl(3);
    end
    featureExtracted_true{k} = [featureExtracted_true{k};ones(1,size(featureExtracted_true{k},2))];
    figure(k+1);
    set(gcf,'visible','off');
    plot(featureExtracted_true{k}(1,:),featureExtracted_true{k}(2,:),'*');title(['2D points in camera',num2str(k)]);
%     saveas(gcf,['/home/jahid/Fei/summer/5Point_new/test/iter',num2str(whole_iter),' 2D points camera',num2str(k),'.jpg']);
end

% Evec = calibrated_fivepoint(featureExtracted_true{1}(:,1:5),featureExtracted_true{2}(:,1:5));
Evec = calibrated_fivepoint_GrLex(featureExtracted_true{1}(:,1:5),featureExtracted_true{2}(:,1:5));
Q1 = featureExtracted_true{1}(:,6);
Q2 = featureExtracted_true{2}(:,6);
Q1 = Q1';
Q2 = Q2';
Q = [Q2(:,1).*Q1(:,1) , ...
     Q2(:,2).*Q1(:,1) , ...
     Q2(:,3).*Q1(:,1) , ... 
     Q2(:,1).*Q1(:,2) , ...
     Q2(:,2).*Q1(:,2) , ...
     Q2(:,3).*Q1(:,2) , ...
     Q2(:,1).*Q1(:,3) , ...
     Q2(:,2).*Q1(:,3) , ...
     Q2(:,3).*Q1(:,3) ] ; 
Evec_error = zeros(size(Evec,2),1);
for i=1:size(Evec,2)
     Evec_error(i) = abs(Q*Evec(:,i));
end
[~,Evec_idx] = min(Evec_error);

% reshape essential matrix to be 3*3
E = reshape(Evec(:,Evec_idx),[3,3]);
[R(:,:,1), t(:,:,1), R(:,:,2), t(:,:,2), R(:,:,3), t(:,:,3), R(:,:,4), t(:,:,4)] = CameraPose(E);
% [R_true,C_true,n1]=TruePose(R,C,featureExtracted_true{1}(1:2,:),featureExtracted_true{2}(1:2,:),eye(3));

% [R_true_2,t_true_2,n2]=TruePose_RANSAC(R,t,featureExtracted_true{1}(:,6),featureExtracted_true{2}(:,6),eye(3));
[R_true_2,t_true_2,n2]=TruePose(R,t,featureExtracted_true{1}(:,6),featureExtracted_true{2}(:,6),eye(3));


figure;
% set(gcf,'visible','off');
Draw3DWorld([R_true_2,t_true_2], FeatureBag_true, W_T_C1);hold on
% text(C_true(1,1,Idx)/C_true(3,1,Idx),C_true(2,1,Idx)/C_true(3,1,Idx),'estimated');
Draw3DWorld(AbsolutePoses_true(:,:,2), FeatureBag_true, W_T_C1);
% saveas(gcf,['/home/jahid/Fei/summer/5Point_new/test/iter',num2str(whole_iter),' camera pose.jpg']);
% close all


%% display
disp(['Iter: ',num2str(whole_iter)]);
% disp(['n1: ',num2str(n1)]);
% disp(['n2: ',num2str(n2)]);
% disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
end