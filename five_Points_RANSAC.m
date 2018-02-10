clear all; close all;
addpath('../');
addpath('../robotics3D/');

for whole_iter = 1:10
NumOfFeatures = 50;
NumOfPoses = 2;

PoseGraphMatrix = randi([0 1], NumOfFeatures, NumOfPoses);
both_see_feat = find(sum(PoseGraphMatrix,2)==2);
% Make sure two cameras observe more than 6 features in common
if length(both_see_feat) < 10
    addMoreFeature = 10 - length(both_see_feat) + 1;
    listIdx = find(sum(PoseGraphMatrix,2) < 2);
    for k = 1:addMoreFeature
        PoseGraphMatrix(listIdx(k),:) = [1 1];
    end
end
both_see_feat = find(sum(PoseGraphMatrix,2)==2);

% maxIter = 30;
% whole_iter = 1;
% for whole_iter = 21:maxIter
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
% saveas(gcf,['/home/jahid/Fei/summer/5Points/test/iter',num2str(whole_iter),' 3D points.jpg']);

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
%     saveas(gcf,['/home/jahid/Fei/summer/5Points/test/iter',num2str(whole_iter),' 2D points camera',num2str(k),'.jpg']);
end

% add noise to feature
sigma_z = 0;
featureExtracted = cell(NumOfPoses, 1);
for k = 1:NumOfPoses
   featureExtracted{k} = featureExtracted_true{k} + sigma_z*randn(size(featureExtracted_true{k}));   
%     figure(4); clf;
%     plot(featureExtracted_true{k}(1,:), featureExtracted_true{k}(2,:), 'b*'); hold on;
%     plot(featureExtracted{k}(1,:), featureExtracted{k}(2,:), 'ro'); hold on;
%     pause;
end

% RANSAC

max_inliners = -1;
for iter = 1:1000
    % randomly select 6 common points, 5 use for computing essential matix,
    % 1 for folllowing test
    feat_select = both_see_feat(randperm(length(both_see_feat),6),:);
    feat_select_1 = PoseGraphMatrix(feat_select,1);
    feat_select_2 = PoseGraphMatrix(feat_select,2);
    % compute essential matrix candidates (up to 10)
%     Evec = calibrated_fivepoint(featureExtracted{1}(:,feat_select_1(1:5)),featureExtracted{2}(:,feat_select_2(1:5)));
    Evec = calibrated_fivepoint_GrLex(featureExtracted{1}(:,feat_select_1(1:5)),featureExtracted{2}(:,feat_select_2(1:5)));
    % check which one among all candidates is correct one
    % by q2'*E*q1 = 0 
    Q1 = featureExtracted{1}(:,feat_select_1(6));
    Q2 = featureExtracted{2}(:,feat_select_2(6));
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
      if ~isempty(Evec)
          E = reshape(Evec(:,Evec_idx),[3,3]);
      else
          E = eye(3);
      end
      % compute 4 candidates for R and t
      [R(:,:,1), t(:,:,1), R(:,:,2), t(:,:,2), R(:,:,3), t(:,:,3), R(:,:,4), t(:,:,4)] = CameraPose(E);
      [R_true,t_true,n]=TruePose(R,t,featureExtracted{1}(:,feat_select_1(6)),featureExtracted{2}(:,feat_select_2(6)),eye(3));
      [num_inliners,error] = count_inliner(R_true,t_true,featureExtracted,PoseGraphMatrix(both_see_feat,:),0.5);
      if num_inliners > max_inliners
          max_inliners = num_inliners;
          R_est = R_true;
          t_est = t_true;
      end
%       disp(['Iter: ',num2str(iter)]);      
end

% error_final(whole_iter) = acos((trace(AbsolutePoses_true(:,1:3,2)*R_est')-1)/2);
figure;
set(gcf,'visible','off');
Draw3DWorld([R_est,t_est], FeatureBag_true, W_T_C1);hold on
Draw3DWorld(AbsolutePoses_true(:,:,2), FeatureBag_true, W_T_C1);
saveas(gcf,['/home/jahid/Fei/summer/5Point_new/test3/iter',num2str(whole_iter),' camera pose.jpg']);
close all
disp(['Iter: ',num2str(whole_iter)]);
end