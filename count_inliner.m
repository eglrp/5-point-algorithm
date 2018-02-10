function [num_inliners,error] = count_inliner(R,t,featureExtracted_true,both_see_feat,thresh)

num_inliners = 0;
error = zeros(length(both_see_feat),1);
for i = 1:length(both_see_feat)
    b1 = featureExtracted_true{1}(:,both_see_feat(i,1));
    b2 = featureExtracted_true{2}(:,both_see_feat(i,2));
    b1 = b1/norm(b1);
    b2 = b2/norm(b2);
    A = [R*b1 -b2 t];
    [U,S,V] = svd(A);
    d=V(:,end)/V(end,end);
    Pf2 = d(1)*R*b1 + t;
    Pf2 = Pf2./Pf2(3,:);
    Pf1 = d(2)*R*b2 + t;
    Pf1 = Pf1./Pf1(3,:);
    error(i) = norm(Pf2-b2)+norm(Pf1-b1);
    if error(i) < thresh
        num_inliners = num_inliners + 1;
    end
end
