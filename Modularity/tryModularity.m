% run PCfdr on the data

clear all,
clc,

load('STC_CTC_N006_ROIs.mat')
dim = ROItcs{1}.hdr.Dimensions;

X1 = double([ROItcs{8}.DCDdata, ROItcs{12}.DCDdata]);
Xpointer1 = [ROItcs{8}.pointer; ROItcs{12}.pointer];


% X2 = double([ROItcs{21}.DCDdata, ROItcs{25}.DCDdata]);
% Xpointer2 = [ROItcs{21}.pointer; ROItcs{25}.pointer];

% idx = 1;
% for i = [1:7,9:11,13:52];
%     R(:,idx) = mean(double(ROItcs{i}.DCDdata),2);
%     idx = idx + 1;
% end

xN1 = length(Xpointer1);
for i = 1:xN1
    vol = zeros(dim);
    vol(Xpointer1(i)) = 1;
    [Xpos1(i) Ypos1(i) Zpos1(i)] = find(vol);
end

for i = 1:xN1
    Dis1(:,i) = sqrt((Xpos1(i) - Xpos1).^2 + (Ypos1(i) - Ypos1).^2 + (Zpos1(i) - Zpos1).^2);
    Dis1(i,i) = 0;
end

% xN2 = length(Xpointer2);
% for i = 1:xN2
%     vol = zeros(dim);
%     vol(Xpointer2(i)) = 1;
%     [Xpos2(i) Ypos2(i) Zpos2(i)] = find(vol);
% end
% 
% for i = 1:xN2
%     Dis2(:,i) = sqrt((Xpos2(i) - Xpos2).^2 + (Ypos2(i) - Ypos2).^2 + (Zpos2(i) - Zpos2).^2);
%     Dis2(i,i) = 0;
% end

% load the Task index for segmentation purpose  TaskIndex
load('STC_CTC_N006_timeindex.mat');

XX1 = X1(TaskIndex==1,:);
% XX2 = X2(TaskIndex==1,:);
% XR = R(TaskIndex==1,:);

%band = max(Xpos) - min(Xpos) + max(Ypos) - min(Ypos) + max(Zpos) - min(Zpos);
CC1 = corrcoef(XX1);
CC1 = abs(CC1);
CC1(CC1<0.4) = 0;
for i = 1:xN1
    CC1(i,i) = 0;
end

% CC2 = corrcoef(XX2);
% CC2 = abs(CC2);
% CC2(CC2<0.4) = 0;
% for i = 1:xN2
%     CC2(i,i) = 0;
% end


% S = cov(XX);
% 
% % S -- Sample covariance matrix
% % r -- number of Neumann series terms
% % e -- level of replacing small / negative eigenvalues to ensure psd
% % cutoff -- trunction of precision matrix to promote sparsity
% r = 5;
% e = 0.01;
% cutoff = 0.001;
% Omega = neumann(S, r, e, cutoff);

[Ci1 Q1]=spa_modularity_und(CC1, Dis1, 500, 2);
% [Ci2 Q2]=spa_modularity_und(CC2, Dis2, 500, 2);

% DT1 = Dis1;
% DT1(DT1>100) = 0;
% DT2 = Dis2;
% DT2(DT2>100) = 0;
% 
% [Ci1 Q1]=modularity_und_2(DT1, 2);
% [Ci2 Q2]=modularity_und_2(DT2, 2);

% [Ci Q]=modularity_und_2(CC,2);
% 
% [Ci Q]=spa_modularity_und(CC, Dis, band, 2);
% 
% Dtmp = Dis;
% Dtmp = 1./Dtmp;
% Dtmp(Dtmp == inf) = 0;
% Dtmp(Dtmp<0.01) = 0;
% [Ci Q]=modularity_und_2(Dtmp,2);
% 
% [Ci Q] = dis_modularity_und(CC, Dis, 500, 2); 
% 
% TTT = Dis;
% TTT(TTT>500) = 0;
% TTT = CC.*logical(TTT);
% [Ci Q]=modularity_und_2(TTT,2);
clf
figure; plot(Xpos1(Ci1==1),Ypos1(Ci1==1),'r+'); hold on
scatter(Xpos1(Ci1==2),Ypos1(Ci1==2),'g'); hold on
xlabel('X position')
ylabel('Y position')
legend('Dorsal-Striatum','Ventral-Striatum')


% plot(Xpos2(Ci2==1),Ypos2(Ci2==1),'b*'); hold on
% plot(Xpos2(Ci2==2),Ypos2(Ci2==2),'m+'); hold on


