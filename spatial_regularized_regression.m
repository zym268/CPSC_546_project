% 
clc,
clear;

load('STC_CTC_N006_ROIs.mat')
dim = ROItcs{1}.hdr.Dimensions;

X = double([ROItcs{8}.DCDdata, ROItcs{12}.DCDdata]);
Xpointer = [ROItcs{8}.pointer; ROItcs{12}.pointer];

idx = 1;
for i = [1:7,9:11,13:17,19:52];
    R(:,idx) = mean(double(ROItcs{i}.DCDdata),2);
    idx = idx + 1;
end

Y = mean(double(ROItcs{18}.DCDdata),2);

% load the Task index for segmentation purpose  TaskIndex
load('STC_CTC_N006_timeindex.mat');

YY = Y(TaskIndex==1,:);
XX = X(TaskIndex==1,:);
XR = R(TaskIndex==1,:);

xN = length(Xpointer);


lambda = 1;
gamma = 30;
XR=[];
xN = size(XX,2);

% initial the distance matrix
D = zeros(xN,xN);

% get the position for each point in X
for i = 1:xN
    vol = zeros(dim);
    vol(Xpointer(i)) = 1;
    [Xpos(i) Ypos(i) Zpos(i)] = find(vol);
end

% generate the distance matrix
for i = 1:xN
%   D(:,i) = sqrt((Xpos(i) - Xpos).^2 + (Ypos(i) - Ypos).^2 + (Zpos(i) - Zpos).^2);
    D(:,i) = sqrt((Xpos(i) - Xpos).^2 + (Ypos(i) - Ypos).^2);
    D(i,i) = inf;
end

% thresholding the distance matrix
minD = min(D);
Con = zeros(xN, xN);
for i = 1:xN
    tmp = D(:,i);
%     tmp(tmp > minD(i)) = 0;
    tmp(tmp > 1.5) = 0;
    Con(:,i) = tmp;
    Con(i,i) = 1;
end
Con = Con + Con';

% first check the spatial continous of the input XX
[H, NewInd] = CheckSpatialCont(Con, 1:xN);
if ~H
    warning('The input voxels are not spatial contiunous! They are assumed to be continuous in the following analysis');
    Con = MakeSpatialCont(D, Con, NewInd);
end

% initialize the number of group
K = xN;
tic;
beta = SpatialReguRegOnX_new( [XX,XR], YY, D, lambda, gamma);
%beta = SpatialReguRegOnX_SPG( [XX,XR], YY, D, lambda, gamma);
%[beta, iteration,history] = SpatialReguRegOnX_ADMM( [XX,XR], YY, D, lambda, gamma);
toc;
normD=0;
for i = 1:length(D)
   for j=i+1:length(D)
      normD=normD+D(i,j);
   end
end

alpha=1;
k=1;
for i = 1:length(beta)
   for j=i+1:length(beta)
       Y(k)= alpha*(beta(i)-beta(j))/norm(beta)*D(i,j)/normD+(1-alpha)*D(i,j)/normD;
       %Y(k)= (beta(i)-beta(j))/norm(beta)*D(i,j);
       %Y(k)= D(i,j);
       k=k+1;
   end
end
Y=Y';

Z = linkage(Y,'ward');
T = cluster(Z,'maxclust',2);
group1=find(T==1);
group2=find(T==2);

% get the position for each point in X
for i = 1:length(group1)
    vol = zeros(dim);
    vol(Xpointer(group1(i))) = 1;
    [Xpos_group1(i) Ypos_group1(i) Zpos(i)] = find(vol);
end
% get the position for each point in X
for i = 1:length(group2)
    vol = zeros(dim);
    vol(Xpointer(group2(i))) = 1;
    [Xpos_group2(i) Ypos_group2(i) Zpos(i)] = find(vol);
end
plot(Xpos_group1, Ypos_group1,'r+');
hold on
scatter(Xpos_group2, Ypos_group2,'g');
xlabel('X position')
ylabel('Y position')
legend('Dorsal-Striatum','Ventral-Striatum')