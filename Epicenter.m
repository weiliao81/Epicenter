clear all;clc

load dataPATWscore.mat;%load w-score matrix,nPat*mReg.
load LabelList.mat;%load cluster assignment 1, 2 mean biotype 1 and 2
load geneexpression.mat; %load AHBA gene expression,  mReg*gGene 

mReg = 308; % number of regions
thr = 0; % threshold for CGE, 0, 0.1, 0.2, 0.3

% split data into two biotype
PATWscore1 = mean(dataPATWscore(find(LabelList == 1),:),1)';
PATWscore2 = mean(dataPATWscore(find(LabelList == 2),:),1)';

% construct CGE matrix
CGE = corr(geneexpression');
CGE(1:size(CGE,1)+1:end) = 0;   % clear diagonal
CGE(CGE < thr) = 0;   % apply for CGE threshold

% prediction of regional w-score using its neighbours'w-score 
for i = 1:mReg
    nodepre1(i,:) = (PATWscore1'*CGE(:,i))/sum(CGE(:,i));
    nodepre2(i,:) = (PATWscore2'*CGE(:,i))/sum(CGE(:,i));
end
[rpre1,~] = corr(nodepre1,PATWscore1);
[rpre2,~] = corr(nodepre2,PATWscore2);
figure; plot(PATWscore1,nodepre1,'o');
figure; plot(PATWscore2,nodepre2,'o');

% sort PATWscore and predicted PATWscore (PATWscore_neighbor1) 
PATWscore_neighbor1 = nodepre1; 
PATWscore_neighbor2 = nodepre2;

[B00,I00] = sort(abs(PATWscore1)); 
[B01,I01] = sort(abs(PATWscore_neighbor1));
  
[B10,I10] = sort(abs(PATWscore2));
[B11,I11] = sort(abs(PATWscore_neighbor2));
    
for j = 1:mReg
    noderank00(j,:) = find(I00 == j);
    noderank01(j,:) = find(I01 == j);
    noderank1(j,:) = (noderank00(j,:) + noderank01(j,:))/2;
    
    noderank10(j,:) = find(I10 == j);
    noderank11(j,:) = find(I11 == j);
    noderank2(j,:) = (noderank10(j,:) + noderank11(j,:))/2;
end

