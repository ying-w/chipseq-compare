% Script to get difference counts first then normalize data
% see normalizeMeanVarStep1.m and normalizeMeanVarStep2.m for more details
%
% copyright 2009 Cenny Taslim
% Please cite: Taslim C, Wu J, Yan P, Singer G, Parvin J, Huang T, Lin S, 
% Huang K. Comparative study on ChIP-seq data: normalization and binding 
% pattern characterization, Bioinformatics, Vol. 25, No. 18. (15 September
% 2009), pp. 2334-2340

source("normalizeMeanVarStep1.m");
source("normalizeMeanVarStep2.m");

clear all
clc

% folder path to both samples.
#dataPath = strcat(cd,filesep,'data',filesep);
dataPath = strcat('./data/');
% name of the reference sample. This is also the filename, i.e.
% <dataPath><refName>\<refName>_<binSize>_chr<i>, will be the full path to
% reference sample. This is also the filename
refName = 'none1';
% name of the sample to be compared against reference. 
sampleName = 'high1';

% size of each bin
binSize = 500;
% parameters for loess regression, 1st number is for normalization with
% respect to mean, 2nd number is for normalization w.r.t. variance
span = [0.6 0.1];
% flag whether to do sequence depth normalization
fSeqDepth = 0;
% specify chromosome to process. 23 = Chromosome MT, 24 = chr X, 25 = chr Y
% to process all chromosome except chr MT "chrSet = [1:22 24:25];"
chrSet = [1:22 24:25]; #change to 25 for Y only

path1 = normalizeMeanVarStep1(dataPath,refName,sampleName,binSize,...
    chrSet,fSeqDepth);
disp('step 2.... RUN IN R');
%disp('step 2....');
%path2 = normalizeMeanVarStep2(path1,span,chrSet);

%-------------------------------------------------------------------------%
% plot the normalized data
% load(path2);
% i = 25
% plot(meanCt(1:chr_size(i),i),varNorm_diffCt(1:chr_size(i),i),'.');
% xlabel('Mean Binding Quantity');
% ylabel('mean-variance Normalized Difference Binding Quantity');
% title('Mean and Variance Normalized');
%-------------------------------------------------------------------------%


