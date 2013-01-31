% Function to do mean and variance normalization.
% !!!! ASSUMPTION: difference between the two samples are affecting a small
% proportion of the genome. Thus if the reference is input DNA, use
% normalizeMeanVarStep2_1.m instead!!
%
% NOTE: chr MT is skipped.
%
% input: 
% dataPath - path to matlab matrix contains un-normalized meanCt,
%            diffCt,ref_smpl and treat_smpl run normalizeMeanVarStep1.m to 
%            get <refName>_<sampleName>_MAdata.mat 
% span     - 2 parameters for lowess regression (1st - mean normalization,
%            2nd - variance normalization), default(0.6,0.1). Must be
%            between 0 and 1. Smaller value gives smoother curve.
% output: 
% .mat file containing normalized counts saved in <matfile>
%
% copyright 2009 Cenny Taslim
% Please cite: Taslim C, Wu J, Yan P, Singer G, Parvin J, Huang T, Lin S, 
% Huang K. Comparative study on ChIP-seq data: normalization and binding 
% pattern characterization, Bioinformatics, Vol. 25, No. 18. (15 September
% 2009), pp. 2334-2340

function [matfile] = normalizeMeanVarStep2(dataPath,span,chrSet)
load(dataPath);
if ~exist('span','var') || length(span) < 2 || span(1) > 1 || span(1) < 0 ...
        || span(2) > 1 || span(2) < 0
    span =[0.6 0.1];
end;
spanMean = span(1);
spanVar = span(2);

if ~exist('chrSet','var')
    % chromosomes to be processed
    % include 23 below to add chromosome MT
    chrSet = [1:22 24:25];
end;

%initialize matrices
diffCt_loessMean=zeros(size(diffCt,1),25);
norm_diffCt=zeros(size(diffCt,1),25);
norm_diffCt_loessMean=zeros(size(diffCt,1),25);
norm_diffCt_loessVar=zeros(size(diffCt,1),25);
varNorm_diffCt_loessVar=zeros(size(diffCt,1),25);
varNorm_diffCt=zeros(size(diffCt,1),25);
max_chr_size=max(chr_size);
varNorm_y=zeros(max_chr_size,25);

for i=chrSet
    % display progress
    disp(strcat({'processing Step 2 normalization chromosome '},num2Chr(i)));
    
    % estimate mean using loess smoother and spanMean of data
    diffCt_loessMean(1:chr_size(i),i)=smooth(meanCt(1:chr_size(i),i),diffCt(1:chr_size(i),i),spanMean,'loess');
    % substract estimated mean from observations points
    norm_diffCt(1:chr_size(i),i)=diffCt(1:chr_size(i),i)-diffCt_loessMean(1:chr_size(i),i);
    % calculate a loess smoother for plotting purpose
    norm_diffCt_loessMean(1:chr_size(i),i)=smooth(meanCt(1:chr_size(i),i),norm_diffCt(1:chr_size(i),i),spanMean,'loess');
    % estimate a running variance using loess smoother and 10% of data
    norm_diffCt_loessVar(1:chr_size(i),i)=smooth(meanCt(1:chr_size(i),i),abs(norm_diffCt(1:chr_size(i),i)),spanVar,'loess');
    % standardize observation points by dividing it with its estimated variance
    % pre-processing avoid NaN (0 divided by 0). Changing all entries of 0
    % to eps
    norm_diffCt(norm_diffCt==0)=eps;
    norm_diffCt_loessVar(norm_diffCt_loessVar==0)=eps;
    varNorm_diffCt(1:chr_size(i),i)=norm_diffCt(1:chr_size(i),i)./norm_diffCt_loessVar(1:chr_size(i),i);
    % loess smoother for plotting purpose
    varNorm_diffCt_loessVar(1:chr_size(i),i)=smooth(meanCt(1:chr_size(i),i),varNorm_diffCt(1:chr_size(i),i),spanMean,'loess');
    % calculate the normalized data in original units
    varNorm_y(1:chr_size(i),i)=varNorm_diffCt(1:chr_size(i),i)+ref_smpl(1:chr_size(i),i);
    % save normalization for each individual chromosome
    tmpPath1 = strsplit(dataPath, filesep);
    outPath = regexprep(dataPath,tmpPath1{end},'');
    filename = regexprep(tmpPath1{end},'_MAdata.mat','_NormzMeanVar_chr');
    matfile = strcat(outPath,filename,num2Chr(i),'.mat');
    save(matfile,'diffCt_loessMean','norm_diffCt','norm_diffCt_loessMean',...
        'norm_diffCt_loessVar','varNorm_diffCt_loessVar','varNorm_diffCt',...
        'varNorm_y','chr_size','meanCt','diffCt', 'ref_smpl', 'treat_smpl',...
        'unScaleMeanCt', 'unScaleDiffCt','unScaleTreat_smpl');
end;
tmpPath1 = strsplit(dataPath, filesep);
outPath = regexprep(dataPath,tmpPath1{end},'');
filename = regexprep(tmpPath1{end},'_MAdata.mat','_NormzMeanVar.mat');
matfile = strcat(outPath,filename);
% create folder if not exist
if(exist(outPath,'dir')==0)
    mkdir(outPath);
end
save(matfile,'diffCt_loessMean','norm_diffCt','norm_diffCt_loessMean',...
    'norm_diffCt_loessVar','varNorm_diffCt_loessVar','varNorm_diffCt',...
    'varNorm_y','chr_size','meanCt','diffCt', 'ref_smpl', 'treat_smpl',...
    'unScaleMeanCt', 'unScaleDiffCt','unScaleTreat_smpl');
end

function strChr = num2Chr(i)
if i <= 22
    strChr = num2str(i);
elseif i ==23
    strChr = 'MT';
elseif i ==24
    strChr = 'X';
else
    strChr = 'Y';
end;
end
