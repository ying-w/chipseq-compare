% function to calculate the difference between segment counts in each
% chromosomes also calculate the mean of counts in each chromosomes
%
% NOTE: Chromosome MT is skipped in this code!
%
% input:
% dataPath   - path to sample data
% refName    - name of reference sample
% sampleName - name of sample you want to compare with the reference
% binSize    - size of each bin
% fSeqDepth  - is a flag if 1 do sequence depth normalization
%
% Each cell in sample data is a count of fragments which falls in a bin of
% size <binSize> in each chromosome, which are saved in a vector called h
% in <dataPath><refName>\<refName>_<binSize>_chr<i>, where i is the
% chromosome index.
%
% For example: HS0268_1000_chr1 is a (1 x n) vector, where n is the number
% of bins of size 1000. HS0268_1000_chr1(1) is the number of fragments in
% genomic location 1 to 1000 bp in chr1. HS0268_1000_chr1(2) is the number
% of fragments in genomic location 1001 to 2000 bp in chromosome 1.
%
% output:
% <refName>_<sampleName>_MAdata.mat save in <matfile>
%
% copyright 2009 Cenny Taslim
% Please cite: Taslim C, Wu J, Yan P, Singer G, Parvin J, Huang T, Lin S,
% Huang K. Comparative study on ChIP-seq data: normalization and binding
% pattern characterization, Bioinformatics, Vol. 25, No. 18. (15 September
% 2009), pp. 2334-2340

function [matfile] = normalizeMeanVarStep1(dataPath,refName,sampleName,binSize,chrSet,fSeqDepth)

fpath0 = strcat(dataPath,refName,filesep,refName,'_',num2str(binSize),'_chr');
fpath1 = strcat(dataPath,sampleName,filesep,sampleName,'_',num2str(binSize),'_chr');

if ~exist('chrSet','var')
    chrSet = [1:22 24:25];
end;

if ~exist('fSeqDepth','var')
    fSeqDepth = 0;
end;

% initialize chr_size to store "real" chromosome size or length
chr_size=zeros(25,1);

for i = chrSet
    
    % display progress
    disp(strcat({'processing Step 1 normalization chromosome '},num2Chr(i)));
    % get data for each chromosome
    filepath0=strcat(fpath0,num2Chr(i),'.mat');
    filepath1=strcat(fpath1,num2Chr(i),'.mat');
    filepath0=strcat(fpath0,num2Chr(i));
    filepath1=strcat(fpath1,num2Chr(i));
    
    #h0=load(filepath0);
    #data0=h0.h';
    #h1=load(filepath1);
    #data1=h1.h';
    data0 = load(filepath0);
    data1 = load(filepath1);
           
    % make the number of bins of both data equal
    % ok! just after effect of previous step which stops counting when
    % there is no counts in last few bins
    if size(data0,1)< size(data1,1)
        data0(size(data0,1)+1:size(data1,1),1)=0;
    elseif size(data0,1) > size(data1,1)
        data1(size(data1,1)+1:size(data0,1),1)=0;
    end;
    % store chromosome length
    chr_size(i) = length(data0);
    
    % save original data (- centromere)
    unScaleData1 = data1;
    if fSeqDepth
        % make the total counts the same for both data
        sum0 = sum(data0);
        sum1 = sum(data1);
        scale = sum0/sum1;
        unScaleData1 = data1;
        data1 = data1.* scale;
    end;
    % make the size of data equal in all chromosome
    if (i ~= chrSet(1)&& size(data0,1)>size(meanCt,1))
        unScaleMeanCt(size(meanCt,1)+1:size(data0,1),1:size(meanCt,2))=0;
        unScaleDiffCt(size(diffCt,1)+1:size(data0,1),1:size(meanCt,2))=0;
        ref_smpl(size(meanCt,1)+1:size(data0,1),1:size(meanCt,2))=0;
        treat_smpl(size(meanCt,1)+1:size(data0,1),1:size(meanCt,2))=0;
        unScaleTreat_smpl(size(meanCt,1)+1:size(data0,1),1:size(meanCt,2))=0;
        meanCt(size(meanCt,1)+1:size(data0,1),1:size(meanCt,2))=0;
        diffCt(size(diffCt,1)+1:size(data0,1),1:size(meanCt,2))=0;
    elseif (i ~= chrSet(1) && size(data0,1)<size(meanCt,1))
        data0(size(data0,1)+1:size(meanCt,1),1)=0;
        data1(size(data1,1)+1:size(meanCt,1),1)=0;
        unScaleData1(size(unScaleData1,1)+1:size(meanCt,1),1)=0;
    end;
    % store the counts of original data (reference and treated sample) the
    % t stands for transposed data
    ref_smpl(:,i)=data0;
    treat_smpl(:,i)=data1;
    unScaleTreat_smpl(:,i)=unScaleData1;
    % store counts of mean and difference between reference and treated
    meanCt(:,i)=(data1+data0)./2;
    diffCt(:,i)= data1-data0;
    unScaleMeanCt(:,i) = (unScaleData1+data0)./2;
    unScaleDiffCt(:,i) = (unScaleData1-data0);
end;
tmpPath0 = strsplit(filepath0, filesep);
outPath = regexprep(filepath0,tmpPath0{end},'');

% create output folder if not exist
if(exist(outPath,'dir')==0)
    mkdir(outPath);
end
matfile = strcat(outPath,refName,'_',sampleName,'_MAdata.mat');
save(matfile, 'meanCt','diffCt', 'ref_smpl', 'treat_smpl','chr_size',...
    'unScaleMeanCt', 'unScaleDiffCt','unScaleTreat_smpl');
% dlmwrite(strcat(outPath,refName,'_',sampleName,'_diffCt.txt'), diffCt) 
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

% test1 is done on STAT1_HeLa vs input chr2
% checked: unScaleDiffCt,unScaleMeanCt, diffCt, meanCt, ref_smpl,
% treat_smpl

