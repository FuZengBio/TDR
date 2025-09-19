% [Projection, Summary] = TDR_time(Data, analyzeTimes, params)
% Code adapted by Adam Zaidel (used Mark Churchland's jPCA code as a start) 
%
% OUTPUTS:
%   Projection is a struct with one element per condition. It contains the following fields:
%       .allTimes                   All the times (exactly the same as Data(1).times)
%       .times                      Those times that were used
%       .tradPCAproj                The traditional PCA projection
%       .tradPCAprojAllTimes        Above but for all times.
%       .TargetedPCAproj            The TDR projection
%       .TargetedPCAprojAllTimes    Above but for all times.
%
%   Summary contains the following fields:
%       .PCs                        The PCs (each column is a PC, each row a neuron)
%       .varCaptEachPC              The data variance captured by each PC
%       .varCaptEachPCT             The data variance captured by each PC over time
%       .varCaptEachTargetedPC      The data variance captured by each targeted PC
%       .varCaptEachTargetedPCT     The data variance captured by each targeted PC over time
%       .B_denoised                 the denoised beta parameters from the regression
%       .meanA                      self-explanatory (so that after adaptation I can use the original meanA)
%       .Btheta                     angle between betas
%
% INPUTS:
% The input 'Data' needs to be a struct, with one entry per condition.
% For a given condition, Data(c).A should hold the data (e.g. firing rates).
% Each column of A corresponds to a neuron, and each row to a timepoint.
%
% Data(c).times is an optional field.  If you provide it, only those entries that match
% 'analyzeTimes' will be used for the analysis. If  analyzeTimes == [], all times will be used.
% If you don't provide it, a '.times' field is created that starts at 1.  'analyzeTimes' then refers to those times.
%
% 'params' is optional, and can contain the following fields:
%   .B              The beta parameters
%   .numPCs         Default is 12. The number of traditional PCs to use
%   .normalize      Default is 'true'.  Whether or not to normalize each neurons response by its FR range.
%   .softenNorm     Default is 10.  Determines how much we undernormalize for low FR neurons.  0 means
%                   complete normalization.  10 means a neuron with FR range of 10 gets mapped to a range of 0.5.
%   .meanSubtract   Default is true (1).  Whether or not we remove the across-condition mean from each neurons rate.
%                   if meanSubtract==2 subtract a different meanA (entered also in params)
%   .meanA          (optional, needed if meanSubtract==2)
%   .PCs            (optional) use these input defined PC vectors

function [Projection, Summary] = TDR_time(Data, analyzeTimes, params)


%% quick bookkeeping

numConds = length(Data);
numTimes = size(Data(1).A,1);
% there is also a 'numAnalyzedTimes' defined below.

%% setting parameters that may or may not have been specified

% 'times' field
% if the user didn't specify a times field, create one that starts with '1'
if ~isfield(Data(1),'times')
    for c = 1:length(Data)
        Data(c).times = 1:numTimes;
    end
end

if exist('analyzeTimes', 'var') && ~isempty(analyzeTimes) && max(diff(analyzeTimes)) > max(diff(Data(1).times))
    disp('error, you can use a subset of times but you may not skip times within that subset');
    Projection = []; Summary = []; return;
end

% the number of PCs to look within
numPCs = 12;
if exist('params', 'var') && isfield(params,'numPCs')
    numPCs = params.numPCs;
end
if rem(numPCs,2)>0
    disp('you MUST ask for an even number of PCs.'); return;
end

% do we normalize
normalize = true;
if exist('params', 'var') && isfield(params,'normalize')
    normalize = params.normalize;
end

% do we soften the normalization (so weak signals stay smallish)
% numbers larger than zero mean soften the norm.
% The default (10) means that 10 spikes a second gets mapped to 0.5, infinity to 1, and zero to zero.
% Beware if you are using data that isn't in terms of spikes/s, as 10 may be a terrible default
softenNorm = 10;
if exist('params', 'var') && isfield(params,'softenNorm')
    softenNorm = params.softenNorm;
end

% do we mean subtract
meanSubtract = 1; % 1=subtract off the across-condition mean from each neurons response; 2=subtract a different meanA (entered in params)
if exist('params', 'var') && isfield(params,'meanSubtract')
    meanSubtract = params.meanSubtract;
end
if length(Data)==1, meanSubtract = 0; end  % cant mean subtract if there is only one condition;

if ~exist('analyzeTimes', 'var') || isempty(analyzeTimes)
    disp('analyzing all times');
    analyzeTimes = Data(1).times;
end

if exist('params', 'var') && isfield(params,'PCs') %AZ
    disp('Using predefined PCs')
    predefined_PCs=1;
else
    predefined_PCs=0;
end

%% figure out which times to analyze and make masks

analyzeIndices = ismember(Data(1).times, analyzeTimes);
if size(analyzeIndices,1) == 1
    analyzeIndices = analyzeIndices';  % orientation matters for the repmat below
end
analyzeMask = repmat(analyzeIndices,numConds,1);  % used to mask bigA
if diff( Data(1).times(analyzeIndices) ) <= 5
    disp('mild warning!!!!: you are using a short time base which might make the computation of the derivative a bit less reliable');
end

% these are used to take the derivative
bunchOtruth = true(sum(analyzeIndices)-1,1);
maskT1 = repmat( [bunchOtruth;false],numConds,1);  % skip the last time for each condition
maskT2 = repmat( [false;bunchOtruth],numConds,1);  % skip the first time for each condition

if sum(analyzeIndices) < 5
    disp('warning, analyzing few or no times');
    disp('if this wasnt your intent, check to be sure that you are asking for times that really exist');
end

%% make a version of A that has all the data from all the conditions.
% in doing so, mean subtract and normalize

bigA = vertcat(Data.A);  % append conditions vertically

% note that normalization is done based on ALL the supplied data, not just what will be analyzed
if normalize  % normalize (incompletely unless asked otherwise)
    ranges = range(bigA);  % For each neuron, the firing rate range across all conditions and times.
    normFactors = (ranges+softenNorm);
    bigA = bsxfun(@times, bigA, 1./normFactors);  % normalize
else
    normFactors = ones(1,size(bigA,2));
end

sumA = 0;
for c = 1:numConds
    sumA = sumA + bsxfun(@times, Data(c).A, 1./normFactors);  % using the same normalization as above
end
meanA = sumA/numConds;
if meanSubtract==1  % subtract off the across-condition mean from each neurons response (AZ: over time)
    bigA = bigA-repmat(meanA,numConds,1);
elseif meanSubtract==2 %subtract a different meanA (entered in params)
    bigA = bigA-repmat(params.meanA,numConds,1);
end

%% now do traditional PCA

smallA = bigA(analyzeMask,:);
if ~predefined_PCs %AZ
    [PCvectors,rawScores,latent] = pca(smallA, 'NumComponents', numPCs);  % apply PCA to the analyzed times
    PCvar=cumsum(latent)./sum(latent);
    disp(sprintf('The first %u PCs explain %.0f%% of the variance',2,100*PCvar(2)));
    disp(sprintf('The first %u PCs explain %.0f%% of the variance',3,100*PCvar(3)));
    disp(sprintf('The first %u PCs explain %.0f%% of the variance',6,100*PCvar(6)));
    disp(sprintf('The first %u PCs explain %.0f%% of the variance',numPCs,100*PCvar(numPCs)));
    %plot(PCvar);
else %use the input defined PC vectors
    PCvectors=params.PCs;
    rawScores=bsxfun(@minus, smallA, mean(smallA)) * PCvectors;
end

meanFReachNeuron = mean(smallA);  % this will be kept for use by future attempts to project onto the PCs

% these are the directions in the high-D space (the PCs themselves)
if numPCs > size(PCvectors,2)
    disp('You asked for more PCs than there are dimensions of data');
    disp('Giving up');
end
PCvectors = PCvectors(:,1:numPCs);  % cut down to the right number of PCs

% CRITICAL STEP
% This is what we are really after: the projection of the data onto the PCs
Ared = rawScores(:,1:numPCs);  % cut down to the right number of PCs

% Some extra steps

% projection of all the data
% princomp subtracts off means automatically so we need to do that too when projecting all the data
bigAred = bsxfun(@minus, bigA, mean(smallA)) * PCvectors(:,1:numPCs); % projection of all the data (not just the analyzed chunk) into the low-D space.
% need to subtract off the mean for smallA, as this was done automatically when computing the PC scores from smallA, and we want that projection to match this one.

% projection of the mean
meanAred = bsxfun(@minus, meanA, mean(smallA)) * PCvectors(:,1:numPCs);  % projection of the across-cond mean (which we subtracted out) into the low-D space.

% will need this later for some indexing
numAnalyzedTimes = size(Ared,1)/numConds;

% for targeted dimensionality reduction
D=zeros(size(bigA,2),size(bigA,2));
for i=1:numPCs
    D = D + PCvectors(:,i) * PCvectors(:,i)';
end

B=[]; BK=[];
for i=1:length(params.B)
    % denoised regression coefficients
    B_denoised{i} = D * params.B{i};
    if size(params.B{i},2)>1 %collape over time (not Kernel method)
        [~,I(i)] = max(sqrt(sum(B_denoised{i}.^2,1)));
        Bmax{i} = B_denoised{i}(:,I(i));
        B=[B Bmax{i}];
    else %single B values (instead of over time) - kernel method
        B=[B params.B{i}];
        %B=[B B_denoised{i}]; %AZ seems to work just as well without denoising 
    end
end
pairs=nchoosek(1:length(params.B),2);
for i=1:size(pairs,1)
    costheta = dot(B(:,pairs(i,1)),B(:,pairs(i,2)))/(norm(B(:,pairs(i,1)))*norm(B(:,pairs(i,2))));
    theta(pairs(i,1),pairs(i,2)) = acos(costheta)*180/pi;
    theta(pairs(i,2),pairs(i,1)) = acos(costheta)*180/pi; %just for symmetry    
end

%orthogonalize
[Q,R] = qr(B);
cor = repmat(sign(sum(Q(:,1:length(params.B)) .* B)),length(Q),1); %orthogonalizing can flip the sign! correct for that...
Targeted_PCvectors=Q(:,1:length(params.B)) .* cor; %task related axes (orthogonalized)
%calculate targeted scores
%Targeted_rawScores=bsxfun(@minus, smallA, mean(smallA)) * Targeted_PCvectors;
Targeted_rawScores=smallA * Targeted_PCvectors; %since I "mean subtract" the whole timecourse with either the current (or the previous block's) mean don't subtract the average mean here (NB same needs to be done below for TargetedPCA_AllTimes)

% %calculate targeted scores (not orthogonolized)
% Targeted_PCvectors=B; %task related axes
% Targeted_rawScores=bsxfun(@minus, smallA, mean(smallA)) * Targeted_PCvectors;

%variance captured
origVar = sum(sum( bsxfun(@minus, smallA, mean(smallA)).^2));
smallAcond=reshape(smallA,numAnalyzedTimes,numConds,[]); %smallA sorted by condition
origVarT = sum(sum( bsxfun(@minus, smallAcond, mean(smallAcond,2)).^2,2),3); %total variance computed seperately at each time

%PC
varCaptEachPC = sum(Ared.^2) / origVar;  % this equals latent(1:numPCs) / sum(latent)
%over time
AredCond=reshape(Ared,numAnalyzedTimes,numConds,[]); %Ared sorted by condition
varCaptEachPCT = bsxfun(@rdivide,sum(AredCond.^2,2),origVarT); %over time
%planes
varCaptEachPlanePC = reshape(varCaptEachPC, 2, numPCs/2);
varCaptEachPlanePC = sum(varCaptEachPlanePC);
%Targeted PC
varCaptEachTargetedPC = sum(Targeted_rawScores.^2) / origVar;
%over time
Targeted_rawScoresCond=reshape(Targeted_rawScores,numAnalyzedTimes,numConds,[]); %Targeted_rawScores sorted by condition
varCaptEachTargetedPCT = bsxfun(@rdivide,sum(Targeted_rawScoresCond.^2,2),origVarT);
%varCaptEachTargetedPCT = bsxfun(@rdivide,sum(Targeted_rawScoresCond.^2,2),1);

tradPCA_AllTimes = bsxfun(@minus, bigA, mean(smallA)) * PCvectors;  % mean center in exactly the same way as for the shorter time period.
%TargetedPCA_AllTimes = bsxfun(@minus, bigA, mean(smallA)) * Targeted_PCvectors;  % mean center in exactly the same way as for the shorter time period (see above).
TargetedPCA_AllTimes = bigA * Targeted_PCvectors;  % mean center in exactly the same way as for the shorter time period (see above).

% Do some annoying output formatting.
% Put things back so we have one entry per condition
index1 = 1;
index2 = 1;
for c = 1:numConds
    index1b = index1 + numAnalyzedTimes -1;  % we will go from index1 to this point
    index2b = index2 + numTimes -1;  % we will go from index2 to this point
    
    Projection(c).allTimes = Data(1).times;
    Projection(c).times = Data(1).times(analyzeIndices);
    
    %PC
    Projection(c).tradPCAproj = Ared(index1:index1b,:);
    Projection(c).tradPCAprojAllTimes = tradPCA_AllTimes(index2:index2b,:);
    
    %targeted PCA
    Projection(c).TargetedPCAproj = Targeted_rawScores(index1:index1b,:);
    Projection(c).TargetedPCAprojAllTimes = TargetedPCA_AllTimes(index2:index2b,:);
    
    index1 = index1+numAnalyzedTimes;
    index2 = index2+numTimes;
end

%% summary output structure
Summary.B_denoised=B_denoised; %the denoised beta parameters from the regression
Summary.meanA=meanA; % meanA (so that after adaptation I can use the original meanA)
Summary.PCs = PCvectors;
Summary.Btheta=theta; % angle between betas
Summary.varCaptEachPC = varCaptEachPC;
Summary.varCaptEachPCT = squeeze(varCaptEachPCT);
Summary.varCaptEachPlanePC = varCaptEachPlanePC;
Summary.varCaptEachTargetedPC = varCaptEachTargetedPC;
Summary.varCaptEachTargetedPCT = squeeze(varCaptEachTargetedPCT);
Summary.acrossCondMeanRemoved = meanSubtract;
Summary.preprocessing.normFactors = normFactors;  % Used for projecting new data from the same neurons into the space
Summary.preprocessing.meanFReachNeuron = meanFReachNeuron; % You should first normalize and then mean subtract using this (the original) mean
% conversely, to come back out, you must add the mean back on and then MULTIPLY by the normFactors
% to undo the normalization.

