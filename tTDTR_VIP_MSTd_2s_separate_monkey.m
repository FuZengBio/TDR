function VAC_monkey_2s_batch()
clear; clc;
%% =========================================================
% unified paths
%% =========================================================
FilePath   = 'Z:\Users\TDR\Data\';

%% =========================================================
% parameters
%% =========================================================
sigma = 30;
half_guas = sigma*3;
dist = (-(half_guas):(half_guas))';
window = 1/(sigma*sqrt(2*pi))*exp(-((dist.^2)/(2*sigma^2)))*1000;

span=10;
t1=101:2000;
REG_t=square_sum(t1,span)/span;
T = numel(REG_t);

% VAC kernel regression params
REG_params_kernel={'VELangle','ACCangle','choice'};
REG_vel_dimK=find(strcmp(REG_params_kernel,'VELangle'));
REG_acc_dimK=find(strcmp(REG_params_kernel,'ACCangle'));
REG_choice_dimK=find(strcmp(REG_params_kernel,'choice'));

% ideal kernel
kernel_del=0;
[gausacc1, gausvel1, gauspos1]=genGaussian(0.13,6,2);
REG_acc_ideal1=(square_sum(gausacc1((t1)-kernel_del)',span)/span);
REG_acc_ideal1=REG_acc_ideal1./sqrt(REG_acc_ideal1 * REG_acc_ideal1');

REG_vel_ideal1=(square_sum(gausvel1((t1)-kernel_del)',span)/span);
REG_vel_ideal1=REG_vel_ideal1./sqrt(REG_vel_ideal1 * REG_vel_ideal1');

REG_pos_ideal1=(square_sum(gauspos1((t1)-kernel_del)',span)/span);
REG_pos_ideal1=REG_pos_ideal1./sqrt(REG_pos_ideal1 * REG_pos_ideal1');

kernel=[];
kernel(REG_choice_dimK,:)=REG_pos_ideal1;
kernel(REG_vel_dimK,:)=REG_vel_ideal1;
kernel(REG_acc_dimK,:)=REG_acc_ideal1;

% axis order for VAC
id=[REG_vel_dimK REG_acc_dimK REG_choice_dimK];
% id=[REG_choice_dimK REG_acc_dimK REG_vel_dimK];
% id=[REG_vel_dimK REG_choice_dimK REG_acc_dimK];

dim2name = containers.Map( ...
    [REG_vel_dimK, REG_acc_dimK, REG_choice_dimK], {'VEL','ACC','CHOICE'} );

axisNamesOrdered = cell(1,numel(id));
for k=1:numel(id)
    axisNamesOrdered{k} = dim2name(id(k));
end

fprintf('VAC axis order: %s\n', strjoin(axisNamesOrdered,' -> '));

% significance params
nPerm    = 1000;
alphaPt  = 0.05;
alphaCl  = 0.05;

stim_type=[1,2]; % 1=Vest, 2=Vis
durName = '2s';

%% =========================================================
% run list (2 datasets, in order)
%% =========================================================
datasets = { ...
    struct('area','VIP',  ...
           'mat','VIP_Data_2s_REG.mat', ...
           'label','monkey_lable_VIP_2s.mat', ...
           'monkeyName',{{'U','C'}}), ...
    struct('area','MSTd', ...
           'mat','MSTd_Data_2s_REG.mat', ...
           'label','monkey_lable_MSTd_2s.mat', ...
           'monkeyName',{{'A2','C'}}) ...
    };

for d = 1:numel(datasets)

    ds = datasets{d};

    fprintf('\n====================================================\n');
    fprintf('Running %d/%d: %s %s\n', d, numel(datasets), ds.area, durName);
    fprintf('====================================================\n');

    matfile   = fullfile(FilePath, ds.mat);
    labelfile = fullfile(FilePath, ds.label);

    if ~exist(matfile, 'file')
        warning('Data file not found: %s. Skip.', matfile);
        continue;
    end
    if ~exist(labelfile, 'file')
        warning('Label file not found: %s. Skip.', labelfile);
        continue;
    end

    run_one_dataset(matfile, labelfile, ds.area, ds.monkeyName, durName, ...
        window, T, REG_t, kernel, ...
        REG_params_kernel, REG_vel_dimK, REG_acc_dimK, REG_choice_dimK, ...
        id, axisNamesOrdered, nPerm, alphaPt, alphaCl, stim_type);

end

fprintf('\nAll datasets finished.\n');

end


%% =========================================================
function run_one_dataset(matfile, labelfile, areaName, monkeyName, durName, ...
    window, T, REG_t, kernel, ...
    REG_params_kernel, REG_vel_dimK, REG_acc_dimK, REG_choice_dimK, ...
    id, axisNamesOrdered, nPerm, alphaPt, alphaCl, stim_type)

load(matfile);
load(labelfile);   % loads monkey_lable

idxM1_all = find(monkey_lable==1);
idxM2_all = find(monkey_lable==2);
if isempty(idxM1_all) || isempty(idxM2_all)
    error('%s: monkey_lable must contain both 1 and 2.', areaName);
end

OutputPath = 'Z:\Users\TDR\results\tTDR_separate_monkey\';
if ~exist(OutputPath, 'dir')
    mkdir(OutputPath);
end

for ss=1:length(stim_type)
    if ss==1
        stimname='Vest';
        VELcell = REG_VELangle(:,1);
        ACCcell = REG_ACCangle(:,1);
        CHOcell = REG_choice(:,1);
    else
        stimname='Vis';
        VELcell = REG_VELangle(:,2);
        ACCcell = REG_ACCangle(:,2);
        CHOcell = REG_choice(:,2);
    end
    fprintf('\n=========== %s | %s %s ===========\n', stimname, areaName, durName);

    %% -----------------------------------------------------
    % Part 1: build trial-wise SDF per neuron
    %% -----------------------------------------------------
    Test_REG_r = cell(size(REG_raster_per_trial,1),1);

    for nCell=1:size(REG_raster_per_trial,1)
        temp_spike=REG_raster_per_trial{nCell,ss};
        sdf=[];

        if isempty(temp_spike)
            Test_REG_r{nCell,1} = [];
            continue;
        end

        for i=1:size(temp_spike,2)
            spike=squeeze(temp_spike(:,i,:));
            [ntrials,~] = size(spike);
            for ii=1:ntrials
                convspike = conv(spike(ii,:), window, 'same');
                convspike = convspike(101:2000);
                tempdata  = downsample(convspike,10);
                sdf((i-1)*ntrials+ii,:) = tempdata;
            end
        end
        Test_REG_r{nCell,1} = sdf;
    end

    %% -----------------------------------------------------
    % Part 2: VAC regression to get REG_b_use
    %% -----------------------------------------------------
    REG_bVesK_Test = cell(numel(REG_params_kernel),1);
    for b=1:numel(REG_params_kernel)
        REG_bVesK_Test{b} = nan(size(Test_REG_r,1),1);
    end

    for Ni=1:size(Test_REG_r,1)

        tempdata = REG_Stim{Ni};
        if isempty(tempdata), continue; end
        TMP_Ves = (tempdata==1);
        if sum(TMP_Ves)<=0, continue; end

        idxTrials = find(TMP_Ves);
        nTr = numel(idxTrials);

        TMP = nan(numel(REG_params_kernel), nTr);
        TMP(REG_vel_dimK,:)    = VELcell{Ni}(idxTrials);
        TMP(REG_acc_dimK,:)    = ACCcell{Ni}(idxTrials);
        TMP(REG_choice_dimK,:) = CHOcell{Ni}(idxTrials);

        TMP_FVesK = [];
        for trl=1:nTr
            TMP_FVesK=[TMP_FVesK repmat(TMP(:,trl),1,sum(REG_t<=2000)) .* kernel]; %#ok<AGROW>
        end

        temp_r = Test_REG_r{Ni};
        if isempty(temp_r), continue; end
        temp_r = temp_r(idxTrials,1:sum(REG_t<=2000));
        r = bsxfun(@minus,temp_r,mean(temp_r));

        if rcond(TMP_FVesK*(TMP_FVesK')) < 1e-12
            continue
        end

        TMP_bK = ((TMP_FVesK*(TMP_FVesK'))^-1)*TMP_FVesK*reshape(r',[],1);
        for b=1:numel(REG_params_kernel)
            REG_bVesK_Test{b}(Ni,1)=TMP_bK(b);
        end
    end

    REG_b_use = REG_bVesK_Test(id);

    %% -----------------------------------------------------
    % Part 3: build condition-mean matrix DataRun
    %% -----------------------------------------------------
    clear Test_PCA_Data
    for nCell=1:size(REG_raster_per_trial,1)
        temp_spike=REG_raster_per_trial{nCell,ss};

        if isempty(temp_spike), continue; end

        for i=1:size(temp_spike,1)
            spike=squeeze(temp_spike(i,:,:));
            ntrials=size(spike,1);
            temp_mean = nan(ntrials,T);

            for ii=1:ntrials
                convspike = conv(spike(ii,:), window, 'same');
                convspike = convspike(101:2000);
                tempdata  = downsample(convspike,10);
                temp_mean(ii,:) = tempdata;
            end

            Test_PCA_Data(1,i).A(:,nCell) = mean(temp_mean,1,'omitnan')';
            Test_PCA_Data(1,i).times      = REG_t;
        end
    end

    DataRun = Test_PCA_Data;
    REG_params_use=REG_params_kernel(id);
    PCA_params = struct();
    PCA_params.B = REG_b_use;
    PCA_params.normalize = false;
    PCA_params.numPCs = 12;
    times = [];

    [~,Summary_all] = TDR_time(DataRun, times, PCA_params);
    Vall = Summary_all.varCaptEachTargetedPCT; %#ok<NASGU>

    %% -----------------------------------------------------
    % Monkey-specific variance on SAME axes + significance
    %% -----------------------------------------------------
    fprintf('[%s] Monkey variance + significance...\n', stimname);

    Nall = size(DataRun(1).A,2);
    idxM1 = find(monkey_lable==1);
    idxM2 = find(monkey_lable==2);

    if isempty(idxM1) || isempty(idxM2)
        warning('One monkey has no cells; skip significance.');
        continue;
    end

    [V_M1, M1proj, tvec] = computeVar_subsetVAC(DataRun, PCA_params, idxM1);
    [V_M2, M2proj, ~]    = computeVar_subsetVAC(DataRun, PCA_params, idxM2);

    D_real = V_M1 - V_M2;

    nAxes = numel(id);
    Tlen  = size(V_M1,1);
    nM1 = numel(idxM1);
    allIdx = (1:Nall)';

    D_null = nan(Tlen, nPerm, nAxes);
    for p=1:nPerm
        rp = randperm(Nall);
        idxM1p = allIdx(rp(1:nM1));
        idxM2p = allIdx(rp(nM1+1:end));

        V_M1p = computeVar_subsetVAC(DataRun, PCA_params, idxM1p);
        V_M2p = computeVar_subsetVAC(DataRun, PCA_params, idxM2p);

        D_null(:,p,:) = V_M1p - V_M2p;
    end

    sig = struct();
    for k=1:nAxes
        axName = axisNamesOrdered{k};
        [mask, thrPoint, clThr] = clusterPerm_twoSided(D_real(:,k), squeeze(D_null(:,:,k)), alphaPt, alphaCl);
        sig.(axName).mask = mask;
        sig.(axName).thrPoint = thrPoint;
        sig.(axName).clThr = clThr;
    end

    %% -----------------------------------------------------
    % Plot
    %% -----------------------------------------------------
    figure('Color','w','Name',[areaName '_' stimname '_VAC_monkey_variance']);
    for k=1:nAxes
        axName = axisNamesOrdered{k};

        subplot(nAxes,1,k); hold on; box on
        plot(tvec, V_M1(:,k), 'k-', 'LineWidth', 2);
        plot(tvec, V_M2(:,k), 'k--','LineWidth', 1.5);
        yl = ylim;

        m = sig.(axName).mask;
        if any(m)
            ymark = yl(1) + 0.03*diff(yl);
            plot(tvec(m), ymark*ones(sum(m),1), '.', 'MarkerSize', 10);
        end

        title(sprintf('%s axis (VAC order): %s(solid) vs %s(dashed)', ...
            axName, monkeyName{1}, monkeyName{2}));
        ylabel('Var explained');
        if k==nAxes, xlabel('Time (ms)'); end
        legend({monkeyName{1}, monkeyName{2}, 'Sig (cluster)'},'Location','best');
    end

    %% -----------------------------------------------------
    % Save outputs
    %% -----------------------------------------------------
    OutMonkey = struct();
    OutMonkey.area     = areaName;
    OutMonkey.duration = durName;
    OutMonkey.stimname = stimname;
    OutMonkey.order    = REG_params_use;
    OutMonkey.M1proj   = M1proj;
    OutMonkey.M2proj   = M2proj;
    OutMonkey.V_M1     = V_M1;
    OutMonkey.V_M2     = V_M2;
    OutMonkey.sig      = sig;

    axisTag = strjoin(axisNamesOrdered,'_');

    SaveFileName = fullfile(OutputPath, ...
        sprintf('%s_%s_%s_%s_vs_%s_%s.mat', ...
        areaName, durName, stimname, monkeyName{1}, monkeyName{2}, axisTag));
    save(SaveFileName, 'OutMonkey','-v7.3');

    set(gcf,'paperpositionmode','auto');
    filename = fullfile(OutputPath, ...
        sprintf('%s_%s_%s_%s_vs_%s_%s.png', ...
        areaName, durName, stimname, monkeyName{1}, monkeyName{2}, axisTag));
    saveas(gcf,filename,'png');
    close;

end

end


% ============================================================
% HELPERS
% ============================================================
function [Vsub, TargetedPCAproj, tvec] = computeVar_subsetVAC(DataRun, PCA_params, idxKeep)
DataSub = DataRun;
for c = 1:numel(DataSub)
    DataSub(c).A = DataRun(c).A(:, idxKeep);
end

PCAp = PCA_params;
for k=1:numel(PCA_params.B)
    bk = PCA_params.B{k};
    PCAp.B{k} = bk(idxKeep,:);
end

[ProjectionRun,Sum] = TDR_time(DataSub, [], PCAp);
Vsub = Sum.varCaptEachTargetedPCT;
TargetedPCAproj = {ProjectionRun.TargetedPCAproj};

if isfield(DataSub(1),'times') && ~isempty(DataSub(1).times)
    tvec = DataSub(1).times(:);
else
    tvec = (1:size(Vsub,1))';
end
end

function [sigMask, thrPoint, clThr] = clusterPerm_twoSided(realCurve, nullCurves, alphaPoint, alphaClust)
absNull = abs(nullCurves(:));
thrPoint = quantile(absNull, 1-alphaPoint);

nPerm = size(nullCurves,2);
maxMass = zeros(nPerm,1);
for p=1:nPerm
    maxMass(p) = maxClusterMass_abs(nullCurves(:,p), thrPoint);
end
clThr = quantile(maxMass, 1-alphaClust);

sigMask = clusterMask_abs(realCurve, thrPoint, clThr);
end

function m = maxClusterMass_abs(x, thr)
ax = abs(x);
mask = ax > thr;
if ~any(mask), m = 0; return; end
d = diff([false; mask(:); false]);
s = find(d==1); e = find(d==-1)-1;
m = 0;
for i=1:numel(s)
    seg = ax(s(i):e(i)) - thr;
    m = max(m, sum(seg,'omitnan'));
end
end

function sigMask = clusterMask_abs(x, thr, clThr)
ax = abs(x);
mask = ax > thr;
sigMask = false(size(x));
if ~any(mask), return; end
d = diff([false; mask(:); false]);
s = find(d==1); e = find(d==-1)-1;
for i=1:numel(s)
    segMass = sum(ax(s(i):e(i)) - thr,'omitnan');
    if segMass > clThr
        sigMask(s(i):e(i)) = true;
    end
end
end