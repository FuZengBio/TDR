function tTDTR_timecourse_decoding()
clear; clc;
%% =========================================================
% unified paths
%% =========================================================
FilePath   = 'Z:\Users\TDR\Data\';

%% =========================================================
% unified settings
%% =========================================================
sigma = 30;
half_guas = sigma * 3;
dist = (-(half_guas):(half_guas))';
window = 1/(sigma*sqrt(2*pi)) * exp(-((dist.^2)/(2*sigma^2))) * 1000;

stim_type = [1,2];

% variance significance
nPermVar      = 1000;
alphaVarPoint = 0.05;
alphaVarClust = 0.05;

% decoding
trainFrac     = 0.7;
nIterReal     = 500;
nIterNull     = 500;   % keep helper logic unchanged
nPermDec      = 1000;
alphaCluster  = 0.05;

span = 10;

%% =========================================================
% unified run list (5 datasets)
%% =========================================================
datasets = { ...
    struct('area','VIP',  'dur','1s', 'mat','VIP_Data_1s_REG.mat',  'idMode','VAC'), ...
    struct('area','VIP',  'dur','2s', 'mat','VIP_Data_2s_REG.mat',  'idMode','VAC'), ...
    struct('area','MSTd', 'dur','1s', 'mat','MSTd_Data_1s_REG.mat', 'idMode','VAC'), ...
    struct('area','MSTd', 'dur','2s', 'mat','MSTd_Data_2s_REG.mat', 'idMode','VAC'), ...
    struct('area','PIVC', 'dur','1s', 'mat','PIVC_Data_1s_REG.mat', 'idMode','VAC') ...
    };

%% =========================================================
% run all
%% =========================================================
for k = 1:numel(datasets)
    ds = datasets{k};

    fprintf('\n====================================================\n');
    fprintf('Running %d/%d: %s %s\n', k, numel(datasets), ds.area, ds.dur);
    fprintf('====================================================\n');

    matfile = fullfile(FilePath, ds.mat);
    if ~exist(matfile, 'file')
        warning('Missing file: %s\nSkip.', matfile);
        continue;
    end

    run_one_dataset_vac_strictcv(matfile, ds.area, ds.dur, ds.idMode,...
        window, stim_type, ...
        nPermVar, alphaVarPoint, alphaVarClust, ...
        trainFrac, nIterReal, nIterNull, nPermDec, alphaCluster, span);
end

fprintf('\nAll datasets finished.\n');

end


%% ========================================================================
function run_one_dataset_vac_strictcv(matfile, areaName, durName, idMode,...
    window, stim_type, ...
    nPermVar, alphaVarPoint, alphaVarClust, ...
    trainFrac, nIterReal, nIterNull, nPermDec, alphaCluster, span)

load(matfile);

%% -----------------------------------------------------
% dataset-specific timing / crop / gaussian kernel
%% -----------------------------------------------------
switch durName
    case '2s'
        crop_idx = 101:2000;
        t1       = 101:2000;
        gaussArg = [0.13, 6, 2];

    case '1s'
        t1 = 101:1000;
        gaussArg = [0.13, 4, 1];

        if strcmp(areaName,'VIP') || strcmp(areaName,'MSTd')
            crop_idx = 1186:2085;
        elseif strcmp(areaName,'PIVC')
            crop_idx = 101:1000;
        else
            error('Unknown areaName for 1s data: %s', areaName);
        end

    otherwise
        error('Unknown durName: %s', durName);
end

REG_t = square_sum(t1, span) / span;

%% regression params
REG_params_stim   = {'choice','FUNangle','ONES'};
REG_choice_dim    = find(strcmp(REG_params_stim,'choice'));
REG_head_dim      = find(strcmp(REG_params_stim,'FUNangle')); %#ok<NASGU>

REG_params_kernel = {'VELangle','ACCangle','choice'};
REG_vel_dimK      = find(strcmp(REG_params_kernel,'VELangle'));
REG_acc_dimK      = find(strcmp(REG_params_kernel,'ACCangle'));
REG_choice_dimK   = find(strcmp(REG_params_kernel,'choice'));

%% kernel
kernel_del = 0;
[gausacc1, gausvel1, gauspos1] = genGaussian(gaussArg(1), gaussArg(2), gaussArg(3));

REG_acc_ideal1 = square_sum(gausacc1((t1)-kernel_del)', span) / span;
REG_acc_ideal1 = REG_acc_ideal1 ./ sqrt(REG_acc_ideal1 * REG_acc_ideal1');

REG_vel_ideal1 = square_sum(gausvel1((t1)-kernel_del)', span) / span;
REG_vel_ideal1 = REG_vel_ideal1 ./ sqrt(REG_vel_ideal1 * REG_vel_ideal1');

REG_pos_ideal1 = square_sum(gauspos1((t1)-kernel_del)', span) / span;
REG_pos_ideal1 = REG_pos_ideal1 ./ sqrt(REG_pos_ideal1 * REG_pos_ideal1');

kernel = [];
kernel(REG_choice_dimK,:) = REG_pos_ideal1;
kernel(REG_vel_dimK,:)    = REG_vel_ideal1;
kernel(REG_acc_dimK,:)    = REG_acc_ideal1;

%% choose axis order
switch upper(idMode)
    case 'VAC'   % [VEL ACC choice]
        id = [REG_vel_dimK REG_acc_dimK REG_choice_dimK];
    case 'VCA'   % [VEL choice ACC]
        id = [REG_vel_dimK REG_choice_dimK REG_acc_dimK];
    case 'CAV'
        id = [REG_choice_dimK REG_acc_dimK REG_vel_dimK];
    otherwise
        error('Unknown idMode: %s', idMode);
end

%% unified output folder
OutputPath = 'Z:\Users\TDR\results\tTDR_timecourse_decoding\';
if ~exist(OutputPath, 'dir')
    mkdir(OutputPath);
end

%% =========================================================
% run stim type
%% =========================================================
for ss = 1:length(stim_type)

    clear Test_REG_r Test_PCA_Data REG_bVesK_Test projection Summary

    %% -----------------------------------------------------
    % Part 1: SDF
    %% -----------------------------------------------------
    for nCell = 1:size(REG_raster_per_trial,1)
        temp_spike = REG_raster_per_trial{nCell,ss};
        sdf = [];

        if isempty(temp_spike)
            Test_REG_r{nCell,1} = [];
            continue;
        end

        for i = 1:size(temp_spike,2)
            spike = double(squeeze(temp_spike(:,i,:)));
            [ntrials,~] = size(spike);

            for ii = 1:ntrials
                convspike = conv(spike(ii,:), window, 'same');
                convspike = convspike(crop_idx);
                tempdata = downsample(convspike,10);
                sdf((i-1)*ntrials+ii,:) = tempdata; %#ok<AGROW>
            end
        end
        Test_REG_r{nCell,1} = sdf;
    end

    %% -----------------------------------------------------
    % Part 2: kernel regression fit
    %% -----------------------------------------------------
    REG_bVesK_Test = cell(length(REG_params_kernel),1);

    for Ni = 1:size(Test_REG_r,1)

        TMP_FVesK = [];
        TMP = [];
        tempdata = REG_Stim{Ni};

        if ~isempty(tempdata)
            TMP_Ves = (tempdata == 1);

            if sum(TMP_Ves) > 0
                for b = 1:length(REG_params_kernel)
                    if ss == 1
                        eval(sprintf('TMP(b,:)=REG_%s{Ni}(find(TMP_Ves));', REG_params_kernel{b}));
                    else
                        visGroup = Ni + size(REG_Stim,1);
                        eval(sprintf('TMP(b,:)=REG_%s{visGroup}(find(TMP_Ves));', REG_params_kernel{b}));
                    end
                end

                for trl = 1:sum(TMP_Ves)
                    TMP_FVesK = [TMP_FVesK repmat(TMP(:,trl),1,sum(REG_t<=max(t1))) .* kernel]; %#ok<AGROW>
                end

                temp_r = Test_REG_r{Ni}(TMP_Ves,1:sum(REG_t<=max(t1)));
                r = bsxfun(@minus,temp_r,mean(temp_r));

                TMP_bK = ((TMP_FVesK*(TMP_FVesK'))^-1) * TMP_FVesK * reshape(r',[],1);

                for b = 1:length(REG_params_kernel)
                    REG_bVesK_Test{b}(Ni,1) = TMP_bK(b);
                end
            else
                for b = 1:length(REG_params_kernel)
                    REG_bVesK_Test{b}(Ni,1) = NaN;
                end
            end
        else
            for b = 1:length(REG_params_kernel)
                REG_bVesK_Test{b}(Ni,1) = NaN;
            end
        end
    end

    %% -----------------------------------------------------
    % Part 3: PCA data
    %% -----------------------------------------------------
    REG_b_use = REG_bVesK_Test(id);

    for nCell = 1:size(REG_raster_per_trial,1)
        temp_spike = REG_raster_per_trial{nCell,ss};

        if isempty(temp_spike)
            continue;
        end

        for i = 1:size(temp_spike,1)
            spike = double(squeeze(temp_spike(i,:,:)));
            [ntrials,~] = size(spike);
            temp_mean = [];

            for ii = 1:ntrials
                convspike = conv(spike(ii,:), window, 'same');
                convspike = convspike(crop_idx);
                tempdata = downsample(convspike,10);
                temp_mean(ii,:) = tempdata; %#ok<AGROW>
            end

            Test_PCA_Data(1,i).A(:,nCell) = mean(temp_mean,1)';
            Test_PCA_Data(1,i).times = REG_t;
        end
    end

    %% -----------------------------------------------------
    % Part 4: TDR
    %% -----------------------------------------------------
    DataRun = Test_PCA_Data;
    PCA_params = [];
    PCA_params.B = REG_b_use;
    PCA_params.normalize = false;
    PCA_params.numPCs = 12;
    times = [];

    [projection, Summary] = TDR_time(DataRun, times, PCA_params); %#ok<ASGLU>

    if ss==1
        stimname = 'Vest';
    else
        stimname = 'Vis';
    end

    fprintf('\n[%s | %s %s] VAC-TDR: variance sig + decoding(sig) ...\n', ...
        stimname, areaName, durName);

    REG_params_use = REG_params_kernel(id);

    %% =========================
    % (A) Variance significance
    %% =========================
    realVar = struct();
    for k = 1:numel(REG_params_use)
        realVar.(REG_params_use{k}) = smooth(Summary.varCaptEachTargetedPCT(:,k),5);
    end

    choiceMap = detectChoiceMapping_vac(REG_choice, ss);
    Nstim = size(REG_Stim,1);

    [sigVar, ~] = VAC_varianceSig_clusterperm_labelshuffle( ...
        Test_REG_r, REG_Stim, ss, REG_params_kernel, id, kernel, REG_t, ...
        DataRun, PCA_params, ...
        choiceMap, ...
        nPermVar, alphaVarPoint, alphaVarClust, ...
        REG_VELangle, REG_ACCangle, REG_choice, Nstim);

    %% =========================
    % (B) Decoding + significance
    %% =========================
    idx_vel = find(id==REG_vel_dimK);
    idx_acc = find(id==REG_acc_dimK);
    idx_ch  = find(id==REG_choice_dimK);

    w_vel = normalizeAxisSafe_local(REG_b_use{idx_vel});
    w_acc = normalizeAxisSafe_local(REG_b_use{idx_acc});
    w_ch  = normalizeAxisSafe_local(REG_b_use{idx_ch});

    all_heads = [];
    for n=1:size(REG_FUNangle,1)
        if ~isempty(REG_FUNangle{n,ss})
            all_heads = [all_heads; REG_FUNangle{n,ss}(:)]; %#ok<AGROW>
        end
    end
    all_heads = all_heads(isfinite(all_heads));

    % keep original STRICTCV behavior: unsigned heading levels
    headLevels = unique(abs(all_heads))';
    headLevels = sort(headLevels);

    cvDec = crossValidatedDecoding_only_local( ...
        Test_REG_r, REG_FUNangle, REG_choice, ...
        w_vel, w_acc, w_ch, ...
        ss, choiceMap, headLevels, ...
        trainFrac, nIterReal, nPermDec, alphaCluster);

    acc_vel = cvDec.acc_vel;
    acc_acc = cvDec.acc_acc;
    acc_ch  = cvDec.acc_ch;

    sig_vel = cvDec.sig_vel;
    sig_acc = cvDec.sig_acc;
    sig_ch  = clusterPermSig_fromChance_binary(acc_ch, 0.5, alphaCluster);

    %% =========================
    % Save
    %% =========================
    KernelData = struct();
    KernelData.area   = areaName;
    KernelData.dur    = durName;
    KernelData.order  = REG_params_use;
    KernelData.Variance.real = realVar;
    KernelData.Variance.sig  = sigVar;

    KernelData.Decoding.Vel    = acc_vel;
    KernelData.Decoding.VelSig = sig_vel;
    KernelData.Decoding.Acc    = acc_acc;
    KernelData.Decoding.AccSig = sig_acc;
    KernelData.Decoding.Ch     = acc_ch;
    KernelData.Decoding.ChSig  = sig_ch;

    param_str = strjoin(REG_params_use, '_');
    save_name_str = sprintf('%s_%s_%s_VAC_Decoding_%s.mat', ...
        areaName, durName, stimname, param_str);
    SaveFileName = fullfile(OutputPath, save_name_str);

    save(SaveFileName, 'KernelData');
    fprintf('Data saved to: %s\n', SaveFileName);

    %% optional plots
    figure('Color','w','Name',[areaName '_' durName '_' stimname ' VAC Decoding']);
    vecY = @(y, m) y * ones(sum(m), 1);

    subplot(3,1,1); hold on; box on
    plot(REG_t, acc_vel, 'k', 'LineWidth', 2);
    yline(1/numel(headLevels), 'k--'); ylim([0 1]); yl=ylim;
    if any(sig_vel), plot(REG_t(sig_vel), vecY(yl(1)+0.05, sig_vel), 'k.', 'MarkerSize',10); end
    title('Heading decoding on VELOCITY axis'); ylabel('Accuracy');

    subplot(3,1,2); hold on; box on
    plot(REG_t, acc_acc, 'k', 'LineWidth', 2);
    yline(1/numel(headLevels), 'k--'); ylim([0 1]); yl=ylim;
    if any(sig_acc), plot(REG_t(sig_acc), vecY(yl(1)+0.05, sig_acc), 'k.', 'MarkerSize',10); end
    title('Heading decoding on ACCELERATION axis'); ylabel('Accuracy');

    subplot(3,1,3); hold on; box on
    plot(REG_t, acc_ch, 'k', 'LineWidth', 2);
    yline(0.5, 'k--'); ylim([0.3 1]); yl=ylim;
    if any(sig_ch), plot(REG_t(sig_ch), vecY(yl(1)+0.05, sig_ch), 'k.', 'MarkerSize',10); end
    title('Choice decoding on CHOICE axis'); ylabel('Accuracy'); xlabel('Time [ms]');

    figure('Color','w','Name',[areaName '_' durName '_' stimname ' VAC Variance']);
    for k=1:numel(REG_params_use)
        vn = REG_params_use{k};
        subplot(numel(REG_params_use),1,k); hold on; box on
        plot(REG_t, realVar.(vn), 'k', 'LineWidth', 2);
        yl = ylim;
        m = sigVar.(vn);
        if any(m)
            plot(REG_t(m), (yl(1)+0.05*diff(yl))*ones(sum(m),1), 'k.', 'MarkerSize',10);
        end
        title(['VAC variance on ' vn ' axis']); ylabel('Var%');
    end
    xlabel('Time [ms]');

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------ HELPERS ---------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Choice mapping detection (-1/+1 or 1/2) ---
function mp = detectChoiceMapping_vac(REG_choice, ss)
allC = [];
N = size(REG_choice,1);
for n=1:N
    c = REG_choice{n,ss};
    if ~isempty(c), allC = [allC; c(:)]; end %#ok<AGROW>
end
allC = allC(isfinite(allC));
u = unique(allC); u = sort(u);

if isequal(u, [-1; 1]) || isequal(u, [-1 1]')
    mp.name = '-1/+1';
    mp.toPM = @(c) c;
elseif isequal(u, [1;2]) || isequal(u, [1 2]')
    mp.name = '1/2 -> -1/+1';
    mp.toPM = @(c) (c==2) - (c==1);
else
    error('Unknown choice labels. Unique values: %s', mat2str(u'));
end
end

% -------------------------------------------------------------------------
% VAC variance significance:
%   - shuffle trial labels within each neuron (kernel predictors)
%   - recompute REG_b_use (VAC axes)
%   - run TDR_time -> get varCaptEachTargetedPCT
%   - cluster permutation on time
% -------------------------------------------------------------------------
function [sigVar, info] = VAC_varianceSig_clusterperm_labelshuffle( ...
    Test_REG_r, REG_Stim, ss, REG_params_kernel, id, kernel, REG_t, ...
    DataRun, PCA_params_real, ...
    choiceMap, ...
    nPerm, alphaPoint, alphaClust, ...
    REG_VELangle, REG_ACCangle, REG_choice, Nstim)

[~, SumReal] = TDR_time(DataRun, [], PCA_params_real);
Vreal = SumReal.varCaptEachTargetedPCT;   % T x nAxes
T = size(Vreal,1);
nAxes = size(Vreal,2);
REG_params_use = REG_params_kernel(id);

nullV = nan(T, nPerm, nAxes);
Ncells = numel(Test_REG_r);

for p = 1:nPerm

    REG_bVesK_null = cell(numel(REG_params_kernel),1);
    for b = 1:numel(REG_params_kernel)
        REG_bVesK_null{b} = nan(Ncells,1);
    end

    for Ni = 1:Ncells

        Xresp = Test_REG_r{Ni};
        if isempty(Xresp)
            continue;
        end

        tempStim = REG_Stim{Ni};
        if isempty(tempStim)
            continue;
        end

        % -------- choose label source --------
        if ss == 1
            velCell = REG_VELangle{Ni};
            accCell = REG_ACCangle{Ni};
            choCell = REG_choice{Ni};
        else
            visGroup = Ni + Nstim;
            velCell = REG_VELangle{visGroup};
            accCell = REG_ACCangle{visGroup};
            choCell = REG_choice{visGroup};
        end

        % -------- force all trial vectors to same length --------
        nResp = size(Xresp,1);
        nStim = numel(tempStim);
        nVel  = numel(velCell);
        nAcc  = numel(accCell);
        nCho  = numel(choCell);

        nAll = min([nResp, nStim, nVel, nAcc, nCho]);

        if isempty(nAll) || nAll < 2
            continue;
        end

        Xresp    = Xresp(1:nAll, 1:T);
        tempStim = tempStim(1:nAll);
        velCell  = velCell(1:nAll);
        accCell  = accCell(1:nAll);
        choCell  = choCell(1:nAll);

        % -------- select vestibular trials (same as original code logic) --------
        idxUse = find(tempStim == 1);
        if numel(idxUse) < 2
            continue;
        end

        Xresp = Xresp(idxUse, :);
        velUse = velCell(idxUse);
        accUse = accCell(idxUse);
        choUse = choCell(idxUse);

        nTr = size(Xresp,1);
        if nTr < 2
            continue;
        end

        % -------- build TMP in the same row order as REG_params_kernel --------
        TMP = nan(numel(REG_params_kernel), nTr);

        for b = 1:numel(REG_params_kernel)
            nm = REG_params_kernel{b};
            switch nm
                case 'VELangle'
                    TMP(b,:) = velUse(:)';
                case 'ACCangle'
                    TMP(b,:) = accUse(:)';
                case 'choice'
                    TMP(b,:) = choUse(:)';
                otherwise
                    error('Unknown kernel param: %s', nm);
            end
        end

        % -------- convert / shuffle labels --------
        chRow = find(strcmp(REG_params_kernel,'choice'));
        if ~isempty(chRow)
            TMP(chRow,:) = choiceMap.toPM(TMP(chRow,:));
            TMP(chRow,:) = TMP(chRow, randperm(nTr));
        end

        velRow = find(strcmp(REG_params_kernel,'VELangle'));
        if ~isempty(velRow)
            TMP(velRow,:) = TMP(velRow, randperm(nTr));
        end

        accRow = find(strcmp(REG_params_kernel,'ACCangle'));
        if ~isempty(accRow)
            TMP(accRow,:) = TMP(accRow, randperm(nTr));
        end

        % -------- build design matrix: [nKernel x (nTr*T)] --------
        TMP_FVesK = zeros(numel(REG_params_kernel), nTr*T);
        for trl = 1:nTr
            cols = (trl-1)*T + (1:T);
            TMP_FVesK(:, cols) = repmat(TMP(:,trl), 1, T) .* kernel(:,1:T);
        end

        % -------- response vector: [(nTr*T) x 1] --------
        r = bsxfun(@minus, Xresp, mean(Xresp,1,'omitnan'));
        y = reshape(r', [], 1);

        % -------- final safety check --------
        if size(TMP_FVesK,2) ~= numel(y)
            warning('Ni=%d perm=%d skipped: design/response mismatch (%d vs %d)', ...
                Ni, p, size(TMP_FVesK,2), numel(y));
            continue;
        end

        % -------- robust solve --------
        TMP_bK = pinv(TMP_FVesK * TMP_FVesK') * TMP_FVesK * y;

        for b = 1:numel(REG_params_kernel)
            REG_bVesK_null{b}(Ni,1) = TMP_bK(b);
        end
    end

    Bnull = REG_bVesK_null(id);
    PCA_params_null = PCA_params_real;
    PCA_params_null.B = Bnull;

    [~, SumNull] = TDR_time(DataRun, [], PCA_params_null);
    nullV(:,p,:) = SumNull.varCaptEachTargetedPCT;
end

sigVar = struct();
info = struct();
info.nullV = nullV;
info.Vreal = Vreal;

for a = 1:nAxes
    vn = REG_params_use{a};
    realCurve = Vreal(:,a);
    nullCurve = squeeze(nullV(:,:,a));   % T x nPerm

    pooledNull = nullCurve(:);
    pooledNull = pooledNull(isfinite(pooledNull));

    if isempty(pooledNull)
        sigVar.(vn) = false(T,1);
        continue;
    end

    pointThr = quantile(pooledNull, 1 - alphaPoint);

    realMask = realCurve > pointThr;
    realClusters = find_clusters_1d(realMask);
    realMass = cluster_mass_1d(realCurve, realClusters, pointThr);

    nullMaxMass = zeros(nPerm,1);
    for pp = 1:nPerm
        cur = nullCurve(:,pp);
        curMask = cur > pointThr;
        curClusters = find_clusters_1d(curMask);
        curMass = cluster_mass_1d(cur, curClusters, pointThr);

        if isempty(curMass)
            nullMaxMass(pp) = 0;
        else
            nullMaxMass(pp) = max(curMass);
        end
    end

    clustThr = quantile(nullMaxMass, 1 - alphaClust);
    sigMask = false(T,1);

    for c = 1:numel(realClusters)
        idx = realClusters{c};
        if realMass(c) > clustThr
            sigMask(idx) = true;
        end
    end

    sigVar.(vn) = sigMask;
end

end

function m = maxClusterMass(x, thr)
mask = x > thr;
if ~any(mask), m = 0; return; end
d = diff([false; mask(:); false]);
s = find(d==1); e = find(d==-1)-1;
m = 0;
for i=1:numel(s)
    seg = x(s(i):e(i)) - thr;
    m = max(m, sum(seg,'omitnan'));
end
end

function sigMask = clusterMask_real(x, thr, clThr)
mask = x > thr;
sigMask = false(size(x));
if ~any(mask), return; end
d = diff([false; mask(:); false]);
s = find(d==1); e = find(d==-1)-1;
for i=1:numel(s)
    segMass = sum(x(s(i):e(i)) - thr,'omitnan');
    if segMass > clThr
        sigMask(s(i):e(i)) = true;
    end
end
end

% -------------------------------------------------------------------------
% build decoding pools with auto choice mapping + train/test split
% -------------------------------------------------------------------------
function sigMask = clusterPermSig_fromNull(accReal, nullAcc, alpha)
thr = quantile(nullAcc(:), 0.95);
nPerm = size(nullAcc,2);
maxMass = zeros(nPerm,1);
for p=1:nPerm
    maxMass(p) = maxClusterMass(nullAcc(:,p), thr);
end
clThr = quantile(maxMass, 1-alpha);
sigMask = clusterMask_real(accReal, thr, clThr);
end

function out = crossValidatedDecoding_only_local( ...
    Test_REG_r, REG_FUNangle, REG_choice, ...
    w_vel, w_acc, w_ch, ...
    ss, choiceMap, headLevels, ...
    trainFrac, nIterReal, nPermDec, alphaCluster)

Ncells = numel(Test_REG_r);

T = [];
for n=1:Ncells
    if ~isempty(Test_REG_r{n})
        T = size(Test_REG_r{n},2);
        break
    end
end
if isempty(T)
    error('No valid Test_REG_r found.')
end

accVel_all = nan(T, nIterReal);
accAcc_all = nan(T, nIterReal);
accCh_all  = nan(T, nIterReal);

for it = 1:nIterReal
    splitInfo = buildTrainTestSplitPerNeuron_noStim_local( ...
        Test_REG_r, REG_FUNangle, REG_choice, ss, trainFrac);

    [Hp_tr, Hp_te, Cp_tr, Cp_te] = buildPoolsFromSplit_noStim_local( ...
        Test_REG_r, REG_FUNangle, REG_choice, ss, choiceMap, headLevels, splitInfo);

    tmp = decode1D_template_pseudopop_robust_local(Hp_tr, Hp_te, w_vel, 1);
    if isempty(tmp), tmp = nan(T,1); end
    accVel_all(:,it) = tmp;

    tmp = decode1D_template_pseudopop_robust_local(Hp_tr, Hp_te, w_acc, 1);
    if isempty(tmp), tmp = nan(T,1); end
    accAcc_all(:,it) = tmp;

    tmp = decode1D_template_pseudopop_robust_local(Cp_tr, Cp_te, w_ch, 1);
    if isempty(tmp), tmp = nan(T,1); end
    accCh_all(:,it) = tmp;
end

acc_vel = mean(accVel_all,2,'omitnan');
acc_acc = mean(accAcc_all,2,'omitnan');
acc_ch  = mean(accCh_all,2,'omitnan');

nullVel = nan(T,nPermDec);
nullAcc = nan(T,nPermDec);
nullCh  = nan(T,nPermDec);

for p = 1:nPermDec
    splitInfo = buildTrainTestSplitPerNeuron_noStim_local( ...
        Test_REG_r, REG_FUNangle, REG_choice, ss, trainFrac);

    % [Hp_tr, Hp_te, Cp_tr, Cp_te] = buildPoolsFromSplit_noStim_local( ...
    %     Test_REG_r, REG_FUNangle, REG_choice, ss, choiceMap, headLevels, splitInfo);
    [Hp_tr, Hp_te, Cp_tr, Cp_te] = buildPoolsFromSplit_noStim_unsigned_local( ...
    Test_REG_r, REG_FUNangle, REG_choice, ss, choiceMap, headLevels, splitInfo);

    [Hp_tr_null, Hp_te_null] = shuffleHeadingPools_local(Hp_tr, Hp_te);
    [Cp_tr_null, Cp_te_null] = shuffleChoicePools_local(Cp_tr, Cp_te);

    tmp = decode1D_template_pseudopop_robust_local(Hp_tr_null, Hp_te_null, w_vel, 1);
    if isempty(tmp), tmp = nan(T,1); end
    nullVel(:,p) = tmp;

    tmp = decode1D_template_pseudopop_robust_local(Hp_tr_null, Hp_te_null, w_acc, 1);
    if isempty(tmp), tmp = nan(T,1); end
    nullAcc(:,p) = tmp;

    tmp = decode1D_template_pseudopop_robust_local(Cp_tr_null, Cp_te_null, w_ch, 1);
    if isempty(tmp), tmp = nan(T,1); end
    nullCh(:,p) = tmp;
end

sig_vel = clusterPermSig_fromNull(acc_vel, nullVel, alphaCluster);
sig_acc = clusterPermSig_fromNull(acc_acc, nullAcc, alphaCluster);
sig_ch  = clusterPermSig_fromNull(acc_ch,  nullCh,  alphaCluster);

out = struct();
out.acc_vel = acc_vel;
out.acc_acc = acc_acc;
out.acc_ch  = acc_ch;
out.sig_vel = sig_vel;
out.sig_acc = sig_acc;
out.sig_ch  = sig_ch;
end

function splitInfo = buildTrainTestSplitPerNeuron_noStim_local( ...
    Test_REG_r, REG_FUNangle, REG_choice, ss, trainFrac)

N = numel(Test_REG_r);
splitInfo = cell(N,1);

for n = 1:N
    Y = Test_REG_r{n};
    if isempty(Y)
        splitInfo{n} = [];
        continue
    end

    h = REG_FUNangle{n,ss};
    c = REG_choice{n,ss};

    if isempty(h) || isempty(c)
        splitInfo{n} = [];
        continue
    end

    L = min([size(Y,1), numel(h), numel(c)]);
    h = h(1:L);
    c = c(1:L);

    if L < 4
        splitInfo{n} = [];
        continue
    end

    trIdx = [];
    teIdx = [];

    heads = unique(h(isfinite(h)))';
    for hh = heads
        idxH = find(abs(h - hh) < 1e-12);

        if numel(idxH) == 1
            trIdx = [trIdx; idxH(:)];
        else
            [trLocal, teLocal] = splitIdx_local(numel(idxH), trainFrac);
            trIdx = [trIdx; idxH(trLocal(:))];
            teIdx = [teIdx; idxH(teLocal(:))];
        end
    end

    trIdx = unique(trIdx(:));
    teIdx = unique(teIdx(:));

    info = struct();
    info.allIdx  = (1:L)';
    info.trIdx   = trIdx;
    info.teIdx   = teIdx;
    info.h_train = h(trIdx);
    info.h_test  = h(teIdx);
    info.c_train = c(trIdx);
    info.c_test  = c(teIdx);

    splitInfo{n} = info;
end
end

function [Hp_tr, Hp_te, Cp_tr, Cp_te] = buildPoolsFromSplit_noStim_local( ...
    Test_REG_r, REG_FUNangle, REG_choice, ss, choiceMap, headLevels, splitInfo)

N = numel(Test_REG_r);
K = numel(headLevels);

Hp_tr = cell(N,K); Hp_te = cell(N,K);
Cp_tr = cell(N,2); Cp_te = cell(N,2);

for n = 1:N
    Y = Test_REG_r{n};
    h = REG_FUNangle{n,ss};
    c = REG_choice{n,ss};

    sp = splitInfo{n};
    if isempty(Y) || isempty(h) || isempty(c) || isempty(sp), continue; end

    L = min([size(Y,1), numel(h), numel(c)]);
    Y = Y(1:L,:);
    h = h(1:L);
    c = c(1:L);
    c_pm = choiceMap.toPM(c);

    trIdx = sp.trIdx(sp.trIdx <= L);
    teIdx = sp.teIdx(sp.teIdx <= L);

    for k = 1:K
        idx_tr = intersect(trIdx(:), find(abs(h-headLevels(k))<1e-12));
        idx_te = intersect(teIdx(:), find(abs(h-headLevels(k))<1e-12));
        Hp_tr{n,k} = Y(idx_tr,:);
        Hp_te{n,k} = Y(idx_te,:);
    end

    idxL_tr = intersect(trIdx(:), find(c_pm==-1));
    idxR_tr = intersect(trIdx(:), find(c_pm==+1));
    idxL_te = intersect(teIdx(:), find(c_pm==-1));
    idxR_te = intersect(teIdx(:), find(c_pm==+1));

    Cp_tr{n,1} = Y(idxL_tr,:);
    Cp_tr{n,2} = Y(idxR_tr,:);
    Cp_te{n,1} = Y(idxL_te,:);
    Cp_te{n,2} = Y(idxR_te,:);
end
end

function acc = decode1D_template_pseudopop_robust_local(poolsTr, poolsTe, axisVec, nIter)

axisVec = axisVec(:);
[N,K] = size(poolsTr);

T = [];
for n=1:N
    for k=1:K
        if ~isempty(poolsTr{n,k})
            T = size(poolsTr{n,k},2);
            break
        end
    end
    if ~isempty(T), break; end
end
if isempty(T)
    acc = nan(0,1);
    return
end

tmpl = nan(K,T);

for k=1:K
    x = zeros(1,T);
    used = 0;
    for n=1:N
        tr = poolsTr{n,k};
        if isempty(tr) || ~isfinite(axisVec(n))
            continue
        end
        mu = mean(tr,1,'omitnan');
        if any(~isfinite(mu)), continue; end
        x = x + axisVec(n) * mu;
        used = used + 1;
    end
    if used > 0
        tmpl(k,:) = x;
    end
end

acc_it = nan(T,nIter);

for it=1:nIter
    X = nan(K,T);

    for k=1:K
        x = zeros(1,T);
        used = 0;
        for n=1:N
            te = poolsTe{n,k};
            if isempty(te) || ~isfinite(axisVec(n))
                continue
            end
            idx = randi(size(te,1));
            y = te(idx,:);
            if any(~isfinite(y)), continue; end
            x = x + axisVec(n) * y;
            used = used + 1;
        end
        if used > 0
            X(k,:) = x;
        end
    end

    a = nan(T,1);
    for t=1:T
        if any(~isfinite(tmpl(:,t))) || any(~isfinite(X(:,t)))
            a(t)=NaN;
            continue
        end

        correct_ct = 0;
        valid_ct = 0;
        for k=1:K
            d = abs(X(k,t)-tmpl(:,t));
            if any(~isfinite(d)), continue; end
            [~,pred] = min(d);
            correct_ct = correct_ct + double(pred==k);
            valid_ct = valid_ct + 1;
        end

        if valid_ct > 0
            a(t) = correct_ct / valid_ct;
        end
    end
    acc_it(:,it) = a;
end

acc = mean(acc_it,2,'omitnan');
end

function [ptr, pte] = shuffleHeadingPools_local(poolsTr, poolsTe)
[N,K] = size(poolsTr);
ptr = poolsTr; 
pte = poolsTe;
for n=1:N
    perm = randperm(K);
    ptr(n,:) = poolsTr(n,perm);
    pte(n,:) = poolsTe(n,perm);
end
end

function [ptr, pte] = shuffleChoicePools_local(poolsTr, poolsTe)
N = size(poolsTr,1);
ptr = poolsTr; 
pte = poolsTe;
for n=1:N
    if rand < 0.5
        ptr{n,1} = poolsTr{n,2}; ptr{n,2} = poolsTr{n,1};
        pte{n,1} = poolsTe{n,2}; pte{n,2} = poolsTe{n,1};
    end
end
end

function w = normalizeAxisSafe_local(w)
w = w(:);
w(isnan(w)) = 0;
nw = norm(w);
if nw > 0
    w = w / nw;
end
end

function [trIdx, teIdx] = splitIdx_local(n, frac)
if n <= 1
    trIdx = 1:n;
    teIdx = [];
    return
elseif n == 2
    trIdx = 1;
    teIdx = 2;
    return
end

rp = randperm(n);

ntr = round(frac*n);
ntr = max(1, min(n-1, ntr));

trIdx = rp(1:ntr);
teIdx = rp(ntr+1:end);

if isempty(teIdx)
    teIdx = trIdx(end);
    trIdx = trIdx(1:end-1);
end
if isempty(trIdx)
    trIdx = teIdx(1);
    teIdx = teIdx(2:end);
end
end

function sigMask = clusterPermSig_fromChance_binary(accReal, chanceLevel, alpha)
% simple cluster test against fixed chance level for binary decoding
thr = chanceLevel;

% build null cluster masses by sign-flipping around chance
T = numel(accReal);
nPerm = 1000;
maxMass = zeros(nPerm,1);

x = accReal - chanceLevel;
x(~isfinite(x)) = 0;

for p = 1:nPerm
    sgn = sign(randn(T,1));   % random sign flip
    xp = x .* sgn;
    maxMass(p) = maxClusterMass(xp + chanceLevel, thr);
end

clThr = quantile(maxMass, 1-alpha);
sigMask = clusterMask_real(accReal, thr, clThr);
end

function [Hp_tr, Hp_te, Cp_tr, Cp_te] = buildPoolsFromSplit_noStim_unsigned_local( ...
    Test_REG_r, REG_FUNangle, REG_choice, ss, choiceMap, headLevels, splitInfo)

% Unsigned-heading version:
% heading pools are grouped by abs(heading), e.g.
%   -16 and +16 go into the same class
% choice pools remain unchanged

N = numel(Test_REG_r);
K = numel(headLevels);

Hp_tr = cell(N,K); Hp_te = cell(N,K);
Cp_tr = cell(N,2); Cp_te = cell(N,2);

for n = 1:N
    Y = Test_REG_r{n};
    h = REG_FUNangle{n,ss};
    c = REG_choice{n,ss};

    sp = splitInfo{n};
    if isempty(Y) || isempty(h) || isempty(c) || isempty(sp), continue; end

    L = min([size(Y,1), numel(h), numel(c)]);
    Y = Y(1:L,:);
    h = h(1:L);
    c = c(1:L);
    c_pm = choiceMap.toPM(c);

    trIdx = sp.trIdx(sp.trIdx <= L);
    teIdx = sp.teIdx(sp.teIdx <= L);

    absH = abs(h);

    % ----------------------------
    % unsigned heading pools
    % ----------------------------
    for k = 1:K
        idx_tr = intersect(trIdx(:), find(abs(absH - headLevels(k)) < 1e-12));
        idx_te = intersect(teIdx(:), find(abs(absH - headLevels(k)) < 1e-12));

        Hp_tr{n,k} = Y(idx_tr,:);
        Hp_te{n,k} = Y(idx_te,:);
    end

    % ----------------------------
    % choice pools unchanged
    % ----------------------------
    idxL_tr = intersect(trIdx(:), find(c_pm==-1));
    idxR_tr = intersect(trIdx(:), find(c_pm==+1));
    idxL_te = intersect(teIdx(:), find(c_pm==-1));
    idxR_te = intersect(teIdx(:), find(c_pm==+1));

    Cp_tr{n,1} = Y(idxL_tr,:);
    Cp_tr{n,2} = Y(idxR_tr,:);
    Cp_te{n,1} = Y(idxL_te,:);
    Cp_te{n,2} = Y(idxR_te,:);
end
end

function clusters = find_clusters_1d(mask)
% Find contiguous true segments in a logical vector
%
% Input:
%   mask : logical or numeric vector
%
% Output:
%   clusters : cell array, each cell is a vector of indices

mask = mask(:)' ~= 0;

d = diff([false, mask, false]);
st = find(d == 1);
ed = find(d == -1) - 1;

clusters = cell(numel(st), 1);
for i = 1:numel(st)
    clusters{i} = st(i):ed(i);
end
end

function mass = cluster_mass_1d(curve, clusters, thr)
% Compute 1D cluster mass above threshold
%
% Input:
%   curve    : [T x 1] or [1 x T]
%   clusters : output from find_clusters_1d
%   thr      : scalar threshold
%
% Output:
%   mass     : [nClusters x 1]

curve = curve(:);
mass = zeros(numel(clusters), 1);

for i = 1:numel(clusters)
    idx = clusters{i};
    mass(i) = sum(curve(idx) - thr, 'omitnan');
end
end