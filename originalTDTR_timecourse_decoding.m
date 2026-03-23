function Run_All_HC_TDR_Decoding()

clear; clc;

%% =========================================================
% paths
%% =========================================================
FilePath   = 'Z:\Users\TDR\Data\';

%% =========================================================
% unified settings
%% =========================================================
sigma     = 30;
half_guas = sigma*3;
dist      = (-(half_guas):(half_guas))';
window    = 1/(sigma*sqrt(2*pi))*exp(-((dist.^2)/(2*sigma^2)))*1000;

stim_type = [1,2];   % 1=Vest, 2=Vis

reverse   = 0;
orderStr  = 'HC';

% decoding
trainFrac      = 0.7;
nIterReal      = 500;
nPermDec       = 1000;
alphaCluster   = 0.05;

% variance significance
nPermVar       = 1000;
alphaVarPoint  = 0.05;
alphaVarClust  = 0.05;

span = 10;

%% =========================================================
% run list (5 datasets)
%% =========================================================
datasets = { ...
    struct('area','VIP',  'dur','1s', 'mat','VIP_Data_1s_REG.mat'), ...
    struct('area','VIP',  'dur','2s', 'mat','VIP_Data_2s_REG.mat'), ...
    struct('area','MSTd', 'dur','1s', 'mat','MSTd_Data_1s_REG.mat'), ...
    struct('area','MSTd', 'dur','2s', 'mat','MSTd_Data_2s_REG.mat'), ...
    struct('area','PIVC', 'dur','1s', 'mat','PIVC_Data_1s_REG.mat') ...
    };

for k = 1:numel(datasets)

    ds = datasets{k};
    fprintf('\n====================================================\n');
    fprintf('Running %d/%d: %s %s\n', k, numel(datasets), ds.area, ds.dur);
    fprintf('====================================================\n');

    matfile = fullfile(FilePath, ds.mat);
    if ~exist(matfile, 'file')
        warning('File not found: %s. Skip.', matfile);
        continue;
    end

    run_one_dataset_HC(matfile, ds.area, ds.dur,...
        window, stim_type, reverse, orderStr, ...
        trainFrac, nIterReal, nPermDec, alphaCluster, ...
        nPermVar, alphaVarPoint, alphaVarClust, span);

end

fprintf('\nAll HC decoding finished.\n');

end

%% =========================================================
function run_one_dataset_HC(matfile, areaName, durName,...
    window, stim_type, reverse, orderStr, ...
    trainFrac, nIterReal, nPermDec, alphaCluster, ...
    nPermVar, alphaVarPoint, alphaVarClust, span)

load(matfile);

%% time window
switch durName
    case '2s'
        t1       = 101:2000;
        crop_idx = 101:2000;

    case '1s'
        t1 = 101:1000;

        if strcmp(areaName,'VIP') || strcmp(areaName,'MSTd')
            crop_idx = 1186:2085;
        elseif strcmp(areaName,'PIVC')
            crop_idx = 101:1000;
        else
            error('Unknown areaName: %s', areaName);
        end

    otherwise
        error('Unknown durName: %s', durName);
end

REG_t = square_sum(t1,span)/span;
T     = numel(REG_t);

REG_params_stim = {'choice','FUNangle','ONES'};
REG_choice_dim  = find(strcmp(REG_params_stim,'choice'));
REG_head_dim    = find(strcmp(REG_params_stim,'FUNangle'));

%% output folder
OutputPath = 'Z:\Users\TDR\results\originalTDR_timecourse_decoding\';
if ~exist(OutputPath,'dir')
    mkdir(OutputPath);
end

for ss = 1:numel(stim_type)

    if ss==1
        stimname='Vest';
    else
        stimname='Vis';
    end
    fprintf('\n=========== %s | %s %s ===========\n', stimname, areaName, durName);

    Ncells = size(REG_raster_per_trial,1);
    nAngle = size(REG_raster_per_trial{1,ss},1);

    %% -----------------------------------------------------
    % Part 1: Build trial-wise SDF
    %% -----------------------------------------------------
    Test_REG_r = cell(Ncells,1);
    for n=1:Ncells
        tmp = REG_raster_per_trial{n,ss};
        sdf = [];

        if isempty(tmp)
            Test_REG_r{n} = [];
            continue;
        end

        for i=1:size(tmp,2)
            blk = double(squeeze(tmp(:,i,:)));
            for a=1:size(blk,1)
                convspk = conv(blk(a,:),window,'same');
                convspk = convspk(crop_idx);
                sdf(end+1,:) = downsample(convspk,10); %#ok<AGROW>
            end
        end
        Test_REG_r{n} = sdf;
    end

    %% -----------------------------------------------------
    % Part 2: Build TDR regression betas
    %% -----------------------------------------------------
    if ~reverse
        id = [REG_head_dim REG_choice_dim];
        axisOrder = {'heading','choice'};
    else
        id = [REG_choice_dim REG_head_dim];
        axisOrder = {'choice','heading'};
    end %#ok<NASGU>

    choiceMap = detectChoiceMapping(REG_choice, ss);
    fprintf('Choice coding detected: %s\n', choiceMap.name);

    REG_bT = nan(Ncells,T,numel(REG_params_stim));

    for n=1:Ncells
        h = REG_FUNangle{n,ss};
        c = REG_choice{n,ss};

        if isempty(h) || isempty(c), continue; end
        if numel(h) ~= numel(c), continue; end

        c_pm = choiceMap.toPM(c);

        if ~reverse
            F = [h(:), c_pm(:), ones(length(h),1)];
        else
            F = [c_pm(:), h(:), ones(length(h),1)];
        end

        X = F;
        Y = Test_REG_r{n};

        if isempty(Y) || size(Y,2) < T, continue; end
        if size(X,1) ~= size(Y,1), continue; end

        B = X \ Y(:,1:T);   % 3 x T
        REG_bT(n,:,id(1)) = B(1,:);
        REG_bT(n,:,id(2)) = B(2,:);
    end

    B_heading = squeeze(REG_bT(:,:,REG_head_dim));
    B_choice  = squeeze(REG_bT(:,:,REG_choice_dim));

    %% -----------------------------------------------------
    % Part 3: Condition means for TDR_time projection
    %% -----------------------------------------------------
    Test_PCA_Data = struct();
    for c=1:nAngle
        Test_PCA_Data(1,c).A = nan(T,Ncells);
        Test_PCA_Data(1,c).times = REG_t;
    end

    for n=1:Ncells
        tmp = REG_raster_per_trial{n,ss};
        if isempty(tmp), continue; end

        for c=1:nAngle
            spk = double(squeeze(tmp(c,:,:)));
            tr  = size(spk,1);
            buf = nan(tr,T);

            for k=1:tr
                cs = conv(spk(k,:),window,'same');
                cs = cs(crop_idx);
                buf(k,:) = downsample(cs,10);
            end
            Test_PCA_Data(1,c).A(:,n) = mean(buf,1,'omitnan')';
        end
    end

    %% -----------------------------------------------------
    % Part 4: TDR_time
    %% -----------------------------------------------------
    PCA_params = struct();
    PCA_params.normalize = false;
    PCA_params.numPCs = 12;

    if strcmpi(orderStr,'HC')
        PCA_params.B = {B_heading, B_choice};
    elseif strcmpi(orderStr,'CH')
        PCA_params.B = {B_choice, B_heading};
    else
        error('orderStr must be HC or CH');
    end

    [projection, Summary] = TDR_time(Test_PCA_Data, [], PCA_params);

    if strcmpi(orderStr,'HC')
        varMap.heading = Summary.varCaptEachTargetedPCT(:,1);
        varMap.choice  = Summary.varCaptEachTargetedPCT(:,2);
    else
        varMap.choice  = Summary.varCaptEachTargetedPCT(:,1);
        varMap.heading = Summary.varCaptEachTargetedPCT(:,2);
    end

    %% -----------------------------------------------------
    % Part 5: Variance significance
    %% -----------------------------------------------------
    fprintf('[%s] Variance significance...\n', stimname);

    sigVarMap = varianceSig_labelShuffle_axes_clusterperm_TDRonly( ...
        Test_REG_r, REG_FUNangle, REG_choice, ss, choiceMap, ...
        Test_PCA_Data, PCA_params, orderStr, reverse, ...
        nPermVar, alphaVarPoint, alphaVarClust);

    %% -----------------------------------------------------
    % Part 6: TDR decoding
    %% -----------------------------------------------------
    fprintf('[%s] Decoding...\n', stimname);

    [wH_tdr, wC_tdr] = manteMaxNormAxes(B_heading, B_choice, orderStr);
    wH_tdr = normalizeAxisSafe_local(wH_tdr);
    wC_tdr = normalizeAxisSafe_local(wC_tdr);

    all_heads = [];
    for n=1:size(REG_FUNangle,1)
        if ~isempty(REG_FUNangle{n,ss})
            all_heads = [all_heads; REG_FUNangle{n,ss}(:)]; %#ok<AGROW>
        end
    end
    all_heads = all_heads(isfinite(all_heads));
    headLevels = unique(all_heads)';
    headLevels = sort(headLevels);

    cvDec = crossValidatedDecoding_TDR_signed_local( ...
        Test_REG_r, REG_FUNangle, REG_choice, ...
        wH_tdr, wC_tdr, ...
        ss, choiceMap, headLevels, ...
        trainFrac, nIterReal, nPermDec, alphaCluster);

    acc_head   = cvDec.acc_head;
    acc_choice = cvDec.acc_choice;

    sig_head   = cvDec.sig_head;
    sig_choice = cvDec.sig_choice;

    chance_head = 1 / numel(headLevels);
    acc_head_norm = (acc_head - chance_head) ./ (1 - chance_head);

    %% -----------------------------------------------------
    % Plot 1: Decoding
    %% -----------------------------------------------------
    f1 = figure('Color','w','Name',[areaName '_' durName '_' stimname '_HC_Decoding']);
    vecY = @(y, mask) y * ones(sum(mask), 1);

    subplot(2,1,1); hold on; box on
    plot(REG_t, acc_head_norm, 'k', 'LineWidth', 2);
    yline(0, 'k--');
    ylim([0 1]); yl=ylim;
    if any(sig_head)
        plot(REG_t(sig_head), vecY(yl(1)+0.03, sig_head), 'k.', 'MarkerSize', 10);
    end
    title('Heading decoding');
    ylabel('Accuracy');

    subplot(2,1,2); hold on; box on
    plot(REG_t, acc_choice, 'k', 'LineWidth', 2);
    yline(0.5, 'k--');
    ylim([0.3 1]); yl=ylim;
    if any(sig_choice)
        plot(REG_t(sig_choice), vecY(yl(1)+0.03, sig_choice), 'k.', 'MarkerSize', 10);
    end
    title('Choice decoding');
    xlabel('Time [ms]');
    ylabel('Accuracy');

    %% -----------------------------------------------------
    % Plot 2: Variance
    %% -----------------------------------------------------
    f2 = figure('Color','w','Name',[areaName '_' durName '_' stimname '_HC_Variance']);
    vnList = axisOrder;

    for k=1:2
        vn = vnList{k};
        subplot(2,1,k); hold on; box on
        SV = smooth(varMap.(vn),5);
        plot(REG_t, SV, 'k', 'LineWidth', 2);

        yl = ylim;
        m = sigVarMap.(vn);
        if any(m)
            plot(REG_t(m), (yl(1)+0.04*diff(yl))*ones(sum(m),1), 'k.', 'MarkerSize',10);
        end
        title([upper(vn) ' variance']);
        ylabel('Variance explained');
    end
    xlabel('Time [ms]');

    %% -----------------------------------------------------
    % Save
    %% -----------------------------------------------------
    Out = struct();
    Out.area      = areaName;
    Out.duration  = durName;
    Out.time      = REG_t;
    Out.stimname  = stimname;
    Out.axisOrder = axisOrder;

    Out.TDR.projection = projection;
    Out.TDR.Summary    = Summary;
    Out.TDR.varMap     = varMap;
    Out.TDR.varSig     = sigVarMap;

    Out.decoding.heading.acc         = acc_head;
    Out.decoding.heading.acc_norm    = acc_head_norm;
    Out.decoding.heading.sig         = sig_head;
    Out.decoding.heading.levels      = headLevels;
    Out.decoding.heading.chance      = chance_head;

    Out.decoding.choice.acc          = acc_choice;
    Out.decoding.choice.sig          = sig_choice;

    param_str = strjoin(axisOrder, '_');
    SaveFileName = fullfile(OutputPath, sprintf('%s_%s_%s_TDR_%s.mat', ...
        areaName, durName, stimname, param_str));
    save(SaveFileName, 'Out');

    saveas(f1, fullfile(OutputPath, sprintf('%s_%s_%s_HC_Decoding.png', areaName, durName, stimname)));
    saveas(f2, fullfile(OutputPath, sprintf('%s_%s_%s_HC_Variance.png', areaName, durName, stimname)));

    close(f1);
    close(f2);

    fprintf('Data saved to: %s\n', SaveFileName);

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mp = detectChoiceMapping(REG_choice, ss)
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

function [wH, wC] = manteMaxNormAxes(BH, BC, orderStr)
eps0=1e-12;
[~,tH] = max(sqrt(sum(BH.^2,1,'omitnan')));
wH = BH(:,tH);
wH = wH/(norm(wH)+eps0);

[~,tC] = max(sqrt(sum(BC.^2,1,'omitnan')));
wC = BC(:,tC);
wC = wC/(norm(wC)+eps0);

if strcmpi(orderStr,'HC')
    wC = wC - wH*(wH'*wC);
    wC = wC/(norm(wC)+eps0);
else
    wH = wH - wC*(wC'*wH);
    wH = wH/(norm(wH)+eps0);
end
end

function sigTDR = varianceSig_labelShuffle_axes_clusterperm_TDRonly( ...
    Test_REG_r, REG_FUNangle, REG_choice, ss, choiceMap, ...
    Test_PCA_Data, PCA_params_real, orderStr, reverse, ...
    nPerm, alphaPoint, alphaClust)

[~, SumReal] = TDR_time(Test_PCA_Data, [], PCA_params_real);
realVarTDR = SumReal.varCaptEachTargetedPCT(:,1:2);
T = size(realVarTDR,1);

nullT1 = nan(T,nPerm);
nullT2 = nan(T,nPerm);

N = numel(Test_REG_r);

for p=1:nPerm
    Hsh = cell(N,1);
    Csh = cell(N,1);
    for n=1:N
        h = REG_FUNangle{n,ss};
        c = REG_choice{n,ss};
        Y = Test_REG_r{n};
        if isempty(h) || isempty(c) || isempty(Y)
            Hsh{n} = h; Csh{n} = c; continue;
        end
        L = min([numel(h), numel(c), size(Y,1)]);
        rp = randperm(L);
        Hsh{n} = h(1:L); Hsh{n} = Hsh{n}(rp);
        Csh{n} = c(1:L); Csh{n} = Csh{n}(rp);
    end

    Bnull = tdrBetas_fromLabels(Test_REG_r, Hsh, Csh, choiceMap, orderStr, reverse, T);
    PCAp = PCA_params_real;
    PCAp.B = Bnull;

    [~, SumNull] = TDR_time(Test_PCA_Data, [], PCAp);
    vNullT = SumNull.varCaptEachTargetedPCT(:,1:2);
    nullT1(:,p) = vNullT(:,1);
    nullT2(:,p) = vNullT(:,2);
end

[sigT1, ~, ~, ~] = clusterPerm_1D(realVarTDR(:,1), nullT1, alphaPoint, alphaClust);
[sigT2, ~, ~, ~] = clusterPerm_1D(realVarTDR(:,2), nullT2, alphaPoint, alphaClust);

sigTDR = struct();
if strcmpi(orderStr,'HC')
    sigTDR.heading = sigT1;
    sigTDR.choice  = sigT2;
else
    sigTDR.choice  = sigT1;
    sigTDR.heading = sigT2;
end
end

function Bcell = tdrBetas_fromLabels(Test_REG_r, Hcell, Ccell, choiceMap, orderStr, reverse, T)

N = numel(Test_REG_r);
B_h = nan(N,T);
B_c = nan(N,T);

for n=1:N
    Y = Test_REG_r{n};
    if isempty(Y), continue; end

    h = Hcell{n};
    c = Ccell{n};
    if isempty(h) || isempty(c), continue; end

    L = min([size(Y,1), numel(h), numel(c)]);
    Y = Y(1:L,1:T);
    h = h(1:L);
    c_pm = choiceMap.toPM(c(1:L));

    if ~reverse
        X = [h(:), c_pm(:), ones(L,1)];
        colH = 1; colC = 2;
    else
        X = [c_pm(:), h(:), ones(L,1)];
        colC = 1; colH = 2;
    end

    if rcond(X'*X) < 1e-12, continue; end
    B = X \ Y;

    B_h(n,:) = B(colH,:);
    B_c(n,:) = B(colC,:);
end

if strcmpi(orderStr,'HC')
    Bcell = {B_h, B_c};
else
    Bcell = {B_c, B_h};
end
end

function out = crossValidatedDecoding_TDR_signed_local( ...
    Test_REG_r, REG_FUNangle, REG_choice, ...
    w_head, w_choice, ...
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

accHead_all   = nan(T, nIterReal);
accChoice_all = nan(T, nIterReal);

for it = 1:nIterReal
    splitInfo = buildTrainTestSplitPerNeuron_noStim_local( ...
        Test_REG_r, REG_FUNangle, REG_choice, ss, trainFrac);

    [Hp_tr, Hp_te, Cp_tr, Cp_te] = buildPoolsFromSplit_noStim_local( ...
        Test_REG_r, REG_FUNangle, REG_choice, ss, choiceMap, headLevels, splitInfo);

    tmp = decode1D_template_pseudopop_robust_local(Hp_tr, Hp_te, w_head, 1);
    if isempty(tmp), tmp = nan(T,1); end
    accHead_all(:,it) = tmp;

    tmp = decode1D_template_pseudopop_robust_local(Cp_tr, Cp_te, w_choice, 1);
    if isempty(tmp), tmp = nan(T,1); end
    accChoice_all(:,it) = tmp;
end

acc_head   = mean(accHead_all,2,'omitnan');
acc_choice = mean(accChoice_all,2,'omitnan');

nullHead   = nan(T,nPermDec);
nullChoice = nan(T,nPermDec);

for p = 1:nPermDec
    splitInfo = buildTrainTestSplitPerNeuron_noStim_local( ...
        Test_REG_r, REG_FUNangle, REG_choice, ss, trainFrac);

    [Hp_tr, Hp_te, Cp_tr, Cp_te] = buildPoolsFromSplit_noStim_local( ...
        Test_REG_r, REG_FUNangle, REG_choice, ss, choiceMap, headLevels, splitInfo);

    [Hp_tr_null, Hp_te_null] = shuffleHeadingPools_local(Hp_tr, Hp_te);
    [Cp_tr_null, Cp_te_null] = shuffleChoicePools_local(Cp_tr, Cp_te);

    tmp = decode1D_template_pseudopop_robust_local(Hp_tr_null, Hp_te_null, w_head, 1);
    if isempty(tmp), tmp = nan(T,1); end
    nullHead(:,p) = tmp;

    tmp = decode1D_template_pseudopop_robust_local(Cp_tr_null, Cp_te_null, w_choice, 1);
    if isempty(tmp), tmp = nan(T,1); end
    nullChoice(:,p) = tmp;
end

sig_head   = clusterPermSig_fromNull(acc_head, nullHead, alphaCluster);
sig_choice = clusterPermSig_fromChance_binary(acc_choice, 0.5, alphaCluster);

out = struct();
out.acc_head   = acc_head;
out.acc_choice = acc_choice;
out.sig_head   = sig_head;
out.sig_choice = sig_choice;
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

function sigMask = clusterPermSig_fromChance_binary(accReal, chanceLevel, alpha)
thr = chanceLevel;

T = numel(accReal);
nPerm = 1000;
maxMass = zeros(nPerm,1);

x = accReal - chanceLevel;
x(~isfinite(x)) = 0;

for p = 1:nPerm
    sgn = sign(randn(T,1));
    xp = x .* sgn;
    maxMass(p) = maxClusterMass(xp + chanceLevel, thr);
end

clThr = quantile(maxMass, 1-alpha);
sigMask = clusterMask_real(accReal, thr, clThr);
end

function [sigMask, thrPoint, clThr, maxMass] = clusterPerm_1D(realCurve, nullCurves, alphaPoint, alphaClust)
thrPoint = quantile(nullCurves(:), 1 - alphaPoint);

nPerm = size(nullCurves,2);
maxMass = zeros(nPerm,1);
for p=1:nPerm
    maxMass(p) = maxClusterMass(nullCurves(:,p), thrPoint);
end
clThr = quantile(maxMass, 1 - alphaClust);

sigMask = clusterMask_real(realCurve, thrPoint, clThr);
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