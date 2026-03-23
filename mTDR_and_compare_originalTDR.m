function TDR_mTDR_jointLR_seqPCA_and_compare()

clear; clc;

%% =========================================================
% unified paths
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

% ---- unified order (follow your previous batch style) ----
reverse   = 0;
orderStr  = 'HC';

% ---- practical joint low-rank settings ----
rankVec = [3 3];
maxIter = 100;
fitTol  = 1e-6;

% ---- display ----
wsProj    = 5;
FDR_q     = 0.01;
Np_pseudo = 1000;
magMode   = 'l2';

% ---- time ----
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

    run_one_dataset_jointLR(matfile, ds.area, ds.dur, ...
        window, stim_type, reverse, orderStr, ...
        rankVec, maxIter, fitTol, ...
        wsProj, FDR_q, Np_pseudo, magMode, span);

end

fprintf('\nAll datasets finished.\n');

end


%% ========================================================================
function run_one_dataset_jointLR(matfile, areaName, durName, ...
    window, stim_type, reverse, orderStr, ...
    rankVec, maxIter, fitTol, ...
    wsProj, FDR_q, Np_pseudo, magMode, span)

load(matfile);

%% --------------------------------------------------------
% dataset-specific timing / crop / heading values
%% --------------------------------------------------------
switch durName
    case '2s'
        t1 = 101:2000;
        crop_idx = 101:2000;
        headingValues = [-9 -3.6 -1.44 -0.58 0 0.58 1.44 3.6 9];

    case '1s'
        t1 = 101:1000;

        if strcmp(areaName,'VIP') || strcmp(areaName,'MSTd')
            crop_idx = 1186:2085;
            headingValues = [-12 -6 -1.5 0 1.5 6 12];
        elseif strcmp(areaName,'PIVC')
            crop_idx = 101:1000;
            headingValues = [-16 -8 -4 -2 -1 1 2 4 8 16];
        else
            error('Unknown areaName for 1s dataset: %s', areaName);
        end

    otherwise
        error('Unknown durName: %s', durName);
end

REG_t = square_sum(t1,span)/span;
T     = numel(REG_t);

REG_params_stim = {'choice','FUNangle','ONES'};
REG_choice_dim  = find(strcmp(REG_params_stim,'choice'));
REG_head_dim    = find(strcmp(REG_params_stim,'FUNangle'));

for ss = 1:numel(stim_type)

    if ss==1, stimname='Vest'; else, stimname='Vis'; end
    fprintf('\n=========== %s | %s %s ===========\n', stimname, areaName, durName);

    Ncells = size(REG_raster_per_trial,1);
    nAngle = size(REG_raster_per_trial{1,ss},1);
    assert(numel(headingValues)==nAngle, 'headingValues size must match nAngle');

    %%%% ============================================================
    % Part 1: Build trial-wise SDF list per neuron
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

    %%%% ============================================================
    % Part 2: TDR regression betas
    if ~reverse
        axisOrder = {'heading','choice'};
    else
        axisOrder = {'choice','heading'};
    end

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
            X = [h(:), c_pm(:), ones(length(h),1)];
            colH = 1; colC = 2;
        else
            X = [c_pm(:), h(:), ones(length(h),1)];
            colC = 1; colH = 2;
        end

        Y = Test_REG_r{n};
        if isempty(Y) || size(Y,2) < T, continue; end
        if size(X,1) ~= size(Y,1), continue; end
        if rcond(X'*X) < 1e-12, continue; end

        B = X \ Y(:,1:T);
        REG_bT(n,:,REG_head_dim)   = B(colH,:);
        REG_bT(n,:,REG_choice_dim) = B(colC,:);
    end

    B_heading = squeeze(REG_bT(:,:,REG_head_dim));
    B_choice  = squeeze(REG_bT(:,:,REG_choice_dim));

    %%%% ============================================================
    % Part 3: Condition means for TDR_time projection
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

    Acond = nan(T,nAngle,Ncells);
    for c=1:nAngle
        Acond(:,c,:) = Test_PCA_Data(1,c).A;
    end
    meanAcrossCond = squeeze(mean(Acond,2,'omitnan'));

    %%%% ============================================================
    % build per-neuron per-heading trial SDF for pseudo-trials
    SDF_trials = buildSDFtrials(REG_raster_per_trial, ss, window, crop_idx, T);

    %%%% ============================================================
    % Part 4: TDR_time
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

    [projection, Summary] = TDR_time(Test_PCA_Data,[],PCA_params);

    [wH_tdr, wC_tdr] = manteMaxNormAxes(B_heading, B_choice, orderStr);

    %%%% ============================================================
    % Part 5: joint low-rank mTDR core
    fprintf('[%s] joint low-rank mTDR fit...\n', stimname);

    [Ytrial, Xdesign, trialHeading] = buildPopulationTensor_HC_practical( ...
        REG_raster_per_trial, REG_FUNangle, REG_choice, ss, window, T, choiceMap, crop_idx);

    optsFit = struct();
    optsFit.maxIter = maxIter;
    optsFit.tol     = fitTol;
    optsFit.verbose = true;

    model = fit_mTDR_joint_lowrank_practical(Ytrial, Xdesign, rankVec, optsFit);

    W_h = model.W{1};
    W_c = model.W{2};

    if strcmpi(orderStr,'HC')
        W_h_orth = orth(W_h);
        W_c_orth = orth(W_c - W_h_orth*(W_h_orth' * W_c));
    else
        W_c_orth = orth(W_c);
        W_h_orth = orth(W_h - W_c_orth*(W_c_orth' * W_h));
    end

    Kh = min(3, size(W_h_orth,2));
    Kc = min(3, size(W_c_orth,2));

    if Kh < 3
        warning('[%s] Heading mTDR rank < 3 (Kh=%d). 3D plots will degrade.', stimname, Kh);
    end
    if Kc < 3
        warning('[%s] Choice mTDR rank < 3 (Kc=%d). 3D plots will degrade.', stimname, Kc);
    end

    %%%% ============================================================
    % Part 6: mTDR condition-mean projections
    proj_mTDR = struct();
    proj_mTDR.heading = nan(T, nAngle, Kh);
    proj_mTDR.choice  = nan(T, nAngle, Kc);

    for c = 1:nAngle
        Ams = Test_PCA_Data(1,c).A - meanAcrossCond;
        proj_mTDR.heading(:,c,:) = Ams * W_h_orth(:,1:Kh);
        proj_mTDR.choice(:,c,:)  = Ams * W_c_orth(:,1:Kc);
    end

    %%%% ============================================================
    % Part 7: seqPCA
    R_seq = struct(); tE = struct(); tM = struct();
    proj_seq = struct();

    if Kh >= 3
        [R_seq.heading, tE.heading, tM.heading] = seqPCA_3axes(proj_mTDR.heading(:,:,1:3));
        proj_seq.heading = applyRotation3(proj_mTDR.heading(:,:,1:3), R_seq.heading);
    else
        R_seq.heading = eye(Kh); tE.heading = NaN; tM.heading = NaN;
        proj_seq.heading = proj_mTDR.heading;
    end

    if Kc >= 3
        [R_seq.choice, tE.choice, tM.choice] = seqPCA_3axes(proj_mTDR.choice(:,:,1:3));
        proj_seq.choice = applyRotation3(proj_mTDR.choice(:,:,1:3), R_seq.choice);
    else
        R_seq.choice = eye(Kc); tE.choice = NaN; tM.choice = NaN;
        proj_seq.choice = proj_mTDR.choice;
    end

    Useq = struct();
    if Kh >= 3
        Useq.heading = W_h_orth(:,1:3) * R_seq.heading;
    else
        Useq.heading = W_h_orth(:,1:Kh) * R_seq.heading;
    end
    if Kc >= 3
        Useq.choice  = W_c_orth(:,1:3) * R_seq.choice;
    else
        Useq.choice  = W_c_orth(:,1:Kc) * R_seq.choice;
    end

    w1D = struct();
    w1D.heading = wH_tdr;
    w1D.choice  = wC_tdr;

    %%%% ============================================================
    % Part 8-9: plots
    for kVar = 1:2
        v = axisOrder{kVar};
        labelStr = upper(v);
        plot_mTDR_3D_traj(REG_t, proj_seq.(v), stimname, labelStr, tE.(v), tM.(v), ss, headingValues);
        plot_seq_axis_curves(REG_t, proj_seq.(v), stimname, labelStr, wsProj, ss, headingValues);
    end

    %%%% ============================================================
    % Part 10: pseudo-trial encoding strength
    fprintf('[%s] Pseudo-trial encoding strength (Np=%d)...\n', stimname, Np_pseudo);

    medMag = struct();
    medMag_signed = struct();

    for kVar = 1:2
        v = axisOrder{kVar};

        [m1D, mm3D, m1Ds, mm3Ds] = pseudoTrialMedianStrength( ...
            SDF_trials, meanAcrossCond, w1D.(v), Useq.(v), Np_pseudo, magMode, headingValues);

        medMag.(v).TDR         = m1D;
        medMag.(v).mTDR        = mm3D;
        medMag_signed.(v).TDR  = m1Ds;
        medMag_signed.(v).mTDR = mm3Ds;
    end

    sig = struct(); pvals = struct();
    for kVar = 1:2
        v = axisOrder{kVar};
        [sig.(v), pvals.(v)] = signrank_FDR_time(medMag.(v).mTDR, medMag.(v).TDR, FDR_q);
    end

    for kVar = 1:2
        v = axisOrder{kVar};
        plot_TDR_mTDR_compare_strongest(REG_t, ...
            medMag_signed.(v).TDR, medMag_signed.(v).mTDR, sig.(v), ...
            stimname, upper(v), headingValues, ss);
    end

    %%%% ============================================================
    % SAVE
    Out = struct();
    Out.time = REG_t;
    Out.stimname = stimname;
    Out.orderStr = orderStr;
    Out.reverse  = reverse;
    Out.axisOrder = axisOrder;

    Out.TDR.projection = projection;
    Out.TDR.Summary    = Summary;

    Out.mTDR.model = model;
    Out.mTDR.U.heading = W_h_orth(:,1:Kh);
    Out.mTDR.U.choice  = W_c_orth(:,1:Kc);

    Out.mTDR.seq.heading.R = R_seq.heading; Out.mTDR.seq.heading.boundary = [tE.heading tM.heading];
    Out.mTDR.seq.choice.R  = R_seq.choice;  Out.mTDR.seq.choice.boundary  = [tE.choice  tM.choice];

    Out.mTDR.seqProj.heading_3D = proj_seq.heading;
    Out.mTDR.seqProj.choice_3D  = proj_seq.choice;

    if size(proj_seq.heading,3) >= 3
        Out.mTDR.seqProj.heading_early  = squeeze(proj_seq.heading(:,:,1));
        Out.mTDR.seqProj.heading_middle = squeeze(proj_seq.heading(:,:,2));
        Out.mTDR.seqProj.heading_late   = squeeze(proj_seq.heading(:,:,3));
    end
    if size(proj_seq.choice,3) >= 3
        Out.mTDR.seqProj.choice_early   = squeeze(proj_seq.choice(:,:,1));
        Out.mTDR.seqProj.choice_middle  = squeeze(proj_seq.choice(:,:,2));
        Out.mTDR.seqProj.choice_late    = squeeze(proj_seq.choice(:,:,3));
    end

    Out.compare.heading.mag_TDR  = medMag_signed.heading.TDR;
    Out.compare.heading.mag_mTDR = medMag_signed.heading.mTDR;
    Out.compare.heading.p   = pvals.heading;
    Out.compare.heading.sig = sig.heading;

    Out.compare.choice.mag_TDR  = medMag_signed.choice.TDR;
    Out.compare.choice.mag_mTDR = medMag_signed.choice.mTDR;
    Out.compare.choice.p   = pvals.choice;
    Out.compare.choice.sig = sig.choice;

    Out.compare.byOrder{1} = Out.compare.(axisOrder{1});
    Out.compare.byOrder{2} = Out.compare.(axisOrder{2});

    OutputPath = 'Z:\Users\TDR\results\mTDR\';
    if ~exist(OutputPath,'dir')
        mkdir(OutputPath);
    end

    param_str = strjoin(axisOrder, '_');
    save_name_str = sprintf('%s_%s_%s_mTDR_jointLR_seqPCA_compare_%s.mat', ...
        areaName, durName, stimname, orderStr);
    SaveFileName = fullfile(OutputPath, save_name_str);

    save(SaveFileName, 'Out');
    fprintf('Data saved to: %s\n', SaveFileName);

    close all
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------- NEW CORE HELPERS -----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ytrial, Xdesign, trialHeading] = buildPopulationTensor_HC_practical( ...
    REG_raster_per_trial, REG_FUNangle, REG_choice, ss, window, T, choiceMap, crop_idx)
% Build trial-wise pseudo-population tensor:
%   Ytrial  : N x T x K
%   Xdesign : K x 2  [heading_z, choice_z]
%   trialHeading : K x 1 original heading value
%
% Practical alignment:
%   For each neuron, flatten trials in the same order as the raster traversal.
%   Use the minimum available trial count across neurons.

Ncells = size(REG_raster_per_trial,1);

cellY = cell(Ncells,1);
cellX = cell(Ncells,1);
cellH = cell(Ncells,1);

for n = 1:Ncells
    ras = REG_raster_per_trial{n,ss};
    h   = REG_FUNangle{n,ss};
    c   = REG_choice{n,ss};

    if isempty(ras) || isempty(h) || isempty(c)
        continue;
    end

    h = h(:);
    c = c(:);
    c_pm = choiceMap.toPM(c);

    Lmeta = min(numel(h), numel(c_pm));

    % zscore within neuron's trial list (practical version)
    hz = h(1:Lmeta);
    cz = c_pm(1:Lmeta);

    hz = (hz - mean(hz,'omitnan')) ./ max(std(hz,[],'omitnan'), eps);
    cz = (cz - mean(cz,'omitnan')) ./ max(std(cz,[],'omitnan'), eps);

    Ylist = [];
    Xlist = [];
    Hlist = [];

    idx = 0;
    nRep = size(ras,2);
    for r = 1:nRep
        blk = double(squeeze(ras(:,r,:))); % nAngle x samples
        for a = 1:size(blk,1)
            idx = idx + 1;
            if idx > Lmeta, break; end

            cs = conv(blk(a,:), window, 'same');
            cs = cs(crop_idx);
            y  = downsample(cs,10);
            y  = y(1:T);

            Ylist(end+1,:) = y; %#ok<AGROW>
            Xlist(end+1,:) = [hz(idx), cz(idx)]; %#ok<AGROW>
            Hlist(end+1,1) = h(idx); %#ok<AGROW>
        end
    end

    if ~isempty(Ylist)
        cellY{n} = Ylist;
        cellX{n} = Xlist;
        cellH{n} = Hlist;
    end
end

% minimum trial count across valid neurons
valid = ~cellfun(@isempty, cellY);
assert(any(valid), 'No valid neurons/trials found for buildPopulationTensor_HC_practical.');

Kmin = min(cellfun(@(x) size(x,1), cellY(valid)));
N = Ncells;

Ytrial = nan(N, T, Kmin);
Xdesign = nan(Kmin, 2);
trialHeading = nan(Kmin,1);

for k = 1:Kmin
    for n = 1:N
        if isempty(cellY{n}), continue; end
        Ytrial(n,:,k) = cellY{n}(k,:);
    end

    % take regressor values from first valid neuron
    firstValid = find(valid,1,'first');
    Xdesign(k,:) = cellX{firstValid}(k,:);
    trialHeading(k) = cellH{firstValid}(k);
end

% replace any remaining NaN neurons with 0 after mean-centering assumption
Ytrial(~isfinite(Ytrial)) = 0;
end

function model = fit_mTDR_joint_lowrank_practical(Y, X, rankVec, opts)
% Practical joint low-rank regression:
%   Y(:,:,k) ~ sum_p X(k,p) * W{p} * S{p}'
%
% Y: N x T x K
% X: K x P
%
% This is a practical ALS fit of the Aoi-style core structure.

[N,T,K] = size(Y);
[K2,P] = size(X);
assert(K==K2, 'Trial count mismatch.');

maxIter = get_opt_local(opts, 'maxIter', 100);
tol     = get_opt_local(opts, 'tol', 1e-6);
verbose = get_opt_local(opts, 'verbose', false);

W = cell(P,1);
S = cell(P,1);

% ----- initialize each variable from crude full-rank regression matrix
for p = 1:P
    xp = X(:,p);
    denom = sum(xp.^2) + eps;

    Bp0 = zeros(N,T);
    for k = 1:K
        Bp0 = Bp0 + xp(k) * Y(:,:,k);
    end
    Bp0 = Bp0 / denom;

    [U,Sv,V] = svd(Bp0, 'econ');
    rp = min([rankVec(p), size(U,2), size(V,2)]);
    W{p} = U(:,1:rp) * Sv(1:rp,1:rp);
    S{p} = V(:,1:rp);
end

prevLoss = inf;

for iter = 1:maxIter

    % ----- update W
    for p = 1:P
        rp = size(S{p},2);

        A = zeros(rp,rp);
        B = zeros(N,rp);

        Gs = S{p}' * S{p};

        for k = 1:K
            xkp = X(k,p);
            Yres = Y(:,:,k);

            for q = 1:P
                if q==p, continue; end
                Yres = Yres - X(k,q) * (W{q} * S{q}');
            end

            A = A + (xkp^2) * Gs;
            B = B + xkp * (Yres * S{p});
        end

        if rcond(A) > 1e-12
            W{p} = B / A;
        end
    end

    % ----- update S
    for p = 1:P
        rp = size(W{p},2);

        A = zeros(rp,rp);
        B = zeros(T,rp);

        Gw = W{p}' * W{p};

        for k = 1:K
            xkp = X(k,p);
            Yres = Y(:,:,k);

            for q = 1:P
                if q==p, continue; end
                Yres = Yres - X(k,q) * (W{q} * S{q}');
            end

            A = A + (xkp^2) * Gw;
            B = B + xkp * (Yres' * W{p});
        end

        if rcond(A) > 1e-12
            S{p} = B / A;   % T x rp
        end
    end

    % ----- compute loss
    loss = 0;
    for k = 1:K
        Yhat = zeros(N,T);
        for p = 1:P
            Yhat = Yhat + X(k,p) * (W{p} * S{p}');
        end
        E = Y(:,:,k) - Yhat;
        loss = loss + sum(E(:).^2);
    end

    if verbose
        fprintf('  iter %d, loss = %.6g\n', iter, loss);
    end

    if isfinite(prevLoss)
        relchg = abs(prevLoss - loss) / max(prevLoss, eps);
        if relchg < tol
            break;
        end
    end
    prevLoss = loss;
end

model = struct();
model.W = W;
model.S = S;
model.B = cell(P,1);
for p = 1:P
    model.B{p} = W{p} * S{p}';   % N x T
end
model.rankVec = rankVec;
model.loss = prevLoss;
model.iter = iter;
end

function val = get_opt_local(opts, name, defaultVal)
if isfield(opts, name)
    val = opts.(name);
else
    val = defaultVal;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------- ORIGINAL HELPERS -----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    mp.toPM = @(c) (c==2) - (c==1);   % 1->-1, 2->+1
else
    error('Unknown choice labels. Unique values: %s', mat2str(u'));
end
end

function [wH, wC] = manteMaxNormAxes(BH, BC, orderStr)
eps0=1e-12;
[~,tH] = max(sqrt(sum(BH.^2,1,'omitnan'))); wH = BH(:,tH); wH = wH/(norm(wH)+eps0);
[~,tC] = max(sqrt(sum(BC.^2,1,'omitnan'))); wC = BC(:,tC); wC = wC/(norm(wC)+eps0);
if strcmpi(orderStr,'HC')
    wC = wC - wH*(wH'*wC); wC = wC/(norm(wC)+eps0);
else
    wH = wH - wC*(wC'*wH); wH = wH/(norm(wH)+eps0);
end
end

function [R, tE, tM] = seqPCA_3axes(X)
[T,~,D] = size(X);
assert(D==3, 'seqPCA_3axes expects 3D input.');

Xc = X;
for t=1:T
    mu = mean(squeeze(X(t,:,:)),1,'omitnan');
    Xc(t,:,:) = squeeze(X(t,:,:)) - mu;
end

[tE, a1] = seq_axis1_from_prefix(Xc, 1, T);
Xr = remove_axis_component(Xc, a1);

if tE < T-2
    [tM, a2] = seq_axis1_from_prefix(Xr, tE+1, T);
else
    tM = T;
    a2 = randn(3,1); a2 = a2 / norm(a2);
end

a1 = a1(:)/norm(a1);
a2 = a2(:);
a2 = a2 - a1*(a1'*a2);
if norm(a2) < 1e-12
    a2 = null(a1'); a2 = a2(:,1);
else
    a2 = a2 / norm(a2);
end
a3 = cross(a1, a2);
if norm(a3) < 1e-12
    tmp = null([a1 a2]'); a3 = tmp(:,1);
end
a3 = a3 / norm(a3);

R = [a1 a2 a3];
end

function [tauBest, a1] = seq_axis1_from_prefix(X, tStart, tEnd)
tauBest = tStart;
bestScore = -Inf;
a1 = [1;0;0];

for tau = tStart:tEnd
    Y = reshape(X(tStart:tau,:,:), [], 3);
    Y = Y(all(isfinite(Y),2),:);
    if size(Y,1) < 10, continue; end

    C = cov(Y,1);
    [V,D] = eig(C);
    [lam,idx] = sort(diag(D),'descend');
    v1 = V(:,idx(1));
    score = lam(1) / max(eps, sum(lam));

    if score > bestScore
        bestScore = score;
        tauBest = tau;
        a1 = v1;
    end
end

a1 = a1(:);
if norm(a1) < 1e-12, a1 = [1;0;0]; end
a1 = a1 / norm(a1);
end

function Xr = remove_axis_component(X, a)
a = a(:)/norm(a);
[T,C,~] = size(X);
Xr = X;
for t=1:T
    M = squeeze(X(t,:,:));
    proj = (M*a)*a';
    Xr(t,:,:) = M - proj;
end
end

function Xrot = applyRotation3(X, R)
[T,C,~] = size(X);
Xrot = nan(T,C,3);
for t=1:T
    M = squeeze(X(t,:,:));
    Xrot(t,:,:) = M * R;
end
end

function plot_mTDR_3D_traj(t, projSeq, stimname, labelStr, tE, tM, ss, headingValues)
if size(projSeq,3) < 3
    figure('Color','w','Name',[stimname ' mTDR traj ' labelStr ' (2D)']);
    hold on; box on; grid on;
    xlabel([labelStr ' axis 1']); ylabel([labelStr ' axis 2']);
    title([stimname ' mTDR ' labelStr ' subspace (rank<3)']);
    C = size(projSeq,2);
    for c=1:C
        X = squeeze(projSeq(:,c,1));
        Y = squeeze(projSeq(:,c,min(2,size(projSeq,3))));
        plot(X,Y,'LineWidth',1.2);
    end
    return
end

figure('Color','w','Name',[stimname ' mTDR traj ' labelStr ' (3D)']);
hold on; box on; grid on;
xlabel([labelStr ' early']); ylabel([labelStr ' middle']); zlabel([labelStr ' late']);
title([stimname ' mTDR ' labelStr ' trajectories in seqPCA coords']);

C = size(projSeq,2);
cols = headingGradientColors(headingValues, ss);

for c=1:C
    X = squeeze(projSeq(:,c,1));
    Y = squeeze(projSeq(:,c,2));
    Z = squeeze(projSeq(:,c,3));
    plot3(X,Y,Z,'LineWidth',1.4,'Color',cols(c,:));
    plot3(X(1),Y(1),Z(1),'o','MarkerSize',6,'MarkerFaceColor',cols(c,:),'MarkerEdgeColor','none');
    plot3(X(end),Y(end),Z(end),'s','MarkerSize',6,'MarkerFaceColor',cols(c,:),'MarkerEdgeColor','none');
end
view(3);

Xall = projSeq(:,:,1); Xall = Xall(:);
Yall = projSeq(:,:,2); Yall = Yall(:);
Zall = projSeq(:,:,3); Zall = Zall(:);

Lx = niceSymLimit(max(abs(Xall(isfinite(Xall)))));
Ly = niceSymLimit(max(abs(Yall(isfinite(Yall)))));
Lz = niceSymLimit(max(abs(Zall(isfinite(Zall)))));

xlim([-Lx Lx]); ylim([-Ly Ly]); zlim([-Lz Lz]);
xticks([-Lx 0 Lx]); yticks([-Ly 0 Ly]); zticks([-Lz 0 Lz]);

if isfinite(tE) && isfinite(tM)
    txt = sprintf('seqPCA boundaries: early<=%g ms, middle<=%g ms', t(tE), t(tM));
    text(0,0,0,txt,'FontSize',10,'Color',[.2 .2 .2]);
end
end

function plot_seq_axis_curves(t, projSeq, stimname, labelStr, ws, ss, headingValues)
[~, nAngle, K] = size(projSeq);

if K >= 3
    axNames = {'early','middle','late'};
else
    axNames = arrayfun(@(k)sprintf('axis%d',k), 1:K, 'UniformOutput', false);
end

cols = headingGradientColors(headingValues, ss);

figure('Color','w','Name',[stimname ' mTDR axis curves ' labelStr]);
tl = tiledlayout(K,1,'Padding','compact','TileSpacing','compact');
title(tl, [stimname ' mTDR ' labelStr ' axis time courses (each heading)']);

for k = 1:K
    nexttile; hold on; box on
    for c = 1:nAngle
        y = squeeze(projSeq(:,c,k));
        if ws > 1
            y = smoothdata(y, 'movmean', ws);
        end
        plot(t, y, 'LineWidth', 1.3, 'Color', cols(c,:));
    end
    yline(0,'Color',[.7 .7 .7]);
    xlabel('Time (ms)'); ylabel('Projection');
    title([upper(labelStr) ' - ' axNames{k} ' axis']);
end

lg = legend(arrayfun(@(v)sprintf('%g°',v), headingValues, 'UniformOutput', false), ...
    'Location','eastoutside');
lg.Title.String = 'Heading';
end

function [sigMask, pvals] = signrank_FDR_time(X_mTDR, X_TDR, q)
[T,~] = size(X_mTDR);
pvals = nan(T,1);
for tt=1:T
    x = X_mTDR(tt,:)';
    y = X_TDR(tt,:)';
    ok = isfinite(x) & isfinite(y);
    if sum(ok) >= 3
        pvals(tt) = signrank(x(ok), y(ok), 'tail','right');
    end
end
sigMask = fdr_bh_mask(pvals, q);
end

function plot_TDR_mTDR_compare_strongest(t, magTDR, magMTDR, sigMask, stimname, labelStr, headingValues, ss)
[~,nAngle] = size(magTDR);
assert(numel(headingValues)==nAngle);

[~, iNeg] = min(headingValues);
[~, iPos] = max(headingValues);

cols = headingGradientColors(headingValues, ss);
cNeg = cols(iNeg,:);
cPos = cols(iPos,:);

figure('Color','w','Name',[stimname ' compare strongest ' labelStr]);
hold on; box on;

plot(t, magTDR(:,iNeg),  '--', 'LineWidth', 1.8, 'Color', cNeg);
plot(t, magMTDR(:,iNeg), '-',  'LineWidth', 2.4, 'Color', cNeg);

plot(t, magTDR(:,iPos),  '--', 'LineWidth', 1.8, 'Color', cPos);
plot(t, magMTDR(:,iPos), '-',  'LineWidth', 2.4, 'Color', cPos);

yl = ylim;
y0 = yl(1) + 0.03*diff(yl);
if any(sigMask)
    plot(t(sigMask), y0*ones(sum(sigMask),1), '|', 'Color', [0.6 0.6 0.6], 'MarkerSize', 6, 'LineWidth', 1.2);
end

title([stimname ' ' labelStr ': strongest headings shown; sig across all headings (FDR q=0.01)']);
xlabel('Time (ms)'); ylabel('Encoding magnitude');

legend({ ...
    sprintf('TDR %g°', headingValues(iNeg)), sprintf('mTDR %g°', headingValues(iNeg)), ...
    sprintf('TDR %g°', headingValues(iPos)), sprintf('mTDR %g°', headingValues(iPos)), ...
    'sig (all headings)'}, ...
    'Location','best');
end

function sig = fdr_bh_mask(p, q)
p = p(:);
sig = false(size(p));
ok = isfinite(p);
m = sum(ok);
if m==0, return; end
ps = sort(p(ok));
th = (1:m)'/m * q;
k = find(ps <= th, 1, 'last');
if isempty(k), return; end
pcrit = ps(k);
sig(ok) = p(ok) <= pcrit;
end

function cols = headingGradientColors(headingVals, ss)
headingVals = headingVals(:);
hmin = min(headingVals);
hmax = max(headingVals);
u = (headingVals - hmin) ./ (hmax - hmin + eps);

if ss==1
    c0 = [0,171,209]/255;
    c1 = [0,19,70]/255;
elseif ss==2
    c0 = [209,0,171]/255;
    c1 = [70,0,19]/255;
else
    error('ss must be 1 (Vest) or 2 (Vis)');
end

cols = (1-u).*c0 + u.*c1;
end

function SDF_trials = buildSDFtrials(REG_raster_per_trial, ss, window, crop_idx, T)
% Build per-neuron, per-heading, per-trial SDF after cropping and downsampling
%
% Output:
%   SDF_trials{n, c} = [nTrials x T]
%
% Inputs:
%   REG_raster_per_trial : cell array, neuron x modality
%   ss                   : modality index
%   window               : Gaussian smoothing kernel
%   crop_idx             : crop indices in original time base
%   T                    : expected number of downsampled time points

Ncells = size(REG_raster_per_trial, 1);

% infer number of heading conditions from first non-empty entry
nAngle = [];
for n = 1:Ncells
    tmp = REG_raster_per_trial{n, ss};
    if ~isempty(tmp)
        nAngle = size(tmp,1);
        break;
    end
end

if isempty(nAngle)
    error('buildSDFtrials: no non-empty REG_raster_per_trial found for ss=%d.', ss);
end

SDF_trials = cell(Ncells, nAngle);

for n = 1:Ncells
    tmp = REG_raster_per_trial{n, ss};

    if isempty(tmp)
        for c = 1:nAngle
            SDF_trials{n,c} = [];
        end
        continue;
    end

    for c = 1:size(tmp,1)
        spk = double(squeeze(tmp(c,:,:)));   % [nTrials x rawTime]
        if isempty(spk)
            SDF_trials{n,c} = [];
            continue;
        end

        nTr = size(spk,1);
        buf = nan(nTr, T);

        for k = 1:nTr
            cs = conv(spk(k,:), window, 'same');
            cs = cs(crop_idx);
            ds = downsample(cs, 10);

            % --- make length exactly T ---
            if numel(ds) > T
                ds = ds(1:T);
            elseif numel(ds) < T
                ds = [ds, nan(1, T-numel(ds))];
            end

            buf(k,:) = ds;
        end

        SDF_trials{n,c} = buf;
    end

    % if some neurons have fewer heading rows than nAngle
    for c = size(tmp,1)+1:nAngle
        SDF_trials{n,c} = [];
    end
end

end

function [med1D, medm3D, med1D_signed, medm3D_signed] = pseudoTrialMedianStrength( ...
    SDF_trials, meanAcrossCond, w1D, Useq, Np, magMode, headingValues)

w1D = w1D(:);
[T, N] = size(meanAcrossCond);
nAngle = size(SDF_trials,2);
K = size(Useq,2);

med1D  = nan(T, nAngle);
medm3D = nan(T, nAngle);
med1D_signed  = nan(T, nAngle);
medm3D_signed = nan(T, nAngle);

okW = isfinite(w1D);
Useq(~isfinite(Useq)) = 0;

for c = 1:nAngle
    mag1D_all  = nan(Np, T);
    magm3D_all = nan(Np, T);

    for p = 1:Np
        A = nan(T, N);

        for n = 1:N
            trials_n = SDF_trials{n,c};
            if isempty(trials_n), continue; end
            rr = randi(size(trials_n,1));
            A(:,n) = trials_n(rr,:)';
        end

        A = A - meanAcrossCond;

        s1 = A(:,okW) * w1D(okW);
        mag1D_all(p,:) = abs(s1)';

        P = A * Useq;
        if strcmpi(magMode,'l1')
            magm = sum(abs(P(:,1:min(3,K))),2);
        else
            magm = sqrt(sum(P(:,1:min(3,K)).^2,2));
        end
        magm3D_all(p,:) = magm';
    end

    m1  = median(mag1D_all,  1,'omitnan')';
    m3  = median(magm3D_all, 1,'omitnan')';

    med1D(:,c)  = m1;
    medm3D(:,c) = m3;

    sgn = sign(headingValues(c));
    if sgn==0, sgn = 1; end
    med1D_signed(:,c)  = sgn * m1;
    medm3D_signed(:,c) = sgn * m3;
end
end

function L = niceSymLimit(maxAbsVal)
if isempty(maxAbsVal) || ~isfinite(maxAbsVal) || maxAbsVal<=0
    L = 1;
    return
end
step = 10;
L = ceil(maxAbsVal/step)*step;
if L==0, L=step; end
end