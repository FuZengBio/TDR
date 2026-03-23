function HC_choiceAxis_Vest_vs_Vis_batch()
% Compute vestibular-vs-visual RAW choice-axis angle
% by fusing the previous three single-dataset versions
%
% Unified data path and result path.

clear; clc;

%% ========================================================================
% unified paths
%% ========================================================================
FilePath   = 'Z:\Users\TDR\Data\';

%% ========================================================================
% unified settings
%% ========================================================================
sigma = 30;
half_guas = sigma * 3;
dist = (-(half_guas):(half_guas))';
window = 1/(sigma*sqrt(2*pi)) * exp(-((dist.^2)/(2*sigma^2))) * 1000;

stim_type = [1, 2];   % 1 = vestibular, 2 = visual
span = 10;

REG_params_kernel = {'VELangle','ACCangle','choice'};
REG_vel_dimK    = find(strcmp(REG_params_kernel,'VELangle'));
REG_acc_dimK    = find(strcmp(REG_params_kernel,'ACCangle'));
REG_choice_dimK = find(strcmp(REG_params_kernel,'choice'));

%% ========================================================================
% run list (5 datasets)
%% ========================================================================
datasets = { ...
    struct('area','VIP',  'dur','1s', 'mat','VIP_Data_1s_REG.mat'), ...
    struct('area','VIP',  'dur','2s', 'mat','VIP_Data_2s_REG.mat'), ...
    struct('area','MSTd', 'dur','1s', 'mat','MSTd_Data_1s_REG.mat'), ...
    struct('area','MSTd', 'dur','2s', 'mat','MSTd_Data_2s_REG.mat'), ...
    struct('area','PIVC', 'dur','1s', 'mat','PIVC_Data_1s_REG.mat') ...
    };

%% ========================================================================
% run all
%% ========================================================================
for kk = 1:numel(datasets)

    ds = datasets{kk};

    fprintf('\n====================================================\n');
    fprintf('Running %d/%d: %s %s\n', kk, numel(datasets), ds.area, ds.dur);
    fprintf('====================================================\n');

    matfile = fullfile(FilePath, ds.mat);
    if ~exist(matfile, 'file')
        warning('File not found: %s. Skip.', matfile);
        continue;
    end

    run_one_dataset(matfile, ds.area, ds.dur,...
        window, stim_type, span, ...
        REG_params_kernel, REG_vel_dimK, REG_acc_dimK, REG_choice_dimK);

end

fprintf('\nAll datasets finished.\n');

end


%% ========================================================================
function run_one_dataset(matfile, areaName, durName,...
    window, stim_type, span, ...
    REG_params_kernel, REG_vel_dimK, REG_acc_dimK, REG_choice_dimK)

load(matfile);

%% ========================================================================
% dataset-specific timing / kernel
%% ========================================================================
switch durName
    case '2s'
        t1 = 101:2000;
        crop_idx = 101:2000;
        gaussArg = [0.13, 6, 2];

    case '1s'
        t1 = 101:1000;
        gaussArg = [0.13, 4, 1];

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

REG_t = square_sum(t1, span) / span;

%% ========================================================================
% build ideal kernels
%% ========================================================================
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

%% ========================================================================
% compute raw choice axis separately for vestibular and visual
%% ========================================================================
nCells = size(REG_raster_per_trial, 1);

ChoiceAxisRaw = cell(1,2);   % {1}=Vest, {2}=Vis
AllBetasRaw   = cell(1,2);

for ss = 1:length(stim_type)

    fprintf('\n---------------------------------------------\n');
    if ss == 1
        fprintf('Processing vestibular condition...\n');
    else
        fprintf('Processing visual condition...\n');
    end

    %% --------------------------------------------------------------------
    % Step 1: recompute SDF
    %% --------------------------------------------------------------------
    clear Test_REG_r
    Test_REG_r = cell(nCells,1);

    for nCell = 1:nCells
        temp_spike = REG_raster_per_trial{nCell, ss};

        if isempty(temp_spike)
            Test_REG_r{nCell,1} = [];
            continue;
        end

        sdf = [];

        for i = 1:size(temp_spike,2)
            spike = double(squeeze(temp_spike(:,i,:)));
            [ntrials, ~] = size(spike);

            for ii = 1:ntrials
                convspike = conv(spike(ii,:), window, 'same');
                convspike = convspike(crop_idx);
                tempdata = downsample(convspike, 10);
                sdf((i-1)*ntrials + ii, :) = tempdata; %#ok<AGROW>
            end
        end

        Test_REG_r{nCell,1} = sdf;
    end

    %% --------------------------------------------------------------------
    % Step 2: kernel regression per neuron
    %% --------------------------------------------------------------------
    REG_bVesK_Test = cell(length(REG_params_kernel),1);
    for b = 1:length(REG_params_kernel)
        REG_bVesK_Test{b} = nan(nCells,1);
    end

    for Ni = 1:nCells

        TMP_FVesK = [];
        TMP = [];

        tempdata = REG_Stim{Ni};

        if isempty(tempdata)
            continue;
        end

        TMP_Ves = (tempdata == 1);

        if sum(TMP_Ves) == 0
            continue;
        end

        for b = 1:length(REG_params_kernel)
            if ss == 1
                eval(sprintf('TMP(b,:) = REG_%s{Ni}(find(TMP_Ves));', REG_params_kernel{b}));
            else
                visGroup = Ni + size(REG_Stim,1);
                eval(sprintf('TMP(b,:) = REG_%s{visGroup}(find(TMP_Ves));', REG_params_kernel{b}));
            end
        end

        for trl = 1:sum(TMP_Ves)
            TMP_FVesK = [TMP_FVesK, repmat(TMP(:,trl),1,sum(REG_t<=max(t1))) .* kernel]; %#ok<AGROW>
        end

        temp_r = Test_REG_r{Ni};

        if isempty(temp_r)
            continue;
        end

        temp_r = temp_r(TMP_Ves, 1:sum(REG_t<=max(t1)));

        if isempty(temp_r)
            continue;
        end

        r = bsxfun(@minus, temp_r, mean(temp_r));

        X = TMP_FVesK * TMP_FVesK';
        y = TMP_FVesK * reshape(r', [], 1);

        if rank(X) < size(X,1)
            TMP_bK = pinv(X) * y;
        else
            TMP_bK = X \ y;
        end

        for b = 1:length(REG_params_kernel)
            REG_bVesK_Test{b}(Ni,1) = TMP_bK(b);
        end
    end

    AllBetasRaw{ss} = REG_bVesK_Test;
    ChoiceAxisRaw{ss} = REG_bVesK_Test{REG_choice_dimK}(:);

    fprintf('Done. Valid choice-beta neurons: %d / %d\n', ...
        sum(~isnan(ChoiceAxisRaw{ss})), nCells);
end

%% ========================================================================
% compare vestibular-choice axis vs visual-choice axis
%% ========================================================================
choice_vest = ChoiceAxisRaw{1};
choice_vis  = ChoiceAxisRaw{2};

valid_id = ~isnan(choice_vest) & ~isnan(choice_vis);

u = choice_vest(valid_id);
v = choice_vis(valid_id);

if isempty(u) || isempty(v)
    error('No shared valid neurons for comparing vestibular and visual choice axes.');
end

if norm(u) == 0 || norm(v) == 0
    error('One of the choice axes has zero norm. Cannot compute angle.');
end

u_unit = u / norm(u);
v_unit = v / norm(v);

cos_theta_raw = dot(u_unit, v_unit);
cos_theta_raw = max(min(cos_theta_raw, 1), -1);
angle_raw = acosd(cos_theta_raw);

cos_theta_unsigned = abs(dot(u_unit, v_unit));
cos_theta_unsigned = max(min(cos_theta_unsigned, 1), -1);
angle_unsigned = acosd(cos_theta_unsigned);

fprintf('\n====================================================\n');
fprintf('%s %s pre-orthogonalization choice-axis comparison:\n', areaName, durName);
fprintf('Number of shared neurons used             = %d\n', length(u));
fprintf('Vestibular vs Visual raw angle (0-180)    = %.4f deg\n', angle_raw);
fprintf('Vestibular vs Visual unsigned angle (0-90)= %.4f deg\n', angle_unsigned);
fprintf('====================================================\n');

%% ========================================================================
% save output
%% ========================================================================
ChoiceAxisComparison = struct();
ChoiceAxisComparison.area = areaName;
ChoiceAxisComparison.dur  = durName;
ChoiceAxisComparison.raw_angle_deg      = angle_raw;
ChoiceAxisComparison.unsigned_angle_deg = angle_unsigned;
ChoiceAxisComparison.n_shared_neurons   = length(u);
ChoiceAxisComparison.choice_vest_raw    = choice_vest;
ChoiceAxisComparison.choice_vis_raw     = choice_vis;
ChoiceAxisComparison.valid_id           = valid_id;
ChoiceAxisComparison.AllBetasRaw        = AllBetasRaw;

OutputPath = 'Z:\Users\TDR\results\tTDR_choiceAngle_Vis_vs_Vest\';
if ~exist(OutputPath, 'dir')
    mkdir(OutputPath);
end

SaveFileName = fullfile(OutputPath, ...
    sprintf('%s_%s_ChoiceAxis_Vest_vs_Vis.mat', areaName, durName));

save(SaveFileName, 'ChoiceAxisComparison');
fprintf('Saved results to:\n%s\n', SaveFileName);

end