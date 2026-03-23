function HC_choiceAxisAngle_Vis_vs_Vest_batch()
% Compute visual vs vestibular CHOICE-axis angle
% under original HC TDR (pre-orthogonalization choice beta)
%
% Output:
%   one .mat file per dataset
%   one figure per dataset

clear; clc;

%% ========================================================================
% Unified paths
% ========================================================================
FilePath   = 'Z:\Users\TDR\Data\';
ResultRoot = 'Z:\Users\TDR\results\originalTDR_timecourse_decoding\';

if ~exist(ResultRoot, 'dir')
    mkdir(ResultRoot);
end

%% ========================================================================
% Dataset list
% ========================================================================
datasets = { ...
    struct('area','VIP',  'dur','1s', 'mat','VIP_Data_1s_REG.mat'), ...
    struct('area','VIP',  'dur','2s', 'mat','VIP_Data_2s_REG.mat'), ...
    struct('area','MSTd', 'dur','1s', 'mat','MSTd_Data_1s_REG.mat'), ...
    struct('area','MSTd', 'dur','2s', 'mat','MSTd_Data_2s_REG.mat'), ...
    struct('area','PIVC', 'dur','1s', 'mat','PIVC_Data_1s_REG.mat') ...
    };

%% ========================================================================
% Shared parameters
% ========================================================================
sigma = 30;
half_guas = sigma * 3;
dist = (-(half_guas):(half_guas))';
window = 1/(sigma*sqrt(2*pi)) * exp(-((dist.^2)/(2*sigma^2))) * 1000;

stim_type = [1, 2];   % 1 = vestibular, 2 = visual
span = 10;

REG_params_stim = {'choice','FUNangle','ONES'};
REG_choice_dim = find(strcmp(REG_params_stim,'choice'));
REG_head_dim   = find(strcmp(REG_params_stim,'FUNangle')); %#ok<NASGU>

reverse = 0;   % keep same logic as your example

%% ========================================================================
% Run all datasets
% ========================================================================
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

    load(matfile);

    %% --------------------------------------------------------------------
    % Dataset-specific timing
    % ---------------------------------------------------------------------
    switch ds.dur
        case '2s'
            t1 = 101:2000;
            crop_idx = 101:2000;

        case '1s'
            t1 = 101:1000;

            if strcmp(ds.area,'VIP') || strcmp(ds.area,'MSTd')
                crop_idx = 1186:2085;
            elseif strcmp(ds.area,'PIVC')
                crop_idx = 101:1000;
            else
                error('Unknown area for 1s dataset: %s', ds.area);
            end

        otherwise
            error('Unknown duration: %s', ds.dur);
    end

    REG_t = square_sum(t1, span) / span;

    %% --------------------------------------------------------------------
    % Store choice-axis matrices
    % Each entry: neurons x time
    % ---------------------------------------------------------------------
    ChoiceAxisHC = cell(1,2);   % {1}=Vest, {2}=Vis

    %% --------------------------------------------------------------------
    % Main loop over modality
    % ---------------------------------------------------------------------
    for ss = 1:length(stim_type)

        fprintf('\n---------------------------------------------\n');
        if ss == 1
            fprintf('Processing vestibular condition...\n');
        else
            fprintf('Processing visual condition...\n');
        end

        %% ----------------------------------------------------------------
        % Part 1: Process spike data into SDF
        % -----------------------------------------------------------------
        clear Test_REG_r
        Test_REG_r = cell(size(REG_raster_per_trial,1),1);

        for nCell = 1:size(REG_raster_per_trial,1)

            temp_spike = REG_raster_per_trial{nCell, ss};

            if isempty(temp_spike)
                Test_REG_r{nCell,1} = [];
                continue
            end

            clear sdf
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

        %% ----------------------------------------------------------------
        % Part 2: Per-timepoint regression
        % Only keep choice beta across neurons over time
        % -----------------------------------------------------------------
        Ncells = size(Test_REG_r,1);

        % determine number of time points
        Ntimepoints = [];
        for nCell = 1:Ncells
            if ~isempty(Test_REG_r{nCell})
                Ntimepoints = size(Test_REG_r{nCell},2);
                break
            end
        end
        if isempty(Ntimepoints)
            error('Could not determine number of time points from Test_REG_r.');
        end

        Nparams = length(REG_params_stim);
        REG_bT = nan(Ncells, Ntimepoints, Nparams);

        for Ni = 1:Ncells

            if isempty(Test_REG_r{Ni})
                continue
            end

            % task variables for this neuron and modality
            heading_trials = REG_FUNangle{Ni, ss};
            choice_trials  = REG_choice{Ni, ss};

            if isempty(heading_trials) || isempty(choice_trials)
                continue
            end

            heading_trials = heading_trials(:);
            choice_trials  = choice_trials(:);

            if reverse == 0
                % columns: heading, choice, constant
                X = [heading_trials, choice_trials, ones(length(heading_trials),1)];
                col_head = 1; %#ok<NASGU>
                col_choice = 2;
            else
                % columns: choice, heading, constant
                X = [choice_trials, heading_trials, ones(length(choice_trials),1)];
                col_choice = 1;
                col_head = 2; %#ok<NASGU>
            end

            % SDF matrix: trials x time
            r_all = Test_REG_r{Ni};

            % make sure trial counts match
            nTrials_reg = size(X,1);
            nTrials_sdf = size(r_all,1);
            nUse = min(nTrials_reg, nTrials_sdf);

            if nUse < 2
                continue
            end

            X_use = X(1:nUse,:);
            r_use = r_all(1:nUse,:);

            for t = 1:Ntimepoints
                y = r_use(:,t);

                if all(isnan(y)) || any(isnan(X_use(:)))
                    continue
                end

                beta_t = X_use \ y;
                REG_bT(Ni, t, REG_choice_dim) = beta_t(col_choice);
            end
        end

        % save ONLY the pre-orthogonalization choice-axis matrix
        % size: neurons x time
        ChoiceAxisHC{ss} = squeeze(REG_bT(:,:,REG_choice_dim));

        fprintf('Done. Choice axis matrix size = [%d neurons x %d time points]\n', ...
            size(ChoiceAxisHC{ss},1), size(ChoiceAxisHC{ss},2));
    end

    %% --------------------------------------------------------------------
    % Part 3: Compare vestibular vs visual choice axes over time
    % ---------------------------------------------------------------------
    choice_vest_mat = ChoiceAxisHC{1};   % neurons x time
    choice_vis_mat  = ChoiceAxisHC{2};   % neurons x time

    if isempty(choice_vest_mat) || isempty(choice_vis_mat)
        warning('ChoiceAxisHC is empty for one or both modalities. Skip %s %s.', ...
            ds.area, ds.dur);
        continue;
    end

    if size(choice_vest_mat,2) ~= size(choice_vis_mat,2)
        warning('Vestibular and visual choice matrices have different time lengths. Skip %s %s.', ...
            ds.area, ds.dur);
        continue;
    end

    nTime = size(choice_vest_mat,2);

    angle_raw_t = nan(1, nTime);         % 0~180 deg
    angle_unsigned_t = nan(1, nTime);    % 0~90 deg
    n_shared_t = nan(1, nTime);

    for t = 1:nTime
        u = choice_vest_mat(:,t);
        v = choice_vis_mat(:,t);

        valid_id = ~isnan(u) & ~isnan(v);
        u = u(valid_id);
        v = v(valid_id);

        n_shared_t(t) = length(u);

        if isempty(u) || norm(u)==0 || norm(v)==0
            continue
        end

        u = u / norm(u);
        v = v / norm(v);

        % raw angle: 0~180 deg
        cos_theta_raw = dot(u,v);
        cos_theta_raw = max(min(cos_theta_raw,1),-1);
        angle_raw_t(t) = acosd(cos_theta_raw);

        % unsigned angle: 0~90 deg
        cos_theta_unsigned = abs(dot(u,v));
        cos_theta_unsigned = max(min(cos_theta_unsigned,1),-1);
        angle_unsigned_t(t) = acosd(cos_theta_unsigned);
    end

    %% --------------------------------------------------------------------
    % Part 4: Summary output
    % ---------------------------------------------------------------------
    fprintf('\n=============================================\n');
    fprintf('%s %s: Visual vs Vestibular CHOICE-axis comparison\n', ds.area, ds.dur);
    fprintf('Number of time points = %d\n', nTime);
    fprintf('Mean unsigned angle   = %.4f deg\n', mean(angle_unsigned_t,'omitnan'));
    fprintf('Min unsigned angle    = %.4f deg\n', min(angle_unsigned_t,[],'omitnan'));
    fprintf('Max unsigned angle    = %.4f deg\n', max(angle_unsigned_t,[],'omitnan'));
    fprintf('=============================================\n');

    %% --------------------------------------------------------------------
    % Part 5: Save output
    % ---------------------------------------------------------------------
    HCChoiceAngle = struct();

    HCChoiceAngle.area = ds.area;
    HCChoiceAngle.dur  = ds.dur;

    HCChoiceAngle.choice_vest_mat = choice_vest_mat;   % neurons x time
    HCChoiceAngle.choice_vis_mat  = choice_vis_mat;    % neurons x time

    HCChoiceAngle.angle_raw_t = angle_raw_t;           % 1 x time
    HCChoiceAngle.angle_unsigned_t = angle_unsigned_t; % 1 x time
    HCChoiceAngle.n_shared_t = n_shared_t;             % 1 x time
    HCChoiceAngle.times = REG_t;                       % 1 x time

    HCChoiceAngle.mean_unsigned_angle_deg = mean(angle_unsigned_t,'omitnan');
    HCChoiceAngle.min_unsigned_angle_deg  = min(angle_unsigned_t,[],'omitnan');
    HCChoiceAngle.max_unsigned_angle_deg  = max(angle_unsigned_t,[],'omitnan');

    OutputPath = fullfile(ResultRoot, [ds.area '_' ds.dur], 'HC_choice_angle');
    if ~exist(OutputPath, 'dir')
        mkdir(OutputPath);
    end

    SaveFileName = fullfile(OutputPath, ...
        sprintf('%s_%s_HC_choiceAxisAngle_Vis_vs_Vest.mat', ds.area, ds.dur));

    save(SaveFileName, 'HCChoiceAngle');
    fprintf('Saved results to:\n%s\n', SaveFileName);

    %% --------------------------------------------------------------------
    % Part 6: Optional plot
    % ---------------------------------------------------------------------
    f = figure('color','w','Name',sprintf('%s_%s_HC_choice_angle', ds.area, ds.dur));
    plot(REG_t, angle_unsigned_t, 'k-', 'LineWidth', 2); hold on
    xlabel('time [ms]');
    ylabel('Visual vs Vestibular choice angle [deg]');
    title(sprintf('%s %s: HC choice-axis similarity (unsigned angle)', ds.area, ds.dur));
    ylim([0 90]);
    box off

    saveas(f, fullfile(OutputPath, ...
        sprintf('%s_%s_HC_choiceAxisAngle_Vis_vs_Vest.png', ds.area, ds.dur)));
    close(f);

end

fprintf('\nAll datasets finished.\n');

end