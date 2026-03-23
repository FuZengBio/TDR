clear; clc;
FilePath = 'Z:\Users\TDR\Data\';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% datasets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
structure_names = {'VIP','MSTd','PIVC'};
duration_names  = {'1s','2s'};

structures = 1:length(structure_names);
durations  = 1:length(duration_names);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fixed parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
span = 10;
reverse = 0;   % 0 = [cue, heading], 1 = [heading, cue]
Kernel_method = 1;

for s = 1:length(structures)
    for d = 1:length(durations)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % skip non-existing dataset: PIVC only has 1s
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(structure_names{s}, 'PIVC') && strcmp(duration_names{d}, '2s')
            continue
        end

        clear REG_raster_per_trial REG_Stim REG_choice REG_FUNangle
        clear Test_REG_r Test_PCA_Data REG_bT projection Summary
        clear STPP SV KernelData

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % load dataset
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(structure_names{s}, 'VIP') && strcmp(duration_names{d}, '1s')
            matfile = fullfile(FilePath, 'VIP_Data_1s_REG.mat');
        elseif strcmp(structure_names{s}, 'VIP') && strcmp(duration_names{d}, '2s')
            matfile = fullfile(FilePath, 'VIP_Data_2s_REG.mat');
        elseif strcmp(structure_names{s}, 'MSTd') && strcmp(duration_names{d}, '1s')
            matfile = fullfile(FilePath, 'MSTd_Data_1s_REG.mat');
        elseif strcmp(structure_names{s}, 'MSTd') && strcmp(duration_names{d}, '2s')
            matfile = fullfile(FilePath, 'MSTd_Data_2s_REG.mat');
        elseif strcmp(structure_names{s}, 'PIVC') && strcmp(duration_names{d}, '1s')
            matfile = fullfile(FilePath, 'PIVC_Data_1s_REG.mat');
        else
            error('Unknown dataset');
        end

        load(matfile);

        fprintf('\n====================================================\n');
        fprintf('Running %s %s\n', structure_names{s}, duration_names{d});
        fprintf('Loading file: %s\n', matfile);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gaussian kernel
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sigma=30;
        half_guas=sigma*3;
        dist=(-(half_guas):(half_guas))';
        window=1/(sigma*sqrt(2*pi))*exp(-((dist.^2)/(2*sigma^2)))*1000;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % dataset-specific settings
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(duration_names{d}, '2s')
            t1 = 101:2000;
            sti_dur = 2000;
            crop_idx = 101:2000;
        else
            t1 = 101:1000;
            sti_dur = 1000;

            if strcmp(structure_names{s}, 'VIP') || strcmp(structure_names{s}, 'MSTd')
                crop_idx = 1186:2085;
            else
                crop_idx = 101:1000;
            end
        end

        REG_t = square_sum(t1, span) / span;
        Ncells = size(REG_raster_per_trial, 1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % heading list and smoothing parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(duration_names{d}, '2s')
            tempAngle   = [-9 -3.6 -1.44 -0.58 0 0.58 1.44 3.6 9];
            smooth_proj = 20;
            smooth_var  = 5;
        else
            if strcmp(structure_names{s}, 'PIVC')
                tempAngle   = [-16 -8 -4 -2 -1 1 2 4 8 16];
                smooth_proj = 10;
                smooth_var  = 10;
            else
                tempAngle   = [-12 -6 -1.5 0 1.5 6 12];
                smooth_proj = 5;
                smooth_var  = 5;
            end
        end

        nAngles = length(tempAngle);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % regressors: cue + heading
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        REG_params_stim = {'Stim','FUNangle'};
        REG_cue_dim  = find(strcmp(REG_params_stim, 'Stim'));
        REG_head_dim = find(strcmp(REG_params_stim, 'FUNangle'));

        if reverse == 0
            id = [REG_cue_dim REG_head_dim];
        else
            id = [REG_head_dim REG_cue_dim];
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Part 1. pooled SDF for regression
        % pool vestibular and visual trials together
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Test_REG_r = cell(Ncells,1);

        for nCell = 1:Ncells
            sdf_all = [];

            for ss = 1:2   % 1 = vestibular, 2 = visual
                temp_spike = REG_raster_per_trial{nCell, ss};

                if isempty(temp_spike)
                    continue
                end

                for i = 1:size(temp_spike, 2)
                    spike = double(squeeze(temp_spike(:, i, :)));   % nAngles x T
                    [ntrials, ~] = size(spike);

                    for ii = 1:ntrials
                        convspike = conv(spike(ii,:), window, 'same');
                        convspike = convspike(crop_idx);
                        tempdata  = downsample(convspike, 10);
                        sdf_all   = [sdf_all; tempdata];
                    end
                end
            end

            Test_REG_r{nCell} = sdf_all;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Part 2. time-resolved regression
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(Test_REG_r) && ~isempty(Test_REG_r{1})
            Ntimepoints = size(Test_REG_r{1}, 2);
        else
            Ntimepoints = numel(REG_t);
        end

        Nparams = length(REG_params_stim);
        REG_bT  = nan(Ncells, Ntimepoints, Nparams);

        for Ni = 1:Ncells

            if isempty(REG_Stim{Ni,1}) || isempty(REG_Stim{Ni,2}) || ...
               isempty(REG_FUNangle{Ni,1}) || isempty(REG_FUNangle{Ni,2})
                continue
            end

            cue_raw        = [REG_Stim{Ni,1}(:); REG_Stim{Ni,2}(:)];
            cue_trials     = (cue_raw == 1) * 2 - 1;   % vest=+1, vis=-1
            heading_trials = [REG_FUNangle{Ni,1}(:); REG_FUNangle{Ni,2}(:)];

            if reverse == 0
                X = [cue_trials, heading_trials, ones(length(heading_trials),1)];
            else
                X = [heading_trials, cue_trials, ones(length(heading_trials),1)];
            end

            r_all = Test_REG_r{Ni};

            if size(r_all,1) ~= numel(heading_trials)
                warning('Cell %d: trial mismatch between sdf (%d) and regressors (%d). Skipped.', ...
                    Ni, size(r_all,1), numel(heading_trials));
                continue
            end

            for t = 1:Ntimepoints
                y = r_all(:,t);
                beta_t = X \ y;

                REG_bT(Ni,t,id(1)) = beta_t(1);
                REG_bT(Ni,t,id(2)) = beta_t(2);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Part 3. build Test_PCA_Data
        % conditions:
        %   1:nAngles = vestibular
        %   nAngles+1:2*nAngles = visual
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Test_PCA_Data = struct([]);

        for nCell = 1:Ncells

            temp_spike_ves = REG_raster_per_trial{nCell,1};
            temp_spike_vis = REG_raster_per_trial{nCell,2};

            if isempty(temp_spike_ves) || isempty(temp_spike_vis)
                continue
            end

            if size(temp_spike_ves,1) ~= nAngles || size(temp_spike_vis,1) ~= nAngles
                warning('Cell %d: angle count mismatch, skipped.', nCell);
                continue
            end

            for ang = 1:nAngles

                % ---------------- vestibular ----------------
                spike_ves = double(squeeze(temp_spike_ves(ang,:,:)));   % nRep x T
                [ntrials_ves, ~] = size(spike_ves);
                temp_mean_ves = zeros(ntrials_ves, numel(REG_t));

                for ii = 1:ntrials_ves
                    convspike = conv(spike_ves(ii,:), window, 'same');
                    convspike = convspike(crop_idx);
                    tempdata  = downsample(convspike, 10);
                    temp_mean_ves(ii,:) = tempdata;
                end

                Test_PCA_Data(1, ang).A(:,nCell) = mean(temp_mean_ves, 1, 'omitnan')';
                Test_PCA_Data(1, ang).times      = REG_t;

                % ---------------- visual ----------------
                spike_vis = double(squeeze(temp_spike_vis(ang,:,:)));   % nRep x T
                [ntrials_vis, ~] = size(spike_vis);
                temp_mean_vis = zeros(ntrials_vis, numel(REG_t));

                for ii = 1:ntrials_vis
                    convspike = conv(spike_vis(ii,:), window, 'same');
                    convspike = convspike(crop_idx);
                    tempdata  = downsample(convspike, 10);
                    temp_mean_vis(ii,:) = tempdata;
                end

                idx_vis = nAngles + ang;
                Test_PCA_Data(1, idx_vis).A(:,nCell) = mean(temp_mean_vis, 1, 'omitnan')';
                Test_PCA_Data(1, idx_vis).times      = REG_t;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Part 4. run TDR
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        beta_first  = squeeze(REG_bT(:,:,id(1)));
        beta_second = squeeze(REG_bT(:,:,id(2)));

        PCA_params.B         = {beta_first, beta_second};
        PCA_params.normalize = false;
        PCA_params.numPCs    = 12;

        DataRun = Test_PCA_Data;
        times   = [];

        fprintf('nCells in DataRun = %d\n', size(Test_PCA_Data(1).A,2));
        fprintf('size(beta_first)  = [%d %d]\n', size(beta_first,1), size(beta_first,2));
        fprintf('size(beta_second) = [%d %d]\n', size(beta_second,1), size(beta_second,2));

        [projection, Summary] = TDR_time(DataRun, times, PCA_params);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Part 5. plot
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Flip_plot = 0;
        plot_times = [];

        f = figure;
        set(gcf, 'Units','normalized', 'color','w', 'position',[0 0.05 1 0.85]);

        plotParams.reusePlot       = 1;
        plotParams.substPCs        = 2;
        plotParams.planes2plot     = 1;
        plotParams.plotPlanEllipse = false;
        plotParams.planMarkerSize  = 3;
        plotParams.lineWidth       = 2;
        plotParams.times           = plot_times;

        REG_params_use         = REG_params_stim(id);
        plotParams.RegNames    = REG_params_use;
        plotParams.params2plot = [2 1];

        graded   = 1;
        cuecolor = 1;
        condsRun = 1:(2*nAngles);
        stimtype = 'All';
        plotParams.colors = plot_colors(condsRun, graded, cuecolor);

        fontsize = 12;
        clear legR legT
        tempAngle2 = [tempAngle tempAngle];

        for i = 1:length(condsRun)
            legR{i} = sprintf('%g%s', round(10*(tempAngle2(i)))/10, '{\circ}');
            legT{i} = sprintf('%g%s', round(10*(tempAngle2(i)))/10, '{\circ}');
        end

        if length(REG_params_use) > 3
            b_comb = [1 2; 2 3; 2 4; 3 4];
        else
            b_comb = nchoosek(1:length(REG_params_use), 2);
        end

        cols = max(length(b_comb), length(REG_params_use) + 1);
        rows = 2;
        enl  = 0;

        % ---------------- phase-space ----------------
        for comb = 1:size(b_comb,1)
            b = b_comb(comb,:);

            if Flip_plot
                plotParams.params2plot = [b(2) b(1)];
            else
                plotParams.params2plot = [b(1) b(2)];
            end

            subplot(rows, cols, comb); hold on;
            POS = get(gca, 'position');
            set(gca, 'position', [POS(1)-enl/2-0.05 POS(2)-enl/2 POS(3)+enl POS(4)+enl]);
            set(gca, 'visible','off', 'fontsize',fontsize);

            [~,~,~,phR] = phaseSpace_ZF(projection, Summary, plotParams);

            if comb == size(b_comb,1)
                l = legend(phR(end:-1:1), legR{end:-1:1}, 'location','east');
                POS = get(l, 'position');
                set(l, 'position', [POS(1)+0.05 POS(2) POS(3) POS(4)]);
            end
            set(findobj('type','line'),'clipping','off');
            axis tight
        end

        % ---------------- projection time course ----------------
        REG_params_names = nameREG(REG_params_use,1);
        clear STPP

        for par = 1:length(REG_params_use)
            subplot(rows, cols, cols + par); hold on;
            POS = get(gca, 'position');
            set(gca, 'position', [POS(1)-0.05 POS(2) POS(3) POS(4)]);
            set(gca, 'fontsize', fontsize);

            for c = 1:length(condsRun)
                STPP{c,par} = smooth(projection(c).TargetedPCAproj(:,par), smooth_proj);
                plot(projection(c).times, STPP{c,par}, '-', ...
                    'color', plotParams.colors{c}, ...
                    'linewidth', plotParams.lineWidth);
            end

            ylabel(sprintf('%s [A.U.]', REG_params_names{par}), 'fontsize', fontsize);
            xlabel('time [ms]');
        end

        % ---------------- variance explained ----------------
        subplot(rows, cols, cols + length(REG_params_use) + 1); hold on;
        POS = get(gca, 'position');
        set(gca, 'position', [POS(1)-enl/2-0.05 POS(2)-enl/2 POS(3)+enl POS(4)+enl]);
        set(gca, 'LineStyleOrder', {'-','--',':','-.'}, ...
            'ColorOrder', [0 0 0], 'color', 'none', 'fontsize', fontsize);
        hold on;

        clear SV
        for I = 1:length(REG_params_use)
            SV(:,I) = smooth(Summary.varCaptEachTargetedPCT(:,I), smooth_var);
        end
        plot(REG_t, SV, 'linewidth', plotParams.lineWidth);

        ylabel('Variance Explained');
        xlabel('time [ms]');
        l = legend(REG_params_names);
        POS = get(l, 'position');
        set(l, 'position', [POS(1)+0.05 POS(2) POS(3) POS(4)]);

        yl = ylim;
        for i = 1:size(b_comb,1)
            text(sti_dur, ...
                (size(b_comb,1)-i+1)*(yl(2)-yl(1))/(2*size(b_comb,1)), ...
                sprintf('{\\angle}(%s,%s)=%u{\\circ}', ...
                REG_params_names{b_comb(i,1)}, ...
                REG_params_names{b_comb(i,2)}, ...
                round(Summary.Btheta(b_comb(i,1), b_comb(i,2)))));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Part 6. save
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear KernelData
        KernelData.projection.TargetedPCAproj = STPP;
        KernelData.VarianceTime               = SV;
        KernelData.projectionNegative         = projection(1).TargetedPCAproj;
        KernelData.projectionPositive         = projection(length(condsRun)).TargetedPCAproj;
        KernelData.order                      = REG_params_use;
        KernelData.UniqueVariance             = Summary.varCaptEachTargetedPC;
        KernelData.VarianceTimecourse         = Summary.varCaptEachTargetedPCT;
        KernelData.area                       = structure_names{s};
        KernelData.duration                   = duration_names{d};
        KernelData.angles                     = tempAngle;
        KernelData.REG_t                      = REG_t;

        SavePath = 'Z:\Users\TDR\results\originalTDR_CueHeading\';
        if ~exist(SavePath, 'dir')
            mkdir(SavePath);
        end

        SaveFileName = fullfile(SavePath, ...
            sprintf('%s_%s_originalCueHeadingTDR.mat', structure_names{s}, duration_names{d}));
        save(SaveFileName, 'KernelData');

        set(gcf, 'paperpositionmode', 'auto');
        FigName = fullfile(SavePath, ...
            sprintf('%s_%s_originalCueHeadingTDR_Kernel_%d.png', ...
            structure_names{s}, duration_names{d}, Kernel_method));
        saveas(gcf, FigName, 'png');

        close(gcf);

        fprintf('Saved to:\n%s\n%s\n', SaveFileName, FigName);
    end
end

fprintf('\nAll datasets finished.\n');