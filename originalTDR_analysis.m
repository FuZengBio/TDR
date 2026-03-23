clear;clc

FilePath='Z:\Users\TDR\Data\';
OutputPath='Z:\Users\TDR\results\originalTDR\';

structures = [1 2 3]; %1=VIP; 2=MSTd; 3=PIVC;
structure_names = {'VIP', 'MSTd', 'PIVC'};

durations = [1 2];    %1=1s; 2=2s
duration_names = {'1s','2s'};

Kernel_method=1;

%%% change order
reverse = 0; %0=heading->choice; 1=choice->heading

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

        stim_type=[1,2];

        for ss=1:length(stim_type)

            clear Test_REG_r Test_PCA_Data REG_bT projection Summary
            clear STPP SV KernelData

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Part 1: recompute sdf from raster
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for nCell=1:size(REG_raster_per_trial,1)
                clear temp_spike
                temp_spike=REG_raster_per_trial{nCell,ss};
                clear sdf

                for i=1:size(temp_spike,2)
                    clear spike
                    spike=double(squeeze(temp_spike(:,i,:)));
                    [ntrials,nsamples] = size(spike); %#ok<ASGLU>

                    for ii = 1:ntrials
                        convspike = conv(spike(ii,:), window, 'same');
                        convspike = convspike(crop_idx);
                        tempdata = downsample(convspike,10);
                        sdf((i-1)*ntrials+ii,:) = tempdata;
                    end
                end
                Test_REG_r{nCell,1}=sdf;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Part 2: setup regression parameters
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            span=10;
            REG_t=square_sum(t1,span)/span;

            REG_params_stim={'choice','FUNangle','ONES'};
            REG_choice_dim=find(strcmp(REG_params_stim,'choice'));
            REG_head_dim=find(strcmp(REG_params_stim,'FUNangle'));

            if reverse==0
                id=[REG_head_dim REG_choice_dim];   % heading -> choice
            elseif reverse==1
                id=[REG_choice_dim REG_head_dim];   % choice -> heading
            else
                error('reverse must be 0 or 1 for original TDR');
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Part 3: regression at each time point
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nCells = size(Test_REG_r,1);
            nTime  = size(Test_REG_r{1},2);
            nParam = length(REG_params_stim);

            REG_bT = nan(nCells,nTime,nParam);

            for Ni=1:nCells

                heading_trials = REG_FUNangle{Ni,ss};
                choice_trials  = REG_choice{Ni,ss};

                % PIVC: recode choice to -1 / +1 as in original code
                if strcmp(structure_names{s}, 'PIVC')
                    choice_trials = 2*(choice_trials - 1.5);
                end

                if reverse==0
                    X = [heading_trials, choice_trials, ones(length(heading_trials),1)];
                else
                    X = [choice_trials, heading_trials, ones(length(heading_trials),1)];
                end

                for tt = 1:nTime
                    y = Test_REG_r{Ni}(:,tt);
                    beta = X \ y;

                    REG_bT(Ni,tt,id(1)) = beta(1);
                    REG_bT(Ni,tt,id(2)) = beta(2);
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Part 4: build Test_PCA_Data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for nCell=1:size(REG_raster_per_trial,1)
                clear temp_spike
                temp_spike=REG_raster_per_trial{nCell,ss};
                clear sdf

                for i=1:size(temp_spike,1)
                    clear spike
                    spike=double(squeeze(temp_spike(i,:,:)));
                    [ntrials,nsamples] = size(spike); %#ok<ASGLU>
                    clear temp_mean

                    for ii = 1:ntrials
                        convspike = conv(spike(ii,:), window, 'same');
                        convspike = convspike(crop_idx);
                        tempdata = downsample(convspike,10);
                        temp_mean(ii,:) = tempdata;
                    end

                    Test_PCA_Data(1,i).A(:,nCell) = mean(temp_mean)';
                    Test_PCA_Data(1,i).times = REG_t;
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Part 5: run TDR
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear DataRun
            DataRun=Test_PCA_Data;
            clear PCA_params

            beta_first  = squeeze(REG_bT(:,:,id(1)));
            beta_second = squeeze(REG_bT(:,:,id(2)));

            PCA_params.B = {beta_first, beta_second};
            PCA_params.normalize=false;
            PCA_params.numPCs = 12;
            times=[];

            fprintf('\nRunning %s %s, ss=%d\n', structure_names{s}, duration_names{d}, ss);
            fprintf('nCells in DataRun = %d\n', size(Test_PCA_Data(1).A,2));
            fprintf('size(beta_first)  = [%d %d]\n', size(beta_first,1), size(beta_first,2));
            fprintf('size(beta_second) = [%d %d]\n', size(beta_second,1), size(beta_second,2));

            [projection,Summary]=TDR_time(DataRun, times, PCA_params);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Part 6: plot settings
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Flip_plot=0;
            plot_times = [];
            f(1)=figure; set(gcf,'Units','normalized','color','w','position',[0 0.05 1 0.85]); %#ok<SAGROW>
            plotParams.reusePlot=1;
            plotParams.substPCs = 2;
            plotParams.planes2plot=1;

            REG_params_use=REG_params_stim(id);
            plotParams.RegNames=REG_params_use;
            plotParams.plotPlanEllipse=false;
            plotParams.planMarkerSize=3;
            plotParams.lineWidth=2;
            plotParams.times=plot_times;
            plotParams.params2plot =[1 2];

            graded=1;
            cuecolor=1;

            if ss==1
                if strcmp(duration_names{d}, '2s')
                    condsRun=1:9;
                elseif strcmp(structure_names{s}, 'PIVC')
                    condsRun=1:10;
                else
                    condsRun=1:7;
                end
                stimtype='Vest';
            elseif ss==2
                if strcmp(duration_names{d}, '2s')
                    condsRun=13:21;
                elseif strcmp(structure_names{s}, 'PIVC')
                    condsRun=13:22;
                else
                    condsRun=13:19;
                end
                stimtype='Vis';
            end

            plotParams.colors=plot_colors(condsRun,graded,cuecolor);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % plot the variance / legends
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fontsize=12;
            clear legR legT

            if strcmp(duration_names{d}, '2s')
                tempAngle=[-9,-3.6,-1.44,-0.58,0,0.58,1.44,3.6,9];
            elseif strcmp(structure_names{s}, 'PIVC')
                tempAngle=[-16, -8, -4, -2, -1, 1, 2, 4, 8, 16];
            else
                tempAngle=[-12, -6, -1.5, 0, 1.5, 6, 12];
            end

            for i=1:length(condsRun)
                legR{i}=sprintf('%g%s',round(10*(tempAngle(i)))/10,'{\circ}');
            end
            for i=1:length(condsRun)
                legT{i}=sprintf('%g%s',round(10*(tempAngle(i)))/10,'{\circ}');
            end

            if length(REG_params_use)>3
                b_comb=[1 2; 2 3; 2 4; 3 4];
            else
                b_comb = nchoosek(1:length(REG_params_use),2);
            end

            cols=max(length(b_comb),length(REG_params_use) +1);
            rows=2;
            enl=0;

            for comb=1:size(b_comb,1)
                b(1)=b_comb(comb,1); b(2)=b_comb(comb,2);
                if Flip_plot
                    plotParams.params2plot=[b(2) b(1)];
                else
                    plotParams.params2plot=[b(1) b(2)];
                end

                subplot(rows,cols,comb); hold on;
                POS=get(gca,'position');
                set(gca,'position',[POS(1)-enl/2-0.05 POS(2)-enl/2 POS(3)+enl POS(4)+enl])
                plotParams.colors=plot_colors(condsRun,graded,cuecolor);
                set(gca,'visible','off','fontsize',fontsize)

                [colorStruct, haxP, vaxP, phR] = phaseSpace(projection, Summary, plotParams); %#ok<ASGLU>
                if comb==size(b_comb,1)
                    l=legend(phR(end:-1:1),legR{end:-1:1},'location','east');
                    POS=get(l,'position');
                    set(l,'position',[POS(1)+0.05 POS(2) POS(3) POS(4)])
                end
                set(findobj('type','line'),'clipping','off')
                axis tight
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % plot Heading / Choice vs. time
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            REG_params_names=nameREG(REG_params_use,1);
            for par=1:length(REG_params_use)
                subplot(rows,cols,cols + par); hold on;
                POS=get(gca,'position');
                set(gca,'position',[POS(1)-0.05 POS(2) POS(3) POS(4)])
                set(gca,'fontsize',fontsize)

                for c=1:length(condsRun)
                    if strcmp(duration_names{d}, '2s')
                        STPP{c,par}=smooth(projection(c).TargetedPCAproj(:,par),20);
                    elseif strcmp(structure_names{s}, 'PIVC')
                        STPP{c,par}=smooth(projection(c).TargetedPCAproj(:,par),10);
                    else
                        STPP{c,par}=smooth(projection(c).TargetedPCAproj(:,par),5);
                    end
                    plot(projection(c).times,STPP{c,par},'-','color',plotParams.colors{c},'linewidth',plotParams.lineWidth);
                    m(c)=max(abs(projection(c).TargetedPCAproj(:,par))); %#ok<SAGROW>
                end
                ylabel(sprintf('%s [A.U.]',REG_params_names{par}),'fontsize',fontsize)
                xlabel('time [ms]')

                if strcmp(duration_names{d}, '2s')
                    lim=[250 200 125];
                elseif strcmp(structure_names{s}, 'PIVC')
                    lim=[250 200 125];
                else
                    lim=[150 100 125];
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % plot variance vs. time
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot(rows,cols,cols + length(REG_params_use) +1); hold on;
            POS=get(gca,'position');
            set(gca,'position',[POS(1)-enl/2-0.05 POS(2)-enl/2 POS(3)+enl POS(4)+enl])
            set(gca,'LineStyleOrder',{'-','--',':','-.'},'ColorOrder',[0 0 0],'color','none','fontsize',fontsize); hold on;

            clear I SV
            for I=1:length(REG_params_use)
                if strcmp(duration_names{d}, '2s')
                    SV(:,I) = smooth(Summary.varCaptEachTargetedPCT(:,I),5);
                elseif strcmp(structure_names{s}, 'PIVC')
                    SV(:,I) = smooth(Summary.varCaptEachTargetedPCT(:,I),5);
                else
                    SV(:,I) = smooth(Summary.varCaptEachTargetedPCT(:,I),5);
                end
            end

            plot(REG_t,SV,'linewidth',plotParams.lineWidth)

            ylabel('Variance Explained')
            xlabel('time [ms]')
            l=legend(REG_params_names);
            POS=get(l,'position');
            set(l,'position',[POS(1)+0.05 POS(2) POS(3) POS(4)])
            yl=ylim;
            for i=1:size(b_comb,1)
                text(sti_dur,(size(b_comb,1)-i+1)*(yl(2)-yl(1))/(2*size(b_comb,1)), ...
                    sprintf('{\\angle}(%s,%s)=%u{\\circ}', ...
                    REG_params_names{b_comb(i,1)}, ...
                    REG_params_names{b_comb(i,2)}, ...
                    round(Summary.Btheta(b_comb(i,1),b_comb(i,2)))))
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % output
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            clear KernelData
            KernelData.projection.TargetedPCAproj = STPP;
            KernelData.VarianceTime = SV;
            KernelData.projectionNegative = projection(1).TargetedPCAproj;
            KernelData.projectionPositive = projection(length(condsRun)).TargetedPCAproj;
            KernelData.order=REG_params_use;
            KernelData.UniqueVariance = Summary.varCaptEachTargetedPC;
            KernelData.VarianceTimecourse = Summary.varCaptEachTargetedPCT;

            if ~exist(OutputPath,'dir')
                mkdir(OutputPath);
            end

            SaveFileName = fullfile(OutputPath, ...
                [stimtype '_' structure_names{s} '_' duration_names{d} '_' cell2mat(REG_params_use) '.mat']);
            save(SaveFileName, 'KernelData');

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % save figure
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            set(gcf,'paperpositionmode','auto');
            filename = fullfile(OutputPath, ...
                [stimtype '_' structure_names{s} '_' duration_names{d} '_Kernel_' num2str(Kernel_method) '_' strjoin(REG_params_use,'_') '.png']);
            saveas(gcf,filename,'png');
            close;
        end
    end
end