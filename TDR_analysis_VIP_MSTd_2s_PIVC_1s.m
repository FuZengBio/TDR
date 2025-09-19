clear;clc
FilePath='Z:\Users\TDR\VIP_MSTd_2s & PIVC_1s\Data\';
structures=[1 2 3]; %1=VIP; 2=MSTd; 3=PIVC;
structure_names = {'VIP', 'MSTd', 'PIVC'};
Kernel_method=1;

%%% change order
reverse = 0; %0; 1=reverse parameters, %2=move acceleration to the end

for s = 1:length(structures)
    matfile = fullfile(FilePath, [structure_names{s} '_Data_REG.mat']);
    load(matfile);

    % gaussian kernel
    sigma=30;
    half_guas=sigma*3;
    dist=(-(half_guas):(half_guas))';
    window=1/(sigma*sqrt(2*pi))*exp(-((dist.^2)/(2*sigma^2)))*1000;

    stim_type=[1,2];
    clear tempdata temp_mean Test_PCA_Data 
    for ss=1:length(stim_type)
        for nCell=1:size(REG_raster_per_trial,1)
            clear temp_spike;temp_spike=REG_raster_per_trial{nCell,ss};
            clear sdf;

            for i=1:size(temp_spike,2)
                clear spike; spike=double(squeeze(temp_spike(:,i,:)));
                [ntrials,nsamples]		= size(spike);
                for ii = 1:ntrials
                    convspike	= conv(spike(ii,:),window,'valid');
                    tempdata=downsample(convspike,10);
                    sdf((i-1)*ntrials+ii,:)	= tempdata;
                end
            end
            Test_REG_r{nCell,1}=sdf;
        end

        % PCA parameters
        span=10; %ms - steps to use for the regression and TDR
        if strcmp(structure_names{s}, 'VIP') || strcmp(structure_names{s}, 'MSTd')
            t1 = 101:2000;
            sti_dur=2000;
            [gausacc1, gausvel1, gauspos1]=genGaussian(0.13,6,2);
        else
            t1 = 101:1000;  % PIVC
            sti_dur=1000;
            [gausacc1, gausvel1, gauspos1]=genGaussian(0.13,4,1);
        end
        REG_t=square_sum(t1,span)/span;
        %% regression params to use for Newsome TDR method (done per stimulus)
        REG_params_stim={'choice','FUNangle','ONES'}; %FUNangle is generic for function angle. typically just the headings
        lin_params=1:(length(REG_params_stim)-1); %exclude ONES for the MEAN FR linear fit (since it is added automatically) and for REG_params_kernel below since I subtract out the mean
        REG_choice_dim=find(strcmp(REG_params_stim,'choice')); %find which dim has the choice parameter
        REG_head_dim=find(strcmp(REG_params_stim,'FUNangle')); %find which dim has the heading parameter
        %% regression params to use for Kernel method
        REG_params_kernel={'VELangle','ACCangle','choice'};
        REG_vel_dimK=find(strcmp(REG_params_kernel,'VELangle'));
        REG_acc_dimK=find(strcmp(REG_params_kernel,'ACCangle'));
        REG_choice_dimK=find(strcmp(REG_params_kernel,'choice'));
        %% create ideal kernel
        kernel_del=0;

        REG_acc_ideal1=(square_sum(gausacc1((t1)-kernel_del)',span)/span); REG_acc_ideal1=REG_acc_ideal1./sqrt(REG_acc_ideal1 * REG_acc_ideal1');
        REG_vel_ideal1=(square_sum(gausvel1((t1)-kernel_del)',span)/span); REG_vel_ideal1=REG_vel_ideal1./sqrt(REG_vel_ideal1 * REG_vel_ideal1');
        REG_pos_ideal1=(square_sum(gauspos1((t1)-kernel_del)',span)/span); REG_pos_ideal1=REG_pos_ideal1./sqrt(REG_pos_ideal1 * REG_pos_ideal1');
        REG_one_ideal1=ones(size(REG_vel_ideal1)); REG_one_ideal1=REG_one_ideal1./sqrt(REG_one_ideal1 * REG_one_ideal1');
        %% generate kernel
        kernel=[];
        kernel(REG_choice_dimK,:)=REG_pos_ideal1; %can replace this with other profiles
        kernel(REG_vel_dimK,:)=REG_vel_ideal1; %first head dim used for velocity profile
        kernel(REG_acc_dimK,:)=REG_acc_ideal1; %since I subtract the mean, I don't need the ONES
        %% raster regression fit (kernel method - fit all T together).
        for Ni=1:size(Test_REG_r,1)
            Ni;
            TMP_FVesK=[]; TMP=[];
            clear tempdata; tempdata = REG_Stim{Ni};
            if ~isempty(tempdata)
                clear TMP_Ves;TMP_Ves=(tempdata==1);
                if sum(TMP_Ves)>0
                    for b=1:length(REG_params_kernel) %NB different params to REG_params_stim
                        if ss==1
                            eval(sprintf('TMP(b,:)=REG_%s{Ni}(find(TMP_Ves));',REG_params_kernel{b}))
                        else
                            clear visGroup;visGroup=Ni+size(REG_Stim,1);
                            eval(sprintf('TMP(b,:)=REG_%s{visGroup}(find(TMP_Ves));',REG_params_kernel{b}))
                        end
                    end
                    for trl=1:sum(TMP_Ves)
                        if strcmp(structure_names{s}, 'VIP') || strcmp(structure_names{s}, 'MSTd')
                            TMP_FVesK=[TMP_FVesK repmat(TMP(:,trl),1,sum(REG_t<=2000)) .* kernel];
                        else
                            TMP_FVesK=[TMP_FVesK repmat(TMP(:,trl),1,sum(REG_t<=1000)) .* kernel];
                        end
                    end
                    clear temp_r;
                    if strcmp(structure_names{s}, 'VIP') || strcmp(structure_names{s}, 'MSTd')
                        temp_r=[Test_REG_r{Ni}(TMP_Ves,1:sum(REG_t<=2000))];
                    else
                        temp_r=[Test_REG_r{Ni}(TMP_Ves,1:sum(REG_t<=1000))];
                    end
                    r=bsxfun(@minus,temp_r,mean(temp_r));
                    TMP_bK=((TMP_FVesK*(TMP_FVesK'))^-1)*TMP_FVesK*reshape(r',[],1);
                    for b=1:length(REG_params_kernel)
                        REG_bVesK_Test{b}(Ni,1)=TMP_bK(b);
                    end
                else
                    REG_bVesK_Test{b}(Ni,1)=NaN;
                end
            else
                for b=1:length(REG_params_kernel)
                    REG_bVesK_Test{b}(Ni,1)=NaN; %only one value since done for whole timecourse at once
                end
            end
        end

        if reverse==0
            id=[REG_vel_dimK REG_acc_dimK REG_choice_dimK];
        elseif reverse ==1
            id=[REG_choice_dimK REG_acc_dimK REG_vel_dimK]; %reverse order
        elseif reverse ==2
            id=[REG_vel_dimK REG_choice_dimK REG_acc_dimK]; %move acceleration to the end
        end

        REG_b_use=REG_bVesK_Test(id);

        for nCell=1:size(REG_raster_per_trial,1)
            clear temp_spike;temp_spike=REG_raster_per_trial{nCell,ss};
            clear sdf;

            for i=1:size(temp_spike,1)
                clear spike; spike=double(squeeze(temp_spike(i,:,:)));
                [ntrials,nsamples]		= size(spike);
                                for ii = 1:ntrials
                    convspike	= conv(spike(ii,:),window,'valid');
                    clear tempdata; tempdata=downsample(convspike,10);
                    if strcmp(structure_names{s}, 'VIP') || strcmp(structure_names{s}, 'MSTd')
                        temp_mean(ii,:)=tempdata(1:190);
                    else
                        temp_mean(ii,:)=tempdata(1:90);
                    end
                end
                Test_PCA_Data(1,i).A(:,nCell) = [mean(temp_mean)]';
                Test_PCA_Data(1,i).times =REG_t;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use these parameters to generate the figures.
        clear DataRun; DataRun=Test_PCA_Data;
        clear PCA_params;

        PCA_params.B=REG_b_use;
        PCA_params.normalize=false;%already Zscored, so don't normalize
        PCA_params.numPCs = 12;  % default 6
        times=[];

        [projection,Summary]=TDR_time(DataRun, times, PCA_params);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Use these parameters to generate the figures.
        Flip_plot=0;
        plot_times = [];
        f(1)=figure; set(gcf,'Units','normalized','color','w','position',[0 0.05 1 0.85]);
        plotParams.reusePlot=1;
        plotParams.substPCs = 2;      %1=use PC projections (rather than jPC projections) 2=use the targeted PCA projections
        plotParams.planes2plot=1;
        % plotParams.useAxes=false;
        %%%%%%%%%%%%%%%%%%%%%%%
        REG_params_use=REG_params_kernel(id);
        plotParams.RegNames=REG_params_use;
        plotParams.plotPlanEllipse=false;
        plotParams.planMarkerSize=3;
        plotParams.lineWidth=2;
        plotParams.times=plot_times;
        plotParams.params2plot =[2 1];

        graded=1; %color R/L on a gradual basis
        cuecolor=1; %keep vestibular blue and visual red
        if ss==1
            if strcmp(structure_names{s}, 'VIP') || strcmp(structure_names{s}, 'MSTd')
                condsRun=1:9;
            else
                condsRun=1:10;
            end
            stimtype='Vest';
        elseif ss==2
            if strcmp(structure_names{s}, 'VIP') || strcmp(structure_names{s}, 'MSTd')
                condsRun=13:21;
            else
                condsRun=13:22;
            end
            stimtype='Vis';
        end

        plotParams.colors=plot_colors(condsRun,graded,cuecolor);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot the variance
        fontsize=12;
        clear legR legT
        tempAngle=[-16, -8, -4, -2, -1, 1, 2, 4, 8, 16];
        % tempAngle=[-16, -8, -2, 2, 8, 16];
        for i=1:length(condsRun) %for legend
            legR{i}=sprintf('%g%s',round(10*(tempAngle(i)))/10,'{\circ}');
        end
        for i=1:length(condsRun) %for legend
            legT{i}=sprintf('%g%s',round(10*(tempAngle(i)))/10,'{\circ}');
        end

        if length(REG_params_use)>3
            b_comb=[1 2; 2 3; 2 4; 3 4];
        else
            b_comb = nchoosek(1:length(REG_params_use),2); %combinations of params
        end
        cols=max(length(b_comb),length(REG_params_use) +1);
        rows=2;
        enl=0;% enlarge axes by...
        for comb=1:size(b_comb,1)
            b(1)=b_comb(comb,1); b(2)=b_comb(comb,2);
            if Flip_plot %which regression parameters to plot (NB can switch x and y here!)
                plotParams.params2plot=[b(2) b(1)];
            else
                plotParams.params2plot=[b(1) b(2)];
            end

            subplot(rows,cols,comb); hold on; %plot BASELINE block
            POS=get(gca,'position'); set(gca,'position',[POS(1)-enl/2-0.05 POS(2)-enl/2 POS(3)+enl POS(4)+enl])
            plotParams.colors=plot_colors(condsRun,graded,cuecolor);
            set(gca,'visible','off','fontsize',fontsize)

            %[colorStruct, haxP, vaxP, phR] = phaseSpace(ProjectionRun, SummaryRun, plotParams);  % makes the plot
            [colorStruct, haxP, vaxP, phR] = phaseSpace_AC(projection, Summary, plotParams);
            if comb==size(b_comb,1)
                l=legend(phR(end:-1:1),legR{end:-1:1},'location','east'); %reorder to get -ve at the bottom
                POS=get(l,'position'); set(l,'position',[POS(1)+0.05 POS(2) POS(3) POS(4)])
            end
            set(findobj('type','line'),'clipping','off')
            axis tight
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot Choice, Acc, Vel vs. time
        REG_params_names=nameREG(REG_params_use,1);
        for par=1:length(REG_params_use)
            subplot(rows,cols,cols + par); hold on;
            POS=get(gca,'position'); set(gca,'position',[POS(1)-0.05 POS(2) POS(3) POS(4)])
            set(gca,'fontsize',fontsize)
            for c=1:length(condsRun)
                STPP{c,par}=smooth(projection(c).TargetedPCAproj(:,par),30);%% SmoothTargetedPCAproj
                plot(projection(c).times,STPP{c,par},'-','color',plotParams.colors{c},'linewidth',plotParams.lineWidth);
                m(c)=max(abs(projection(c).TargetedPCAproj(:,par))); %maximum value on plot
            end
            ylabel(sprintf('%s [A.U.]',REG_params_names{par}),'fontsize',fontsize)
            xlabel('time [ms]')
            %     lim(par)=25*ceil(max(m)/25);
            lim=[100 100 125]; %just dicate the lims
            %     ylim([-lim(par) lim(par)]); set(gca,'ytick',[-lim(par) 0 lim(par)])
            %     daspect([1000 300 1])
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 3D plot on top row and then variance explained on second row
        subplot(rows,cols,length(REG_params_use) +1); hold on;
        if length(REG_params_use)>=3
            iv=find(strcmp(REG_params_use,'VELangle')); ia=find(strcmp(REG_params_use,'ACCangle')); ic=find(strcmp(REG_params_use,'choice'));
            limva=max([lim(iv) lim(ia)]);
            [xx,yy] = meshgrid([-limva limva], [-limva limva]);
            surf(xx,yy,zeros(size(xx)))
            colormap gray
            shading flat
            alpha 0.7
            %12
            plot3(projection(1).TargetedPCAproj(:,iv),projection(1).TargetedPCAproj(:,ia),projection(1).TargetedPCAproj(:,ic),'linewidth',3,'color',plotParams.colors{1},'clipping','off')
            stem3(projection(1).TargetedPCAproj(:,iv),projection(1).TargetedPCAproj(:,ia),projection(1).TargetedPCAproj(:,ic),'color',plotParams.colors{1},'MarkerSize',1,'clipping','off')
            %patch(projection(1).TargetedPCAproj(:,iv),projection(1).TargetedPCAproj(:,ia),zeros(size(projection(1).TargetedPCAproj(:,ic))),'FaceAlpha',0.5,'FaceColor',[0 0 0])
            plot3(projection(1).TargetedPCAproj(:,iv),projection(1).TargetedPCAproj(:,ia),zeros(size(projection(1).TargetedPCAproj(:,ic))),'color',[1 1 1]*0.5)
            %-12
            plot3(projection(length(condsRun)).TargetedPCAproj(:,iv),projection(length(condsRun)).TargetedPCAproj(:,ia),projection(length(condsRun)).TargetedPCAproj(:,ic),'linewidth',3,'color',plotParams.colors{length(condsRun)},'clipping','off')
            stem3(projection(length(condsRun)).TargetedPCAproj(:,iv),projection(length(condsRun)).TargetedPCAproj(:,ia),projection(length(condsRun)).TargetedPCAproj(:,ic),'color',plotParams.colors{length(condsRun)},'MarkerSize',1,'clipping','off')
            %patch(projection(length(condsRun)).TargetedPCAproj(:,iv),projection(length(condsRun)).TargetedPCAproj(:,ia),zeros(size(projection(length(condsRun)).TargetedPCAproj(:,ic))),'FaceAlpha',0,'FaceColor',[0 0 0])
            plot3(projection(length(condsRun)).TargetedPCAproj(:,iv),projection(length(condsRun)).TargetedPCAproj(:,ia),zeros(size(projection(length(condsRun)).TargetedPCAproj(:,ic))),'color',[1 1 1]*0.5)
            view(-65,45)
            grid on
            xlabel(sprintf('%s [A.U.]',REG_params_names{iv}),'fontsize',fontsize)
            ylabel(sprintf('%s [A.U.]',REG_params_names{ia}),'fontsize',fontsize)
            zlabel(sprintf('%s [A.U.]',REG_params_names{ic}),'fontsize',fontsize)

            xlim([-limva limva]); ylim([-limva limva]); zlim([-lim(ic) lim(ic)]); set(gca,'xtick',[-limva 0 limva],'ytick',[-limva 0 limva],'ztick',[-lim(ic) 0 lim(ic)])
            daspect([1 1 1])
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plot the variance vs. time
        subplot(rows,cols,cols + length(REG_params_use) +1); hold on;
        POS=get(gca,'position'); set(gca,'position',[POS(1)-enl/2-0.05 POS(2)-enl/2 POS(3)+enl POS(4)+enl])
        set(gca,'LineStyleOrder',{'-','--',':','-.'},'ColorOrder',[0 0 0],'color','none','fontsize',fontsize); hold on;
        clear I SV
        for I=1:3
            SV(:,I) = smooth(Summary.varCaptEachTargetedPCT(:,I),30);%% SmoothVariance
        end
        plot(REG_t,SV,'linewidth',plotParams.lineWidth)

        ylabel('Variance Explained')
        xlabel('time [ms]')
        l=legend(REG_params_names);
        POS=get(l,'position');
        set(l,'position',[POS(1)+0.05 POS(2) POS(3) POS(4)])
        yl=ylim;
        for i=1:size(b_comb,1)
            text(sti_dur,(size(b_comb,1)-i+1)*(yl(2)-yl(1))/(2*size(b_comb,1)),sprintf('{\\angle}(%s,%s)=%u{\\circ}',REG_params_names{b_comb(i,1)},REG_params_names{b_comb(i,2)},round(Summary.Btheta(b_comb(i,1),b_comb(i,2)))))
        end

        OutputPath = 'Z:\Users\TDR\VIP_MSTd_2s & PIVC_1s\results\';
        set(gcf,'paperpositionmode','auto');
        filename = fullfile(OutputPath, sprintf('%s_%s_Kernel%u_%s.tif', structure_names{s}, stimtype, Kernel_method, strjoin(REG_params_use, '_')));
        saveas(gcf, filename);
        close;
    end
end

