clear;clc
WORK_DIR =  'Z:\Users\TDR\VIP_MSTd_1s\';
DATA_DIR =  'Z:\Users\TDR\VIP_MSTd_1s\Data\';
RESULTS_DIR = 'Z:\Users\TDR\VIP_MSTd_1s\results\';
load(strcat(DATA_DIR,'PCAData'))
load(strcat(DATA_DIR,'RegData'))
load(strcat(DATA_DIR,'Cells'));
N_CELL=length(CELL_monk);
plot_times = [];
times=[];
graded=1;
fontsize=12;
cuecolor=1; %keep vestibular blue and visual red

Kernel_method=0; %1=use my kernel (0= Newsome method)
reverse = 0; %0; 1=reverse parameters, %2=move acceleration to the end
TYPES=[1 2]; %1: vestibular; 2: visual;
structures=[2 3]; %2=VIP; 3=MSTd;
STRUCTURE_name={'' 'VIP' 'MSTd'};
trial_durs=1000;

TYPES_name={'Ves' 'Vis'};
for TYPE=TYPES
    for trial_dur=trial_durs %1s or 2s
        for structure=structures
            if trial_dur==1000
                monk=[19 24];
                blocks = 1;
            else
                if structure==2
                    structure=7; %7=VIP (2s)
                    monk=[5 14]; %Aihua's 2s VIP2 data;
                elseif structure==3
                    structure=8; %8=MST (2s)
                    monk=[2 4]; %Yong's 2s MST2 data
                end
            end
            global DataDir;
            DataDir{19,2}=DataDir{24,2}; %Alvin VIP
            DataDir{19,3}=DataDir{24,2}; %Alvin MST
            blockRun=blocks;
            head_cut=25;
            condsRun=find(abs(DataDir{monk(1),structure})<head_cut);
            crit=1;
            VesVisOnly=TYPE; %vestibular or visual
            Flip_plot=1; %flip x and y to consistently keep choice on x axis
            if VesVisOnly==1 %only do Vest conditions
                condsRun=condsRun(condsRun<=11);
            elseif VesVisOnly==2 %only do Vis conditions
                condsRun=condsRun((condsRun<=22) & (condsRun>11));
            end

            for block=blocks
                Dmean{block}=0;

                Dstd{block}=1;
            end
            PCA_params.normalize=false;
            mo=ismember(CELL_monk,monk)';

            %for struc=structure
            for block=blocks
                st=CELL_rec_structure(:,block)==structure;
                idx{block}= st(:) & mo(:) & crit(:);
                for con=condsRun
                    Data{structure,block}(con).A=(Data_PCA(con,block).A(REG_t<=trial_dur,:)-Dmean{block})./Dstd{block};
                    Data{structure,block}(con).A(:,~idx{block})=NaN;
                    Data{structure,block}(con).times=Data_PCA(con,block).times(REG_t<=trial_dur);
                end
            end

            %RUN PC on this data (average over trials). condsRun defined above
            DataRun=Data{structure,blockRun}(condsRun);

            %remove any cells that weren't present across ALL conditions (in RUN and TEST)
            no_cell=zeros(1,N_CELL);
            for i=1:length(DataRun) %loop through RUN conditions
                no_cell=no_cell+any(isnan(DataRun(i).A));
            end

            for i=1:length(DataRun) %loop through RUN conditions
                DataRun(i).A=DataRun(i).A(:,~no_cell);
            end

            N=sum(~no_cell); %number of cells
            fprintf('N=%u cells used\n',N)

            %calculate regression variables for targeted dimensionality reduction
            if Kernel_method
                if reverse==0
                    id=[REG_RL_dimK REG_vel_dimK REG_acc_dimK REG_choice_dimK];
                elseif reverse ==1
                    id=[REG_choice_dimK REG_acc_dimK REG_vel_dimK]; %reverse order
                elseif reverse ==2
                    id=[REG_RL_dimK REG_vel_dimK REG_choice_dimK REG_acc_dimK]; %move acceleration to the end
                end
            else
                if ~reverse
                    id=[REG_head_dim REG_choice_dim];
                else
                    id=[REG_choice_dim REG_head_dim];
                end
            end
            if VesVisOnly==1 %orthogonalize based on Vest conditions
                if Kernel_method
                    REG_b_use=REG_bVesK(id,:); %single beta (instead of over time. kernel method)
                    REG_params_use=REG_params_kernel(id);
                else
                    REG_b_use=REG_bVes(id,:);
                    REG_params_use=REG_params_stim(id);
                end
            elseif VesVisOnly==2 %orthogonalize based on Vis conditions
                if Kernel_method
                    REG_b_use=REG_bVisK(id,:); %single beta (instead of over time. kernel method)
                    REG_params_use=REG_params_kernel(id);
                else
                    REG_b_use=REG_bVis(id,:);
                    REG_params_use=REG_params_stim(id);
                end
            else %orthogonlize based on default in REG_params or something else specified
                id=[1 2];
                %NB currently I just use ones for kernel - not fully implemented
                if Kernel_method
                    REG_b_use=REG_bK(id,:); %single beta (instead of over time. kernel method)
                else
                    REG_b_use=REG_b(id,:);
                end
                REG_params_use=REG_params(id);
            end
            REG_params_namesO=nameREG(REG_params_use,1);
            REG_params_names=nameREG(REG_params_use,0);
            Brun=[];
            for b=1:size(REG_b_use,1)
                Brun{b}=REG_b_use{b,blockRun}(~no_cell,:);
            end

            %RUN
            PCA_params.numPCs = 12;  % default 6
            PCA_params.B=Brun;
            [ProjectionRun, SummaryRun] = TDR_time(DataRun, times, PCA_params);
            Brun_denoised=SummaryRun.B_denoised;
            varCaptEachPC(structure,TYPE,:)=SummaryRun.varCaptEachPC;
            varCaptEachTargetedPCT{structure,TYPE}=SummaryRun.varCaptEachTargetedPCT; %save for plotting across structures/cues
            varCaptEachTargetedPC(structure,TYPE,:)=SummaryRun.varCaptEachTargetedPC;
            varCaptEachTargetedPC_sum(structure,TYPE) = sum(varCaptEachTargetedPC(structure,TYPE,:));


            %plot population projection onto targeted axes
            f(1)=figure; set(gcf,'Units','normalized','color','w','position',[0 0.05 1 0.85]);
            plotParams.reusePlot=1;
            plotParams.substPCs = 2;      %1=use PC projections (rather than jPC projections) 2=use the targeted PCA projections
            plotParams.planes2plot=1;
            %plotParams.useAxes=false;
            plotParams.RegNames=REG_params_use;
            plotParams.plotPlanEllipse=false;
            plotParams.planMarkerSize=3;
            plotParams.lineWidth=2;
            plotParams.times=plot_times;
            clear legR legT
            for i=1:length(condsRun) %for legend
                legR{i}=sprintf('%g%s',round(10*(DataDir{monk(1),structure}(condsRun(i))))/10,'{\circ}');
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

                %plot BASELINE phasespace plot
                [colorStruct, haxP, vaxP, phR] = phaseSpace(ProjectionRun, SummaryRun, plotParams);  % makes the plot
                if comb==size(b_comb,1)
                    l=legend(phR(end:-1:1),legR{end:-1:1},'location','east'); %reorder to get -ve at the bottom
                    POS=get(l,'position'); set(l,'position',[POS(1)+0.05 POS(2) POS(3) POS(4)])
                end
                set(findobj('type','line'),'clipping','off')
                axis tight

            end

            for par=1:length(REG_params_use)

                subplot(rows,cols,cols + par); hold on;
                POS=get(gca,'position'); set(gca,'position',[POS(1)-enl/2-0.05 POS(2)-enl/2 POS(3)+enl POS(4)+enl])
                set(gca,'fontsize',fontsize)

                for c=1:length(condsRun)
                    plot(ProjectionRun(c).times,ProjectionRun(c).TargetedPCAproj(:,par),'-','color',plotParams.colors{c},'linewidth',plotParams.lineWidth);
                    m(c)=max(abs(ProjectionRun(c).TargetedPCAproj(:,par))); %maximum value on plot
                end

                ylabel(sprintf('%s [A.U.]',REG_params_namesO{par}),'fontsize',fontsize)
                xlabel('time [ms]')
                if ~Kernel_method
                    lim=[150 100];
                    ylim([-lim(par) lim(par)]); set(gca,'ytick',[-lim(par) 0 lim(par)]);  %for kernel==0
                else
                    lim(par)=25*ceil(max(m)/25);
                    if reverse
                        lim=[75 75 100]; %just dicate the lims
                    else
                        lim=[120 60 60]; %just dicate the lims
                    end
                    ylim([-lim(par) lim(par)]); set(gca,'ytick',[-lim(par) 0 lim(par)])
                end
                daspect([1000 300 1])
            end

            %plot 3D plot on top row and then variance explained on second row
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
                plot3(ProjectionRun(1).TargetedPCAproj(:,iv),ProjectionRun(1).TargetedPCAproj(:,ia),ProjectionRun(1).TargetedPCAproj(:,ic),'linewidth',3,'color',plotParams.colors{1},'clipping','off')
                stem3(ProjectionRun(1).TargetedPCAproj(:,iv),ProjectionRun(1).TargetedPCAproj(:,ia),ProjectionRun(1).TargetedPCAproj(:,ic),'color',plotParams.colors{1},'MarkerSize',1,'clipping','off')
                %patch(ProjectionRun(1).TargetedPCAproj(:,iv),ProjectionRun(1).TargetedPCAproj(:,ia),zeros(size(ProjectionRun(1).TargetedPCAproj(:,ic))),'FaceAlpha',0.5,'FaceColor',[0 0 0])
                plot3(ProjectionRun(1).TargetedPCAproj(:,iv),ProjectionRun(1).TargetedPCAproj(:,ia),zeros(size(ProjectionRun(1).TargetedPCAproj(:,ic))),'color',[1 1 1]*0.5)
                %-12
                plot3(ProjectionRun(length(condsRun)).TargetedPCAproj(:,iv),ProjectionRun(length(condsRun)).TargetedPCAproj(:,ia),ProjectionRun(length(condsRun)).TargetedPCAproj(:,ic),'linewidth',3,'color',plotParams.colors{length(condsRun)},'clipping','off')
                stem3(ProjectionRun(length(condsRun)).TargetedPCAproj(:,iv),ProjectionRun(length(condsRun)).TargetedPCAproj(:,ia),ProjectionRun(length(condsRun)).TargetedPCAproj(:,ic),'color',plotParams.colors{length(condsRun)},'MarkerSize',1,'clipping','off')
                %patch(ProjectionRun(length(condsRun)).TargetedPCAproj(:,iv),ProjectionRun(length(condsRun)).TargetedPCAproj(:,ia),zeros(size(ProjectionRun(length(condsRun)).TargetedPCAproj(:,ic))),'FaceAlpha',0,'FaceColor',[0 0 0])
                plot3(ProjectionRun(length(condsRun)).TargetedPCAproj(:,iv),ProjectionRun(length(condsRun)).TargetedPCAproj(:,ia),zeros(size(ProjectionRun(length(condsRun)).TargetedPCAproj(:,ic))),'color',[1 1 1]*0.5)
                view(-65,45)
                grid on
                xlabel(sprintf('%s [A.U.]',REG_params_namesO{iv}),'fontsize',fontsize)
                ylabel(sprintf('%s [A.U.]',REG_params_namesO{ia}),'fontsize',fontsize)
                zlabel(sprintf('%s [A.U.]',REG_params_namesO{ic}),'fontsize',fontsize)

                xlim([-limva limva]); ylim([-limva limva]); zlim([-lim(ic) lim(ic)]); set(gca,'xtick',[-limva 0 limva],'ytick',[-limva 0 limva],'ztick',[-lim(ic) 0 lim(ic)])
                daspect([1 1 1])
            end
            subplot(rows,cols,cols + length(REG_params_use) +1); hold on;
            POS=get(gca,'position'); set(gca,'position',[POS(1)-enl/2-0.05 POS(2)-enl/2 POS(3)+enl POS(4)+enl])
            set(gca,'LineStyleOrder',{'-','--',':','-.'},'ColorOrder',[0 0 0],'color','none','fontsize',fontsize); hold on;
            plot(REG_t(REG_t<=trial_dur),SummaryRun.varCaptEachTargetedPCT,'linewidth',plotParams.lineWidth)
            ylabel('Variance Explained')
            xlabel('time [ms]')
            l=legend(REG_params_namesO);
            POS=get(l,'position');
            set(l,'position',[POS(1)+0.05 POS(2) POS(3) POS(4)])
            yl=ylim;
            for i=1:size(b_comb,1)
                text(trial_dur,(size(b_comb,1)-i+1)*(yl(2)-yl(1))/(2*size(b_comb,1)),sprintf('{\\angle}(%s,%s)=%u{\\circ}',REG_params_names{b_comb(i,1)},REG_params_names{b_comb(i,2)},round(SummaryRun.Btheta(b_comb(i,1),b_comb(i,2)))))
            end
            set(gcf,'paperpositionmode','auto');
            filename = fullfile(RESULTS_DIR, sprintf('PCAspace_%s_%s_Kernel%u_%s.tif', STRUCTURE_name{structure}, TYPES_name{TYPE}, Kernel_method,strjoin(REG_params_use, '_')));
            print(gcf, '-dtiff', '-r300', filename)
            close
        end
    end
end