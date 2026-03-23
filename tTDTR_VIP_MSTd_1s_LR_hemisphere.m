initialize_monkey
load(strcat(RESULTS_DIR,'PCAData'))
load(strcat(RESULTS_DIR,'RegData'))
load(strcat(RESULTS_DIR,'Cells'));
load(strcat(RESULTS_DIR,'LR_hemisphere'));  % LR_hemisphere: N_CELL x 1, 1=Left,2=Right

N_CELL=length(CELL_monk);
plot_times = [];
times=[];

MOVIE=0;
Zscore=0;
graded=1;
cuecolor=1;
Kernel_method=1;
fontsize=12;
reverse = 0;%0:VAC; 1=reverse parameters, %2=move acceleration to the end

TYPES=[1 2]; %1=Ves,2=Vis
structures=[2 3]; %2=VIP; 3=MSTd;
trial_durs=1000;

TYPES_name={'Ves' 'Vis' ''};
STRUCTURE_name={'' 'VIP' 'MSTd'};
% ---- significance params 
nPerm    = 1000;
alphaPt  = 0.05;
alphaCl  = 0.05;

for TYPE=TYPES
    for trial_dur=trial_durs
        for structure=structures

            if trial_dur==1000
                monk=[19 24];
                blocks = [1 3];
            else
                if structure==2
                    structure=7; monk=[5 14];
                elseif structure==3
                    structure=8; monk=[2 4];
                else
                    break
                end
                blocks=6;
            end

            blockRun=blocks(1);
            condsRun=find(abs(DataDir{monk(1),structure})<head_cut);

            crit=1;
            SigShift=~isnan(CELL_shiftVes(:,Nseg)) | ~isnan(CELL_shiftVis(:,Nseg)); 

            if TYPE==3
                VesVisOnly=3; Flip_plot=0;
            elseif TYPE==1 || TYPE==2
                VesVisOnly=TYPE; Flip_plot=1;
            else
                VesVisOnly=0; Flip_plot=0;
            end

            if VesVisOnly==1
                condsRun=condsRun(condsRun<=11);
            elseif VesVisOnly==2
                condsRun=condsRun((condsRun<=22) & (condsRun>11));
            elseif VesVisOnly==3 || ismember(structure,4)
                condsRun=condsRun((condsRun<=22));
            end

            for block=blocks
                Dmean{block}=0;
                Dstd{block}=1;
            end

            PCA_params.normalize=false;
            mo=ismember(CELL_monk,monk)';

            for block=blocks
                st=CELL_rec_structure(:,block)==structure;
                idxAll{block}= st(:) & mo(:) & crit(:);   
                for con=condsRun
                    Data{structure,block}(con).A=(Data_PCA(con,block).A(REG_t<=trial_dur,:)-Dmean{block})./Dstd{block};
                    Data{structure,block}(con).A(:,~idxAll{block})=NaN;
                    Data{structure,block}(con).times=Data_PCA(con,block).times(REG_t<=trial_dur);
                end
            end

            % RUN PC on this data
            DataRun=Data{structure,blockRun}(condsRun);

            % remove cells not present across all conds
            no_cell=zeros(1,N_CELL);
            for i=1:length(DataRun)
                no_cell=no_cell+any(isnan(DataRun(i).A));
            end
            for i=1:length(DataRun)
                DataRun(i).A=DataRun(i).A(:,~no_cell);
            end

            % N = number of kept cells
            N=sum(~no_cell);
            fprintf('\nTYPE=%d, struct=%d, dur=%d: N=%u cells used\n', TYPE, structure, trial_dur, N)

            % ------------------------------------------------------------
            % Build id / Brun
            % ------------------------------------------------------------
            if Kernel_method
                if reverse==0
                    id=[REG_RL_dimK REG_vel_dimK REG_acc_dimK REG_choice_dimK];
                elseif reverse==1
                    id=[REG_choice_dimK REG_acc_dimK REG_vel_dimK];
                elseif reverse==2
                    id=[REG_RL_dimK REG_vel_dimK REG_choice_dimK REG_acc_dimK];
                end
            else
                if ~reverse
                    id=[REG_head_dim REG_choice_dim];
                else
                    id=[REG_choice_dim REG_head_dim];
                end
            end

            id = id(~isnan(id) & id>0);

            if VesVisOnly==1
                if Kernel_method
                    REG_b_use=REG_bVesK(id,:);
                    REG_params_use=REG_params_kernel(id);
                else
                    REG_b_use=REG_bVes(id,:);
                    REG_params_use=REG_params_stim(id);
                end
            elseif VesVisOnly==2
                if Kernel_method
                    REG_b_use=REG_bVisK(id,:);
                    REG_params_use=REG_params_kernel(id);
                else
                    REG_b_use=REG_bVis(id,:);
                    REG_params_use=REG_params_stim(id);
                end
            else
                id=[1 2];
                if Kernel_method
                    REG_b_use=REG_bK(id,:);
                else
                    REG_b_use=REG_b(id,:);
                end
                REG_params_use=REG_params(id);
            end

            Brun=[];
            for b=1:size(REG_b_use,1)
                Brun{b}=REG_b_use{b,blockRun}(~no_cell,:);  % N x 1
            end
            PCA_params.numPCs = 12;
            PCA_params.B=Brun;

            % RUN (optional)
            [ProjectionRun, SummaryRun] = TDR_time(DataRun, times, PCA_params); 

            kept_global_idx = find(~no_cell);          % size N, map back to original cell index
            hemi_kept = LR_hemisphere(kept_global_idx); % N x 1

            idxL = find(hemi_kept==1);
            idxR = find(hemi_kept==2);

            if isempty(idxL) || isempty(idxR)
                warning('No both hemispheres after filtering; skip significance.');
                continue;
            end

            [V_L, Lproj, tvec] = computeVar_subsetVAC(DataRun, PCA_params, idxL);
            [V_R, Rproj,~]    = computeVar_subsetVAC(DataRun, PCA_params, idxR);
            D_real = V_L - V_R;  % T x nAxes

            Tlen = size(V_L,1);
            nAxes = size(V_L,2);

            nL = numel(idxL);
            allIdx = (1:N)';

            D_null = nan(Tlen, nPerm, nAxes);
            for p=1:nPerm
                rp = randperm(N);
                idxLp = allIdx(rp(1:nL));
                idxRp = allIdx(rp(nL+1:end));

                V_Lp = computeVar_subsetVAC(DataRun, PCA_params, idxLp);
                V_Rp = computeVar_subsetVAC(DataRun, PCA_params, idxRp);

                D_null(:,p,:) = V_Lp - V_Rp;
            end

            axisNamesOrdered = REG_params_use; 
            if ischar(axisNamesOrdered), axisNamesOrdered = cellstr(axisNamesOrdered); end
            sig = struct();
            for k=1:nAxes
                axName = axisNamesOrdered{k};
                [mask, thrPoint, clThr] = clusterPerm_twoSided(D_real(:,k), squeeze(D_null(:,:,k)), alphaPt, alphaCl);
                sig.(axName).mask = mask;
                sig.(axName).thrPoint = thrPoint;
                sig.(axName).clThr = clThr;
            end
            figure('Color','w','Name',sprintf('%s struct=%d dur=%d hemi variance', TYPES_name{TYPE}, structure, trial_dur));
            for k=1:nAxes
                axName = axisNamesOrdered{k};

                subplot(nAxes,1,k); hold on; box on
                plot(tvec, V_L(:,k), 'k-',  'LineWidth',2);
                plot(tvec, V_R(:,k), 'k--', 'LineWidth',1.5);

                yl = ylim;
                m = sig.(axName).mask;
                if any(m)
                    ymark = yl(1) + 0.03*diff(yl);
                    plot(tvec(m), ymark*ones(sum(m),1), '.', 'MarkerSize',10);
                end

                title(sprintf('%s: Left(solid) vs Right(dashed)', axName));
                ylabel('Var explained');
                if k==nAxes, xlabel('Time'); end
                legend({'Left','Right','Sig (cluster)'},'Location','best');
            end


            clear OutHemi;
            % save
            OutHemi.order=REG_params_use;
            OutHemi.Lproj=Lproj;
            OutHemi.Rproj=Rproj;
            OutHemi.V_left = V_L;
            OutHemi.V_right= V_R;
            OutHemi.sig = sig;

            OutputPath='Z:\Users\TDR\results\tTDR_LR_hem\';
            SaveFileName=[OutputPath STRUCTURE_name{structure} '_', TYPES_name{TYPE}, '_', cell2mat(REG_params_use)];
            save(SaveFileName, 'OutHemi','-v7.3');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Step10: save the figures
            set(gcf,'paperpositionmode','auto');
            filename=[OutputPath STRUCTURE_name{structure} '_',TYPES_name{TYPE}, '_', cell2mat(REG_params_use)];
            saveas(gcf,filename,'png');
            close;
        end
    end
end


% ============================================================
% HELPERS (copy-paste to the end of your script)
% ============================================================
function [Vsub, TargetedPCAproj, tvec] = computeVar_subsetVAC(DataRun, PCA_params, idxKeep)
% subset neuron columns in DataRun, subset rows in PCA_params.B axes
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
% two-sided cluster permutation on |curve|
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
