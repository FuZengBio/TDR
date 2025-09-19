span=10; %ms - steps to use for the regression and TDR
t1=[101:1000]; %times [ms] to use for regression and PCA, downsampled by 'span' below
REG_t=square_sum(t1,span)/span;

%regression params to use for Newsome TDR method (done per stimulus)
REG_params_stim={'choice','FUNangle','ONES'}; %FUNangle is generic for function angle. typically just the headings
lin_params=1:(length(REG_params_stim)-1); %exclude ONES for the MEAN FR linear fit (since it is added automatically) and for REG_params_kernel below since I subtract out the mean
REG_choice_dim=find(strcmp(REG_params_stim,'choice')); %find which dim has the choice parameter
REG_head_dim=find(strcmp(REG_params_stim,'FUNangle')); %find which dim has the heading parameter

% regression params to use for Kernel method
REG_params_kernel={'choice','VELangle','ACCangle'};
REG_vel_dimK=find(strcmp(REG_params_kernel,'VELangle'));
REG_acc_dimK=find(strcmp(REG_params_kernel,'ACCangle'));
REG_choice_dimK=find(strcmp(REG_params_kernel,'choice'));

%create ideal kernel
kernel_del=0; %(between 0 and 100) delay the kernel by kernel_del (ms) estimated delay between stim and phsyiology?
[gausacc1, gausvel1, gauspos1]=genGaussian(0.13,4,1); %using the function provided from Jing and Jian that calculates the actual trajectories. Values corerspond to t=REG_t
REG_acc_ideal1=(square_sum(gausacc1((t1)-kernel_del)',span)/span); REG_acc_ideal1=REG_acc_ideal1./sqrt(REG_acc_ideal1 * REG_acc_ideal1');
REG_vel_ideal1=(square_sum(gausvel1((t1)-kernel_del)',span)/span); REG_vel_ideal1=REG_vel_ideal1./sqrt(REG_vel_ideal1 * REG_vel_ideal1');
REG_pos_ideal1=(square_sum(gauspos1((t1)-kernel_del)',span)/span); REG_pos_ideal1=REG_pos_ideal1./sqrt(REG_pos_ideal1 * REG_pos_ideal1');
REG_one_ideal1=ones(size(REG_vel_ideal1)); REG_one_ideal1=REG_one_ideal1./sqrt(REG_one_ideal1 * REG_one_ideal1');
% can use choice kernal directly from the data (CURRENTLY SPECIFIC TO VIP! not included here)

%generate kernel
kernel=[];
kernel(REG_choice_dimK,:)=REG_pos_ideal1; %can replace this with other profiles
kernel(REG_vel_dimK,:)=REG_vel_ideal1; %first head dim used for velocity profile
kernel(REG_acc_dimK,:)=REG_acc_ideal1; %since I subtract the mean, I don't need the ONES

%load the data and trial parameters per cell here
load REG_r; %REG_r{Ni}= ; %load the sdf per cell here for cell Ni. size: number_of_trials X time (90 if duration=1000, since it is in steps of 10ms (span) and starting at 100 until 1000s)
load REG_choice; %REG_choice{Ni} = ; %load the choice per trial here
load REG_FUNangle; %REG_FUNangle{Ni} = ; %load the heading per trial here
load REG_ACCangle; %REG_ACCangle{Ni} = ; %load the heading per trial here
load REG_VELangle; %REG_VELangle{Ni} = ; %load the heading per trial here
load REG_ONES; %REG_ONES{Ni} = ; %ones the size of # of trials
load REG_Stim; %REG_Stim{Ni} = ; %1 for ves, -1 for vis (currently only works for unisensory)
%TMP_Vis = visual trials;
%TMP_Ves = vestibular trials;

%build F matrix from the trial parameters
Ni=1
blockrun=1
b=1

TMP_FVis=[];
for b=1:length(REG_params_stim)
    eval(sprintf('TMP_FVis=[TMP_FVis REG_%s{Ni}(find(TMP_Vis))];',REG_params_stim{b}))
end
TMP_FVisT=TMP_FVis'; %transpose
TMP_b=((TMP_FVisT*(TMP_FVisT'))^-1)*TMP_FVisT*REG_r{Ni}(find(TMP_Vis),:);
for b=1:length(REG_params_stim)
    REG_bVis{b}(Ni,REG_t<=1000)=TMP_b(b,:);
end

%raster regression fit (kernel method - fit all T together).
%NB: I subtract the mean and thus don't need a ONES term
TMP_FVisK=[]; TMP=[];
for b=1:length(REG_params_kernel) %NB different params to REG_params_stim
    clear tempdata; tempdata = REG_choice{Ni,blockrun};
    clear TMP_Vis;TMP_Vis=(tempdata==-1);
    eval(sprintf('TMP(b,:)=REG_%s{Ni,blockrun}(find(TMP_Vis));',REG_params_kernel{b}))
end
for trl=1:sum(TMP_Vis)
    TMP_FVisK=[TMP_FVisK repmat(TMP(:,trl),1,sum(REG_t<=1000)) .* kernel];
end
r=bsxfun(@minus,REG_r{Ni}(find(TMP_Vis),:),mean(REG_r{Ni}(find(TMP_Vis),:))); %subtract mean sdf
TMP_bK=((TMP_FVisK*(TMP_FVisK'))^-1)*TMP_FVisK*reshape(r',[],1);
for b=1:length(REG_params_kernel)
    REG_bVisK{b,blockrun}(Ni,1)=TMP_bK(b); %only one value since done for whole timecourse at once
end

