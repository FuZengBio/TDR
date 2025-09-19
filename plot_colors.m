function [ colors ] = plot_colors( conds,graded,cuecolor )

if (any(ismember(conds,[1:11])) && any(ismember(conds,[12:22]))) || cuecolor %if there is a visual and vestibular condition, use cue-colors blue and red
    for i=1:length(conds)  %define colors
        cond=conds(i);
        if ismember(cond,[1:11]) %vestibular
            if length(conds)==4 & graded, colors{i}=([0 0 1] + [0 1 0]*(i)/3) *(i)/3 ;  %0 degree headings (testing)
            %elseif graded, colors{i}=([0 0 1] + [0 1 0]*(cond)/11) *(cond)/11 ; %blue to cyan
            elseif graded, colors{i}=([0 0 1] + [0 1 0]*(12-(cond))/11) *(12-(cond))/11 ; %blue to cyan (flip color order)
            else colors{i}=([0 0 1]) *0.7; end %blue / cyan
        elseif ismember(cond,[12:22]) %visual
            if length(conds)==4 & graded, co
                lors{i}=([1 0 0] + [0 0 1]*((i)-2)/3) *((i)-2)/3;  %0 degree headings (testing)
            %elseif graded, colors{i}=([1 0 0] + [0 0 1]*(cond-11)/11) *(cond-11)/11; %red to magenta
            elseif graded, colors{i}=([1 0 0] + [0 0 1]*(12-(cond-11))/11) *(12-(cond-11))/11; %red to magenta (flip color order)
            else colors{i}=([1 0 0]) *0.7; end %red / magenta
        elseif ismember(cond,[23:33])  %combined
            if length(conds)==4 & graded, colors{i}=([0 1 0] + [1 0 0]*((i)-4)/3) *((i)-4)/3;  %0 degree headings (testing)
            %elseif graded, colors{i}=([0 1 0] + [1 0 0]*(cond-22)/11) *(cond-22)/11; %green to yellow
            elseif graded, colors{i}=([0 1 0] + [1 0 0]*(12-(cond-22))/11) *(12-(cond-22))/11; %green to yellow (flip color order)
            else colors{i}=([0 1 0]) *0.7; end %green / yellow
        else
            colors{i}=[0 0 0];
        end
    end
else %define colors by heading only (not stim)
    %     Nc=4; %how extreme I want the colors to be
    %     conds=mod(conds-1,11)+1; %heading
    %     scat_col=cool(Nc*length(conds)-(Nc-1));
    %     for i=1:length(conds)  %define colors
    %         colors{i}=scat_col(1+Nc*(i-1),:);
    %     end
    
    if length(conds)==7
        colors{1}=[1 0.5 0.5]*0.5;
        colors{2}=[1 0.5 0.5]*0.8;
        colors{3}=[1 0 0];
        colors{4}=[0 1 0]*0.8;
        colors{5}=[0 0 1];
        colors{6}=[0.5 0.5 1]*0.8;
        colors{7}=[0.5 0.5 1]*0.5;
    end
    
    if length(conds)==9 %aihua's data
        colors{1}=[1 0.5 0.5]*0.3;
        colors{2}=[1 0.5 0.5]*0.6;
        colors{3}=[1 0.5 0.5]*0.9;
        colors{4}=[1 0 0];
        colors{5}=[0 1 0]*0.8;
        colors{6}=[0 0 1];
        colors{7}=[0.5 0.5 1]*0.9;
        colors{8}=[0.5 0.5 1]*0.6;
        colors{9}=[0.5 0.5 1]*0.3;
    end
    
    if length(conds)==3
        colors{1}=[1 0 0];
        colors{2}=[0 1 0]*0.8;
        colors{3}=[0 0 1];
    end
end