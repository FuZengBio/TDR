function [ Out ] = nameREG( RegNames,markortho )

if nargin<2
    markortho=1;
end
for i=1:length(RegNames)
    if strcmp(RegNames{i},'X') || strcmp(RegNames{i},'FUNangle')
        Out{i}='Heading';
    elseif  strcmp(RegNames{i},'RL')
        Out{i}='Heading Sign';
    elseif  strcmp(RegNames{i},'VELangle')
        Out{i}='Vel.';
    elseif  strcmp(RegNames{i},'ACCangle')
        Out{i}='Accel.';
    elseif  strcmp(RegNames{i},'choice')
        Out{i}='Choice';
    elseif  strcmp(RegNames{i},'Stim')
        Out{i}='Cue';
    end
    if i>1 & markortho
        Out{i}=sprintf('%s %s',Out{i},'(Ortho)');
    end
end

