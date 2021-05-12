%%
clear all;
% close all;
run ../localdef.m

addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_fooof))
addpath(genpath(path_LSCPtools));

files=dir([data_path filesep 'CIcfe_ft_*.mat']);

% table=readtable([save_path 'CTET_behav_res.txt']);
% load('../ICA_Artifacts.mat')


%%
nFc=0;
all_ITI_info=[];
for nF=1:length(files)
    File_Name=files(nF).name;
    File_Name=File_Name(10:end);
    septag=findstr(File_Name,'_');
    SubN=str2num(File_Name(1:septag(1)-1));
    SessN=str2num(File_Name(septag(3)-1));
    DrugC=(File_Name(septag(3)+1:septag(3)+3));
    
    table=readtable([save_path 'CTET_behav_' File_Name(1:end-4) '.txt']);
    NT_ITI(nF)=mode(table.ITI);
    TG_ITI(nF)=mode(table.ITI(table.ITI>mode(table.ITI)+0.1));
    
    fprintf('... processing %s\n', File_Name(1:end-4))
    if length(unique(table.BlockN))~=10
        continue;
    end
    if isempty(all_ITI_info) 
        all_ITI_info(1,:)=nan(1,9);
        all_ITI_info(end,[1 1+match_str(ColorsDlabels,DrugC) 5+match_str(ColorsDlabels,DrugC)])=[SubN NT_ITI(nF) TG_ITI(nF)];
    elseif isempty(all_ITI_info(all_ITI_info(:,1)==SubN,:))
        all_ITI_info=[all_ITI_info ; nan(1,9)];
        all_ITI_info(end,[1 1+match_str(ColorsDlabels,DrugC) 5+match_str(ColorsDlabels,DrugC)])=[SubN NT_ITI(nF) TG_ITI(nF)];
    else
        all_ITI_info(all_ITI_info(:,1)==SubN,[1+match_str(ColorsDlabels,DrugC) 5+match_str(ColorsDlabels,DrugC)])=[NT_ITI(nF) TG_ITI(nF)];
    end
end

%% Get Rid of incomplete sessions
IncompleSubjects=sum(isnan(all_ITI_info(:,2)),2);
fprintf('... ... %g subjects with incomplete sessions\n',sum(IncompleSubjects>0))
%%
IncorrectITI=nanmax(all_ITI_info(:,2:5),[],2);
fprintf('... ... %g subjects with incorrect ITI\n',sum(IncorrectITI>0.9))

fprintf('... ... %g subjects with incomplete sessions OR incorrect ITI\n',sum(IncompleSubjects>0 | IncorrectITI>0.9))

FinalSubjectsID=all_ITI_info(IncompleSubjects==0 & IncorrectITI<0.9,1)';