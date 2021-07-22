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

%% Layout
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

load ../EEG_Cleaning.mat
redo=0;

table_SW=readtable([save_path filesep 'CTET_SWdetection_thr90_byE_P2P_behav_TF_vec_v6.txt']);
table_SW.ERP=nan(size(table_SW,1),1);
%%
winTime=[0.05 0.3];

nFc=0;
all_ERP_NT=[];
all_ERP_TG=[];
all_DrugC=[];
for nF=1:length(ListSubjectsID)
    for nDrug=1:4
        thisF=find_trials({files.name},sprintf('CIcfe_ft_%02d_ctet_session._%s.mat',ListSubjectsID(nF),ColorsDlabels{nDrug}));
        if ~isempty(thisF)
            File_Name=files(thisF).name;
            File_Name=File_Name(10:end);
            septag=findstr(File_Name,'_');
            SubN=str2num(File_Name(1:septag(1)-1));
            SessN=str2num(File_Name(septag(3)-1));
            DrugC=(File_Name(septag(3)+1:septag(3)+3));
            
            if ~ismember(SubN,ListSubjectsID) || strcmp(File_Name(1:end-4),'37_ctet_session2_ATM')
                fprintf('... %s not in subject list\n',File_Name);
                continue;
            end
            
            table=readtable([save_path 'CTET_behav_' File_Name(1:end-4) '.txt']);
            
            if nDrug==1
                nFc=nFc+1;
            end
            fprintf('... processing %s\n',File_Name);
            load([data_path filesep 'CIcfe_ft_' File_Name(1:end-4)]);
            
            %         cfg=[];
            %         cfg.reref           = 'yes';
            %         cfg.refchannel      = 'all';
            %         cfg.demean          = 'yes';
            %         cfg.baselinewindow  = [-0.2 0];
            %         data_clean = ft_preprocessing(cfg,data_clean);
            
            thisF=match_str(EEGCleaning.File,['fe_ft_' File_Name(1:end-4) '.mat']);
            eval(['badTrials=[' EEGCleaning.Trials{thisF} '];']);
            table(ismember(table.TrialN,badTrials),:)=[];
            table.StimType(find(diff(table.BlockN)==1 | diff(table.StimType)==1))=NaN;
            
            for nB=1:10
                cfgerp2        = [];
                cfgerp2.trials = find(table.StimType==0 & table.BlockN==nB)+1; cfgerp2.trials(cfgerp2.trials>length(data_clean.trial))=[];
                av_data_NT_offset = ft_timelockanalysis(cfgerp2, data_clean);
                cfgerp2.trials = find(table.StimType==1 & table.BlockN==nB)+1; cfgerp2.trials(cfgerp2.trials>length(data_clean.trial))=[];
                av_data_TG_offset = ft_timelockanalysis(cfgerp2, data_clean);
                
                ERP_NT_offset=mean(av_data_NT_offset.avg(:,av_data_NT_offset.time>winTime(1) & av_data_NT_offset.time<winTime(2)),2)-mean(av_data_NT_offset.avg(:,av_data_NT_offset.time>-0.2 & av_data_NT_offset.time<0),2);
                ERP_TG_offset=mean(av_data_TG_offset.avg(:,av_data_TG_offset.time>winTime(1) & av_data_TG_offset.time<winTime(2)),2)-mean(av_data_TG_offset.avg(:,av_data_TG_offset.time>-0.2 & av_data_TG_offset.time<0),2);
                
                for nE=1:length(ERP_NT_offset)
                    if sum(table_SW.BlockN==nB & table_SW.SubID==SubN & table_SW.SessN==SessN & ismember(table_SW.Elec,data_clean.label{nE}))~=0
                        table_SW.ERP(table_SW.BlockN==nB & table_SW.SubID==SubN & table_SW.SessN==SessN & ismember(table_SW.Elec,data_clean.label{nE}))=ERP_TG_offset(nE)-ERP_NT_offset(nE);
                    end
                end
            end
        end
    end
end
%%
table_SW.SubID=categorical(table_SW.SubID);
table_SW.SessN=categorical(table_SW.SessN);
table_SW.Elec=categorical(table_SW.Elec);
table_SW.Drug=categorical(table_SW.Drug);
table_SW.Drug=reordercats(table_SW.Drug,[4 1 2 3]);

%%
subIDs=unique(table_SW.SubID);
sessIDs=unique(table_SW.SessN);
Elecs=unique(table_SW.Elec);
table_avSW=table_SW(table_SW.BlockN==1,:);
myCols=[5:9 11:14 16:18];
for i=1:length(subIDs)
    for j=1:length(sessIDs)
        for k=1:length(Elecs)
            table_avSW(table_avSW.SubID==subIDs(i) & table_avSW.SessN==sessIDs(j) & table_avSW.Elec==Elecs(k),[5:9 11:14 16:18])=...
                array2table(nanmean(table2array(table_SW(table_SW.SubID==subIDs(i) & table_SW.SessN==sessIDs(j) & table_SW.Elec==Elecs(k),[5:9 11:14 16:18])),1));
        end
    end
end


%%
for nF=1:length(ListSubjectsID)
    for nDrug=1:4
        thisF=find_trials({files.name},sprintf('CIcfe_ft_%02d_ctet_session._%s.mat',ListSubjectsID(nF),ColorsDlabels{nDrug}));
        if ~isempty(thisF)
            File_Name=files(thisF).name;
            File_Name=File_Name(10:end);
            septag=findstr(File_Name,'_');
            SubN=str2num(File_Name(1:septag(1)-1));
            SessN=str2num(File_Name(septag(3)-1));
            DrugC=(File_Name(septag(3)+1:septag(3)+3));
            
            if ~ismember(SubN,ListSubjectsID) %|| strcmp(File_Name(1:end-4),'37_ctet_session2_ATM')
                fprintf('... %s not in subject list\n',File_Name);
                continue;
            end
            
            table=readtable([save_path 'CTET_behav_' File_Name(1:end-4) '.txt']);
            fprintf('Sub: %2.0f Sess: %g Drug: %s Target: %3.0f Non-Target: %4.0f ITI-TG: %1.3f ITI-NT: %1.3f \n',SubN,SessN,DrugC,sum(table.StimType==1),sum(table.StimType==0),nanmedian(table.ITI(table.StimType==1)),nanmedian(table.ITI(table.StimType==0)));
        end
    end
end
