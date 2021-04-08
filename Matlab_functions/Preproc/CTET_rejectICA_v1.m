%%
clear all;
close all;
run ../localdef.m

addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_fooof))
addpath(genpath(path_LSCPtools));

files=dir([data_path filesep 'Icfe_ft_*.mat']);

table=readtable([save_path 'CTET_behav_res.txt']);
load('../ICA_Artifacts.mat')

%%
nFc=0;
for nF=1:length(files)
    File_Name=files(nF).name;
    File_Name=File_Name(9:end);
    septag=findstr(File_Name,'_');
    SubN=str2num(File_Name(1:septag(1)-1));
    SessN=str2num(File_Name(septag(3)-1));
    DrugC=(File_Name(septag(3)+1:septag(3)+3));
    if length(unique(table.BlockN(table.SubID==SubN & table.SessN==SessN)))~=10
        continue;
    end
    nFc=nFc+1;
    if exist([data_path filesep 'CIf' File_Name])==0
         fprintf('... processing %s\n',File_Name);
   load([data_path filesep 'Icfe_ft_' File_Name(1:end-4)]);
        
        %%% Reject bad component
        cfg = [];
        this_line=match_str(ICA_Artifacts.File,files(nF).name);
        these_components=ICA_Artifacts.OcularArtifact{this_line};
        eval(sprintf('cfg.component = [%s];',these_components)); % to be removed component(s)
        data_clean = ft_rejectcomponent(cfg, comp, data);
        
        cfg=[];
        cfg.reref           = 'yes';
        cfg.refchannel      = 'all';
        cfg.demean          = 'yes';
        cfg.baselinewindow  = [-0.1 0];
        data_clean = ft_preprocessing(cfg,data_clean);
        data = ft_preprocessing(cfg,data);
        
        save([data_path filesep 'CIf' File_Name],'data_clean');
    end
    %     cfgerp = [];
    %     av_data_preICA = ft_timelockanalysis(cfgerp, data);
    %     av_data_postICA = ft_timelockanalysis(cfgerp, data_clean);
    %     figure; plot(av_data_preICA.time,av_data_preICA.avg(match_str(av_data_preICA.label,'Oz'),:)','k')
    %     hold on
    %     plot(av_data_preICA.time,av_data_postICA.avg(match_str(av_data_preICA.label,'Oz'),:)','r')
end

