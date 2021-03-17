%%
clear all;
close all;
run ../localdef.m

addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_fooof))
addpath(genpath(path_LSCPtools));

files=dir([data_path filesep '*.bdf']);

table=readtable([save_path 'CTET_behav_res.txt']);

%%
nFc=0;
for nF=1:length(files)
    File_Name=files(nF).name;
    fprintf('... processing %s\n',File_Name);
    septag=findstr(File_Name,'_');
    SubN=str2num(File_Name(1:septag(1)-1));
    SessN=str2num(File_Name(septag(3)-1));
    DrugC=(File_Name(septag(3)+1:septag(3)+3));
    hdr=ft_read_header([data_path filesep File_Name]);
    if hdr.Fs~=1024 || length(unique(table.BlockN(table.SubID==SubN & table.SessN==SessN)))~=10
        continue;
    end
    nFc=nFc+1;
    
    if exist([data_path filesep 'Icfe_ft_' File_Name(1:end-4) '.mat'])==0
    load([data_path filesep 'cfe_ft_' File_Name(1:end-4)]);
    
    cfg=[];
    cfg.reref      = 'yes';
    cfg.refchannel = 'all';
    data = ft_preprocessing(cfg,data);
    
    %%% run ICA
    rankICA = rank(data.trial{1,1});
    cfg        = [];
    cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB
    cfg.numcomponent = rankICA;
    comp = ft_componentanalysis(cfg, data);
    save([data_path filesep 'Icfe_ft_' File_Name(1:end-4)],'data','comp','rankICA');
    end
%     clear dat
%     cfg2=[];
%     av_data = ft_timelockanalysis(cfg2, data);
end

