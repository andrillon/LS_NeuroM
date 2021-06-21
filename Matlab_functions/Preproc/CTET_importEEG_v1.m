%%
clear all;
close all;
run ../localdef.m

addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_fooof))
addpath(genpath(path_LSCPtools));

files=dir([data_path filesep '*.bdf']);

% table=readtable([save_path 'CTET_behav_res.txt']);

%% Layout
% table_SW=readtable([save_path filesep 'CTET_SWdetection_thr90_byE_P2P_behav_vec_full_v3.txt']);
cfg = [];
cfg.layout = 'biosemi64.lay';
% cfg.channel=unique(table_SW.Elec);
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);


%%
nFc=0;
redo=1;
for nF=40:length(files)
    File_Name=files(nF).name;
    fprintf('... processing %s\n',File_Name);
    septag=findstr(File_Name,'_');
    SubN=str2num(File_Name(1:septag(1)-1));
    SessN=str2num(File_Name(septag(3)-1));
    DrugC=(File_Name(septag(3)+1:septag(3)+3));
    hdr=ft_read_header([data_path filesep File_Name]);
    
    if exist([save_path 'CTET_behav_' File_Name(1:end-4) '.txt'])==0
        continue;
    end
    table=readtable([save_path 'CTET_behav_' File_Name(1:end-4) '.txt']);

    if hdr.Fs~=1024 || length(unique(table.BlockN))~=10
        continue;
    end
    nFc=nFc+1;
    
    if redo==1 || exist([data_path filesep 'fe_ft_' File_Name(1:end-4) '.mat'])==0
        
        
        %%% Define epochs
        cfg=[];
        cfg.trialfun            = 'CTET_trialfun';
        cfg.table               = table;
        cfg.SubID               = SubN;
        cfg.SessN               = SessN;
        cfg.dataset             = [data_path filesep File_Name];
        cfg.trialdef.prestim    = 0.5;
        cfg.trialdef.poststim   = 1.8;
        cfg = ft_definetrial(cfg);
        
        cfg.channel        = hdr.label(match_str(hdr.chantype,'eeg'));
        cfg.demean         = 'yes';
        cfg.lpfilter       = 'yes';        % enable high-pass filtering
        cfg.lpfilttype     = 'but';
        cfg.lpfiltord      = 4;
        cfg.lpfreq         = 40;
        cfg.hpfilter       = 'no';        % enable high-pass filtering
%         cfg.hpfilttype     = 'but';
%         cfg.hpfiltord      = 4;
%         cfg.hpfreq         = 0.1;
        cfg.dftfilter      = 'yes';        % enable notch filtering to eliminate power line noise
        cfg.dftfreq        = [50]; % set up the frequencies for notch filtering
        dat                   = ft_preprocessing(cfg); % read raw data
        
        cfg.resamplefs      = 256;
        cfg.detrend         = 'no';
        cfg.demean          = 'yes';
        cfg.baselinewindow  = [-0.5 0];
        data = ft_resampledata(cfg, dat);
        save([data_path filesep 'fe_ft_' File_Name(1:end-4)],'data');
    end
%     clear dat
%     cfg2=[];
%     av_data = ft_timelockanalysis(cfg2, data);
end

