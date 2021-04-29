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

%% Layout
table_SW=readtable([save_path filesep 'CTET_SWdetection_thr90_byE_P2P_behav_vec_full_v3.txt']);
cfg = [];
cfg.layout = 'biosemi64.lay';
% cfg.channel=unique(table_SW.Elec);
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

load ../EEG_Cleaning.mat
load('../ICA_Artifacts.mat')

%%
redo=1;
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
    
    if redo==1 || exist([data_path filesep 'CIcfe_blocks_ft_' File_Name(1:end-4) '.mat'])==0
        
        %%% interporate channels
        thisF=match_str(EEGCleaning.File,['fe_ft_' File_Name(1:end-4) '.mat']);
        badChannelsStr=EEGCleaning.Channels{thisF};
        badChannels=[];
        if ~isempty(badChannelsStr)
            if sum(badChannelsStr==',')~=0
                comasIdx=find(badChannelsStr==',');
                badChannels=cell(0,0);
                badChannels{1}=badChannelsStr(1:comasIdx(1)-1);
                for k=1:length(comasIdx)-1
                    badChannels{k+1}=badChannelsStr(comasIdx(k)+1:comasIdx(k+1)-1);
                end
                badChannels{length(comasIdx)+1}=badChannelsStr(comasIdx(length(comasIdx))+1:end);
                for k=1:length(badChannels)
                    badChannels{k}(badChannels{k}==' ')=[];
                end
            else
                badChannels{1}=char(badChannelsStr);
                badChannels{1}(badChannels{1}==' ')=[];
            end
        else
            badChannels=[];
        end
        % clean from channels not included in layout
        if ~isempty(badChannels)
            remove=nan(1,length(badChannels));
            for k=1:length(badChannels)
                if ~ismember(badChannels{k},layout.label)
                    remove(k)=1;
                else
                    remove(k)=0;
                end
            end
            badChannels(remove==1)=[];
        end
        
        
        
        %%% Define epochs
        cfg=[];
        cfg.trialfun            = 'CTET_blockfun';
        cfg.table               = table;
        cfg.SubID               = SubN;
        cfg.SessN               = SessN;
        cfg.dataset             = [data_path filesep File_Name];
        cfg.trialdef.prestim    = 2;
        cfg.trialdef.poststim   = 1;
        cfg = ft_definetrial(cfg);
        
        cfg.channel        = hdr.label(match_str(hdr.chantype,'eeg'));
        cfg.demean         = 'yes';
        cfg.lpfilter       = 'yes';        % enable high-pass filtering
        cfg.lpfilttype     = 'but';
        cfg.lpfiltord      = 4;
        cfg.lpfreq         = 40;
        cfg.hpfilter       = 'yes';        % enable high-pass filtering
        cfg.hpfilttype     = 'but';
        cfg.hpfiltord      = 4;
        cfg.hpfreq         = 0.1;
        cfg.dftfilter      = 'yes';        % enable notch filtering to eliminate power line noise
        cfg.dftfreq        = [50]; % set up the frequencies for notch filtering
        dat                   = ft_preprocessing(cfg); % read raw data
        
        cfg.resamplefs      = 256;
        cfg.detrend         = 'no';
        cfg.demean          = 'yes';
        cfg.baselinewindow  = [-0.5 0];
        data = ft_resampledata(cfg, dat);
        
        if ~isempty(badChannels)
            fprintf('... ... interpolating %g channels\n',length(badChannels))
            % find neighbours
            cfg=[];
            cfg.method        = 'triangulation';
            cfg.layout        = layout;
            cfg.feedback      = 'no';
            cfg.channel = layout.label;
            [neighbours] = ft_prepare_neighbours(cfg);
            
            % interpolate channels
            cfg=[];
            cfg.method         = 'weighted';
            cfg.badchannel     = badChannels;
            cfg.missingchannel = [];
            cfg.neighbours     = neighbours;
            cfg.trials         = 'all';
            cfg.layout         = layout;
            cfg.channel = layout.label;
            [data] = ft_channelrepair(cfg, data);
        end
        
        cfg=[];
        cfg.reref           = 'yes';
        cfg.refchannel      = 'all';
        cfg.demean          = 'yes';
        cfg.baselinewindow  = [-2 0];
        data_clean = ft_preprocessing(cfg,data);
        
        %%% Reject bad component
        load([data_path filesep 'Icfe_ft_' files(nF).name(1:end-4) '.mat'], 'comp');
        cfg = [];
        this_line=match_str(ICA_Artifacts.File,['Icfe_ft_' files(nF).name(1:end-4) '.mat']);
        these_components=ICA_Artifacts.OcularArtifact{this_line};
        eval(sprintf('cfg.component = [%s];',these_components)); % to be removed component(s)
        data_clean = ft_rejectcomponent(cfg, comp, data_clean);
        
        
        save([data_path filesep 'CIcfe_blocks_ft_' File_Name(1:end-4)],'data_clean');
    end
    %     clear dat
    %     cfg2=[];
    %     av_data = ft_timelockanalysis(cfg2, data);
end

