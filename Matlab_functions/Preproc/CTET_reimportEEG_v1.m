%%
clear all;
close all;
run ../localdef.m

addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_fooof))
addpath(genpath(path_LSCPtools));

files=dir([data_path filesep 'fe_ft_*.mat']);

table=readtable([save_path 'CTET_behav_res.txt']);

%% Layout
table_SW=readtable([save_path filesep 'CTET_SWdetection_thr90_byE_P2P_behav_vec_full_v3.txt']);
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel=unique(table_SW.Elec);
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

load ../EEG_Cleaning.mat
%%
nFc=0;
for nF=1:length(files)
    File_Name=files(nF).name;
    File_Name=File_Name(7:end);
    fprintf('... processing %s\n',File_Name);
    septag=findstr(File_Name,'_');
    SubN=str2num(File_Name(1:septag(1)-1));
    SessN=str2num(File_Name(septag(3)-1));
    DrugC=(File_Name(septag(3)+1:septag(3)+3));
%     hdr=ft_read_header([data_path filesep File_Name]);
    if length(unique(table.BlockN(table.SubID==SubN & table.SessN==SessN)))~=10
        continue;
    end
    nFc=nFc+1;
    
    if exist([data_path filesep 'cfe_ft_' File_Name(1:end-4) '.mat'])==0
        %%% Define epochs
        load([data_path filesep 'fe_ft_' File_Name(1:end-4)])
        
        %%% take out trial
        thisF=match_str(EEGCleaning.File,['fe_ft_' File_Name(1:end-4) '.mat']);
        eval(['badTrials=[' EEGCleaning.Trials{thisF} '];']);
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
        cfg=[];
        cfg.trials          = setdiff(1:length(data.trial),badTrials);
        data = ft_preprocessing(cfg, data);
        

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
        
        save([data_path filesep 'cfe_ft_' File_Name(1:end-4)],'data');
    end
    %     clear dat
    %     cfg2=[];
    %     av_data = ft_timelockanalysis(cfg2, data);
end

