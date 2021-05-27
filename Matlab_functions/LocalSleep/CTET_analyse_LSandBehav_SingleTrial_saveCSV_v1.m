%%
clear all
close all

run ../localdef.m
addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_LSCPtools));
addpath(genpath([pwd filesep '..']));

% data_path='/Volumes/shared/R-MNHS-SPP/Bellgrove-data/Jess Barnes EEG Backup Data/EEG_CTET/';
files=dir([data_path filesep '*.bdf']);
filesPLA=dir([data_path filesep '*PLA.bdf']);

table=readtable([save_path 'CTET_behav_res.txt']);

%% Get the thresholds
sw_thr=[];
for nF=1:length(filesPLA)
    File_Name=filesPLA(nF).name;
    fprintf('... processing %s\n',File_Name);
    septag=findstr(File_Name,'_');
    SubN=str2num(File_Name(1:septag(1)-1));
    SessN=str2num(File_Name(septag(3)-1));
    DrugC=(File_Name(septag(3)+1:septag(3)+3));
    
    if exist([save_path filesep 'SW_detection' filesep 'PH_CTET_allSW_' File_Name(1:end-4) '.mat'])==0
        continue;
    end
    load([save_path filesep 'SW_detection' filesep 'PH_CTET_allSW_' File_Name(1:end-4)]); %,'all_Waves','hdr')
    
    %%% clean detection
    paramSW.prticle_Thr=90; % 80 or 90 or 95
    paramSW.LimFrqW=[1 4]; % [1 4] or [4 10]
    paramSW.AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=[];
    paramSW.art_ampl=150;
    paramSW.max_posampl=75;
    paramSW.max_Freq=7;
    paramSW.thr_byEl=1;
    
    all_Waves=double(all_Waves);
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./hdr.Fs);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>paramSW.max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
    all_Waves(all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
    
    slow_Waves=[];
    for nE=1:64
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
        
        if ~isempty(paramSW.fixThr)
            thr_Wave=paramSW.fixThr;
        else
            if paramSW.thr_byEl
            thr_Wave=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
            else
            thr_Wave=prctile(all_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
            end
        end
        sw_thr=[sw_thr ; [SubN thr_Wave nE]];
        
    end
end

%%
res_mat=[];
drug_cond=[];
all_SWflag=nan(size(table,1),65);

for nF=1:length(files)
    File_Name=files(nF).name;
    fprintf('... processing %s\n',File_Name);
    septag=findstr(File_Name,'_');
    SubN=str2num(File_Name(1:septag(1)-1));
    SessN=str2num(File_Name(septag(3)-1));
    DrugC=(File_Name(septag(3)+1:septag(3)+3));
    
    if exist([save_path filesep 'SW_detection' filesep 'PH_CTET_allSW_' File_Name(1:end-4) '.mat'])==0
        continue;
    end
    load([save_path filesep 'SW_detection' filesep 'PH_CTET_allSW_' File_Name(1:end-4)]); %,'all_Waves','hdr')
    
    %%% clean detection
    paramSW.prticle_Thr=90; % 80 or 90 or 95
    paramSW.LimFrqW=[1 4]; % [1 4] or [4 10]
    paramSW.AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=[];
    paramSW.art_ampl=150;
    paramSW.max_posampl=75;
    paramSW.max_Freq=7;
    
    all_Waves=double(all_Waves);
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./hdr.Fs);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>paramSW.max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
    all_Waves(all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
    
    slow_Waves=[];
    for nE=1:64
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
        
        thr_Wave=sw_thr(sw_thr(:,1)==SubN & sw_thr(:,3)==nE,2);
        if isempty(thr_Wave) || length(thr_Wave)>1
            continue;
        end
        slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave,:)];
    end
    %     save([save_path filesep 'SW_detection' filesep 'PH_CTET_SW_' File_Name(1:end-4)],'slow_Waves','hdr','paramSW')
    temp_table=table(table.SubID==SubN & table.SessN==SessN,:);
    temp_all_SWflag=zeros(size(temp_table,1),65);
    if isempty(temp_table)
        continue;
    end
    for nTr=1:size(temp_table,1)
        temp_trial=temp_table(nTr,:);
        beg_sample=(temp_trial.Sample)+0*hdr.Fs;
        end_sample=(temp_trial.Sample)+1.5*hdr.Fs;
        temp_slow_Waves=slow_Waves(:,5)>beg_sample & slow_Waves(:,7)<end_sample;
        temp_all_SWflag(nTr,slow_Waves(temp_slow_Waves,3))=1;
        temp_all_SWflag(nTr,65)=sum(temp_slow_Waves);
    end
    all_SWflag(table.SubID==SubN & table.SessN==SessN,:)=temp_all_SWflag;
end
%%
load ../lsneurom_biosemi_cleanlayout.mat
myElecs=ismember(hdr.label,layout.label);
all_SWflag2=all_SWflag(:,myElecs);

VarNames=[];
myLabels=hdr.label(ismember(hdr.label,layout.label));
for nE=1:length(myLabels)
VarNames{nE}=myLabels{nE};
end
% VarNames{nE+1}='ALL';
table_SW=array2table(all_SWflag2,'VariableNames',VarNames);

table_BehavSW=[table table_SW];
table_BehavSW.Treatment=categorical(table_BehavSW.Treatment);
table_BehavSW.Treatment=reordercats(table_BehavSW.Treatment,[4 1 2 3]);
table_BehavSW(find(ismember(table_BehavSW.BlockN,'<undefined>')),:)=[];
table_BehavSW.BlockN=categorical(table_BehavSW.BlockN);
writetable(table_BehavSW,'/Users/tand0009/Data/CTET_Dockree/CTET_SWflag_perTrial_byElec_v2.txt');

