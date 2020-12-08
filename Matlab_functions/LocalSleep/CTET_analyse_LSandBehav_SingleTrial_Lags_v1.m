%%
clear all
close all

path_fieldtrip='/Users/tand0009/Work/local/fieldtrip/';
path_localsleep='/Users/tand0009/WorkGit/projects/inprogress/wanderIM/localsleep';
addpath(path_fieldtrip);
addpath(path_localsleep);
ft_defaults;

path_LSCPtools='/Users/tand0009/WorkGit/LSCPtools/';
addpath(genpath(path_LSCPtools));

data_path='/Volumes/tLab_BackUp1/Monash/CTET_Dockree/EEG_CTET/';
save_path='/Users/tand0009/Data/CTET_Dockree/';
% data_path='/Volumes/shared/R-MNHS-SPP/Bellgrove-data/Jess Barnes EEG Backup Data/EEG_CTET/';
files=dir([data_path filesep '*.bdf']);
filesPLA=dir([data_path filesep '*PLA.bdf']);

table=readtable('/Users/tand0009/Data/CTET_Dockree/CTET_behav_res.txt');

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
all_SWorder=nan(size(table,1),65);

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
    temp_all_SWflag=nan(size(temp_table,1),65);
    temp_all_SWflag2=nan(size(temp_table,1),65);
    if isempty(temp_table)
        continue;
    end
    for nTr=1:size(temp_table,1)
        temp_trial=temp_table(nTr,:);
        beg_sample=(temp_trial.Sample)+0*hdr.Fs;
        end_sample=(temp_trial.Sample)+1.5*hdr.Fs;
        temp_slow_Waves=slow_Waves(:,5)>beg_sample & slow_Waves(:,5)<end_sample;
        tempElecs=(slow_Waves(temp_slow_Waves,3));
        tempIdx=find(temp_slow_Waves);
        uniqueElecs=unique(tempElecs);
        if length(uniqueElecs)<=1
            continue;
        end
        tempLags=nan(1,length(uniqueElecs));
        for nEl=1:length(uniqueElecs)
            tempLags(nEl)=min(slow_Waves(tempIdx(tempElecs==uniqueElecs(nEl)),8));
        end
        temp_all_SWflag(nTr,65)=(nanmin(tempLags)-beg_sample)./hdr.Fs;
        tempLags=(tempLags-nanmin(tempLags))./hdr.Fs;
        uniqueElecs(tempLags>0.5)=[];
        tempLags(tempLags>0.5)=[];
        temp_all_SWflag(nTr,uniqueElecs)=tempLags;
        [~,order]=sort(tempLags);
        temp_all_SWflag2(nTr,uniqueElecs(order))=1:length(order);
    end
    all_SWflag(table.SubID==SubN & table.SessN==SessN,:)=temp_all_SWflag;
    all_SWorder(table.SubID==SubN & table.SessN==SessN,:)=temp_all_SWflag2;
end
%%
load ../lsneurom_biosemi_cleanlayout.mat
myElecs=ismember(hdr.label,layout.label);
all_SWflag2=all_SWflag(:,myElecs);
all_SWorder2=all_SWorder(:,myElecs);

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
% writetable(table_BehavSW,'/Users/tand0009/Data/CTET_Dockree/CTET_SWlags_perTrial_byElec.txt');

%%
figure;
temp_topo=nanmean(all_SWflag2,1)';
simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);

%%
figure;
Drugs={'PLA','ATM','CIT','MPH'};
for nD=1:4
    subplot(1,4,nD);
    temp_topo=nanmean(all_SWflag2(table_BehavSW.Treatment==Drugs{nD},:),1)';
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    title(Drugs{nD});
%     caxis([0.42 0.5])
end

%%
figure;
Drugs={'PLA','ATM','CIT','MPH'};
for nD=1:4
    subplot(1,4,nD);
    temp_topo=nanmean(all_SWorder2(table_BehavSW.Treatment==Drugs{nD},:),1)';
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    title(Drugs{nD});
%     caxis([7 10])
end

%%
figure;
Drugs={'PLA','ATM','CIT','MPH'};
for nD=1:4
    subplot(1,4,nD);
    temp=all_SWorder2(table_BehavSW.Treatment==Drugs{nD},:);
    temp2=double(temp==1);
%     temp2(isnan(temp))=nan;
    temp_topo=nanmean(temp2,1)';
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    title(Drugs{nD});
%     caxis([0.1 0.28])
end