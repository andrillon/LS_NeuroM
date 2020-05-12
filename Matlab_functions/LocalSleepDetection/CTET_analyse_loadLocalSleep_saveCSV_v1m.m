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
            thr_Wave=prctile(all_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
        end
        sw_thr=[sw_thr ; [SubN thr_Wave nE]];
        
    end
end

%%
res_mat=[];
drug_cond=[];
all_slow_Waves=[];
all_drug_types=[];
all_slow_Waves_vec=[];
all_drug_types_vec=[];
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
    nout=histc(slow_Waves(:,3),1:64);
    all_slow_Waves=[all_slow_Waves ; [nF SubN SessN nout'/(hdr.nSamples/hdr.Fs/60)]];
    all_drug_types=[all_drug_types ; {DrugC}];
    all_slow_Waves_vec=[all_slow_Waves_vec ; [repmat([nF SubN SessN],64,1) (1:64)' nout/(hdr.nSamples/hdr.Fs/60)]];
    all_drug_types_vec=[all_drug_types_vec ; repmat({DrugC},64,1)];
end

%%
figure;
Drugs={'PLA','ATM','CIT','MPH'};
for nDrug=1:4
    
    tempSess=(nanmean(all_slow_Waves(ismember(all_drug_types,Drugs{nDrug}),4:67),2));
    
    simpleBarPlot(nDrug,tempSess,'k',0.9,'r',[],3);
end
set(gca,'XTick',1:4,'XTickLabel',Drugs)
ylim([5 7])
%%
table_SW=array2table(all_slow_Waves_vec,'VariableNames',{'FileN','SubID','SessN','Elec','SWdens'});
table_SW.SubID=categorical(table_SW.SubID);
table_SW.SessN=ordinal(table_SW.SessN);
table_SW.Elec=categorical(table_SW.Elec);
for nE=1:64
    table_SW.Elec(table_SW.Elec==num2str(nE))=hdr.label{nE};
end
table_SW.Elec=removecats(table_SW.Elec);
table_SW.Drug=all_drug_types_vec;
table_SW.Drug=categorical(table_SW.Drug);
table_SW.Drug=reordercats(table_SW.Drug,[4 1 2 3]);
mdl0=fitlme(table_SW,'SWdens~1+(1|SubID)');
mdl1=fitlme(table_SW,'SWdens~1+Elec+(1|SubID)');
mdl2=fitlme(table_SW,'SWdens~1+Elec+Drug+(1|SubID)');
mdl3=fitlme(table_SW,'SWdens~1+Elec*Drug+(1|SubID)');

writetable(table_SW,'/Users/tand0009/Data/CTET_Dockree/CTET_SWdetection_thr90_allE_P2P_vec.txt');

%%
path_fieldtrip='/Users/tand0009/Work/local/fieldtrip/';
addpath(path_fieldtrip);
ft_defaults;

cfg = [];
cfg.layout = 'biosemi64.lay';
layout=ft_prepare_layout(cfg);


limMax=16;
figure;
for nDrug=1:4
    subplot(2,4,nDrug);
    temp_topo=(nanmean(all_slow_Waves(ismember(all_drug_types,Drugs{nDrug}),4:67),1));
    %     temp_topo=(nanmean(all_slow_Waves(:,4:67),1));
    simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
    title(Drugs{nDrug}); h=colorbar;  ylabel(h, 'waves/min')
    caxis([0 1]*limMax)
    h=colorbar;
    %     set(h,'Position',[0.85 0.7 0.04 0.2])
    format_fig;
    %     set(h,'FontSize',22);
end
for nDrug2=2:4
    subplot(2,4,4+nDrug2);
    temp_topo=(nanmean(all_slow_Waves(ismember(all_drug_types,Drugs{nDrug2}),4:67),1))-...
        (nanmean(all_slow_Waves(ismember(all_drug_types,Drugs{1}),4:67),1));
    simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
    title([Drugs{nDrug2} '-' Drugs{1}]); h=colorbar;  ylabel(h, 'waves/min')
    caxis([-1 1]*max(abs(temp_topo)))
    h=colorbar;
    format_fig;
end
