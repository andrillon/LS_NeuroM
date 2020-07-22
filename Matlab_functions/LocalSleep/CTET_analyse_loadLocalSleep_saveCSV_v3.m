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
tableblock=readtable('/Users/tand0009/Data/CTET_Dockree/CTET_behav_resblock.txt');

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
    paramSW.byElec=1;
    
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
            if paramSW.byElec
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
all_slow_Waves=[];
all_drug_types=[];
all_slow_Waves_vec=[];
all_drug_types_vec=[];
all_slow_Waves_vec2=[];
all_drug_types_vec2=[];
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
    for nbl=1:10
        temp_table=table(table.SubID==SubN & table.SessN==SessN & ismember(table.BlockN,num2str(nbl)),:);
        temp_table2=tableblock(tableblock.SubID==SubN & tableblock.SessN==SessN & tableblock.BlockN==nbl,:);
        if isempty(temp_table)
            continue;
        end
        min_sample=min(temp_table.Sample);
        max_sample=max(temp_table.Sample);
        temp_slow_Waves=slow_Waves(slow_Waves(:,5)>min_sample & slow_Waves(:,5)<max_sample,:);
        nout=histc(temp_slow_Waves(:,3),1:64);
        temp_P2P=[];
        temp_NegSl=[];
        temp_PosSl=[];
        for nE=1:64
            temp_P2P(nE)=nanmean(temp_slow_Waves(temp_slow_Waves(:,3)==nE,4));
            temp_NegSl(nE)=nanmean(temp_slow_Waves(temp_slow_Waves(:,3)==nE,12));
            temp_PosSl(nE)=nanmean(temp_slow_Waves(temp_slow_Waves(:,3)==nE,13));
        end
        nout(match_str(hdr.label,'Iz'))=NaN;
        temp_P2P(match_str(hdr.label,'Iz'))=NaN;
        temp_NegSl(match_str(hdr.label,'Iz'))=NaN;
        temp_PosSl(match_str(hdr.label,'Iz'))=NaN;
        
        all_slow_Waves=[all_slow_Waves ; [nF SubN SessN nbl nout'/((max_sample-min_sample)/hdr.Fs/60)]];
        all_drug_types=[all_drug_types ; {DrugC}];
        all_slow_Waves_vec=[all_slow_Waves_vec ; [repmat([nF SubN SessN nbl table2array(temp_table2(1,4:8))],64,1) (1:64)' nout/((max_sample-min_sample)/hdr.Fs/60) temp_P2P' temp_NegSl' temp_PosSl']];
        all_drug_types_vec=[all_drug_types_vec ; repmat({DrugC},64,1)];
        all_slow_Waves_vec2=[all_slow_Waves_vec2 ; [repmat([nF SubN SessN nbl table2array(temp_table2(1,4:8))],1,1) 0 nanmean(nout/((max_sample-min_sample)/hdr.Fs/60)) nanmean(temp_P2P) nanmean(temp_NegSl) nanmean(temp_PosSl)]];
        all_drug_types_vec2=[all_drug_types_vec2 ; repmat({DrugC},1,1)];
    end
end

% %%
% figure;
% Drugs={'PLA','ATM','CIT','MPH'};
% for nDrug=1:4
%
%     tempSess=(nanmean(all_slow_Waves(ismember(all_drug_types,Drugs{nDrug}),4:67),2));
%
%     simpleBarPlot(nDrug,tempSess,'k',0.9,'r',[],3);
% end
% set(gca,'XTick',1:4,'XTickLabel',Drugs)
% ylim([5 7])
%%
table_Thr=array2table(sw_thr,'VariableNames',{'FileN','Thr','Elec'});
table_Thr.Elec=categorical(table_Thr.Elec);
for nE=1:64
    table_Thr.Elec(table_Thr.Elec==num2str(nE))=hdr.label{nE};
end
table_Thr(ismember(table_Thr.Elec,{'Iz','TP7','TP8','P9','P10'}),:)=[];
table_Thr.Elec=removecats(table_Thr.Elec);

table_SW=array2table(all_slow_Waves_vec,'VariableNames',{'FileN','SubID','SessN','BlockN','CR','FA','Hit','Miss','Hit_RT','Elec','SWdens','P2P','NegSlope','PosSlope'});
table_SW.SubID=categorical(table_SW.SubID);
table_SW.SessN=categorical(table_SW.SessN);
table_SW.BlockN=ordinal(table_SW.BlockN);
table_SW.Elec=categorical(table_SW.Elec);
for nE=1:64
    table_SW.Elec(table_SW.Elec==num2str(nE))=hdr.label{nE};
end
table_SW.Drug=all_drug_types_vec;
table_SW(ismember(table_SW.Elec,{'Iz','TP7','TP8','P9','P10'}),:)=[];
table_SW.Elec=removecats(table_SW.Elec);
table_SW.Drug=categorical(table_SW.Drug);
table_SW.Drug=reordercats(table_SW.Drug,[4 1 2 3]);

%%% clean for extreme SWdens values
table_SW.SWdens(table_SW.SWdens>(mean(table_SW.SWdens)+5*std(table_SW.SWdens)))=NaN;

mdl0=fitlme(table_SW,'SWdens~1+(1|SubID)');
mdl1=fitlme(table_SW,'SWdens~1+Elec+(1|SubID)');
mdl2=fitlme(table_SW,'SWdens~1+Elec+BlockN+(1|SubID)');
mdl3=fitlme(table_SW,'SWdens~1+Elec*Drug+(1|SubID)');

myS=unique(table_SW.SubID);
myD=unique(table_SW.Drug);
for nS=1:length(myS)
    for nD=1:length(myD)
        nSession(nS,nD)=sum(table_SW.SubID==myS(nS) & table_SW.Drug==myD(nD))~=0;
    end
end
FullSubID=myS(sum(nSession,2)==4);
table_SW2=table_SW(ismember(table_SW.SubID,FullSubID),:);
if paramSW.byElec
    writetable(table_SW,'/Users/tand0009/Data/CTET_Dockree/CTET_SWdetection_thr90_byE_P2P_behav_vec_v3.txt');
    writetable(table_SW2,'/Users/tand0009/Data/CTET_Dockree/CTET_SWdetection_thr90_byE_P2P_behav_vec_full_v3.txt');
else
    writetable(table_SW,'/Users/tand0009/Data/CTET_Dockree/CTET_SWdetection_thr90_allE_P2P_behav_vec_v3.txt');
    writetable(table_SW2,'/Users/tand0009/Data/CTET_Dockree/CTET_SWdetection_thr90_allE_P2P_behav_vec_full_v3.txt');
end
table_allSW=table_SW;
%%
table_SW=array2table(all_slow_Waves_vec2,'VariableNames',{'FileN','SubID','SessN','BlockN','CR','FA','Hit','Miss','Hit_RT','Elec','SWdens','P2P','NegSlope','PosSlope'});
table_SW.SubID=categorical(table_SW.SubID);
table_SW.SessN=categorical(table_SW.SessN);
table_SW.BlockN=ordinal(table_SW.BlockN);
table_SW.Elec=categorical(table_SW.Elec);
for nE=1:64
    table_SW.Elec(table_SW.Elec==num2str(nE))=hdr.label{nE};
end
table_SW.Drug=all_drug_types_vec2;
table_SW(ismember(table_SW.Elec,{'Iz','TP7','TP8','P9','P10'}),:)=[];
table_SW.Elec=removecats(table_SW.Elec);
table_SW.Drug=categorical(table_SW.Drug);
table_SW.Drug=reordercats(table_SW.Drug,[4 1 2 3]);

%%% clean for extreme SWdens values
table_SW.SWdens(table_SW.SWdens>(mean(table_SW.SWdens)+5*std(table_SW.SWdens)))=NaN;

mdl0=fitlme(table_SW,'SWdens~1+(1|SubID)');
mdl1=fitlme(table_SW,'SWdens~1+Elec+(1|SubID)');
mdl2=fitlme(table_SW,'SWdens~1+Elec+BlockN+(1|SubID)');
mdl3=fitlme(table_SW,'SWdens~1+Elec*Drug+(1|SubID)');
table_SW2=table_SW(ismember(table_SW.SubID,FullSubID),:);
if paramSW.byElec
    writetable(table_SW,'/Users/tand0009/Data/CTET_Dockree/CTET_SWdetection_thr90_byE_P2P_avDens_behav_vec_v3.txt');
    writetable(table_SW2,'/Users/tand0009/Data/CTET_Dockree/CTET_SWdetection_thr90_byE_P2P_avDens_behav_vec_full_v3.txt');
    
else
    writetable(table_SW,'/Users/tand0009/Data/CTET_Dockree/CTET_SWdetection_thr90_allE_P2P_avDens_behav_vec_v3.txt');
    writetable(table_SW2,'/Users/tand0009/Data/CTET_Dockree/CTET_SWdetection_thr90_allE_P2P_avDens_behav_vec_full_v3.txt');
end

%%
path_fieldtrip='/Users/tand0009/Work/local/fieldtrip/';
addpath(path_fieldtrip);
ft_defaults;

cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel={'Fp1','AF7','AF3','F1','F3','F5','F7','FT7','FC5','FC3','FC1','C1','C3','C5','T7','CP5','CP3','CP1','P1','P3','P5','P7','PO7','PO3','O1','Oz','POz','Pz','CPz','Fpz','Fp2','AF8','AF4','AFz','Fz','F2','F4','F6','F8','FT8','FC6','FC4','FC2','FCz','Cz','C2','C4','C6','T8','CP6','CP4','CP2','P2','P4','P6','P8','PO8','PO4','O2'}; %setdiff(hdr.label,{'Iz','TP7','TP8'});
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

Drugs={'PLA','CIT','MPH','ATM'};
limMax=16;
figure;
for nDrug=1:4
    subplot(1,4,nDrug);
    temp_topo=[];
    for nE=1:length(layout.label)-2
        temp_topo(nE)=nanmean(table_allSW.SWdens(table_allSW.Drug==Drugs{nDrug} & table_allSW.Elec==layout.label{nE}),1);
    end
    %     temp_topo=(nanmean(all_slow_Waves(:,4:67),1));
    simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
    title(Drugs{nDrug}); h=colorbar;  ylabel(h, 'waves/min')
    caxis([2 6])
    h=colorbar;
    %     set(h,'Position',[0.85 0.7 0.04 0.2])
    format_fig;
    %     set(h,'FontSize',22);
end

figure;
temp_topo=[];
for nE=1:length(layout.label)-2
    temp_topo(nE)=nanmean(table_Thr.Thr(table_Thr.Elec==layout.label{nE}),1);
end
%     temp_topo=(nanmean(all_slow_Waves(:,4:67),1));
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title('Thr'); h=colorbar;  ylabel(h, 'waves/min')
% caxis([2 6])
h=colorbar;
%     set(h,'Position',[0.85 0.7 0.04 0.2])
format_fig;