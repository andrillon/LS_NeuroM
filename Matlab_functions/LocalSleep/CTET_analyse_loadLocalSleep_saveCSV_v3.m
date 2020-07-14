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
        all_slow_Waves=[all_slow_Waves ; [nF SubN SessN nbl nout'/((max_sample-min_sample)/hdr.Fs/60)]];
        all_drug_types=[all_drug_types ; {DrugC}];
        all_slow_Waves_vec=[all_slow_Waves_vec ; [repmat([nF SubN SessN nbl table2array(temp_table2(1,4:8))],64,1) (1:64)' nout/((max_sample-min_sample)/hdr.Fs/60) temp_P2P' temp_NegSl' temp_PosSl']];
        all_drug_types_vec=[all_drug_types_vec ; repmat({DrugC},64,1)];
        all_slow_Waves_vec2=[all_slow_Waves_vec2 ; [repmat([nF SubN SessN nbl table2array(temp_table2(1,4:8))],1,1) 0 mean(nout/((max_sample-min_sample)/hdr.Fs/60)) nanmean(temp_P2P) nanmean(temp_NegSl) nanmean(temp_PosSl)]];
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
table_SW=array2table(all_slow_Waves_vec,'VariableNames',{'FileN','SubID','SessN','BlockN','CR','FA','Hit','Miss','Hit_RT','Elec','SWdens','P2P','NegSlope','PosSlope'});
table_SW.SubID=categorical(table_SW.SubID);
table_SW.SessN=categorical(table_SW.SessN);
table_SW.BlockN=ordinal(table_SW.BlockN);
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
writetable(table_SW,'/Users/tand0009/Data/CTET_Dockree/CTET_SWdetection_thr90_allE_P2P_behav_vec_v2.txt');
table_SW2=table_SW(ismember(table_SW.SubID,FullSubID),:);
writetable(table_SW2,'/Users/tand0009/Data/CTET_Dockree/CTET_SWdetection_thr90_allE_P2P_behav_vec_full_v2.txt');

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
table_SW.Elec=removecats(table_SW.Elec);
table_SW.Drug=all_drug_types_vec2;
table_SW.Drug=categorical(table_SW.Drug);
table_SW.Drug=reordercats(table_SW.Drug,[4 1 2 3]);
mdl0=fitlme(table_SW,'SWdens~1+(1|SubID)');
mdl1=fitlme(table_SW,'SWdens~1+Elec+(1|SubID)');
mdl2=fitlme(table_SW,'SWdens~1+Elec+BlockN+(1|SubID)');
mdl3=fitlme(table_SW,'SWdens~1+Elec*Drug+(1|SubID)');

writetable(table_SW,'/Users/tand0009/Data/CTET_Dockree/CTET_SWdetection_thr90_allE_P2P_avDens_behav_vec_v2.txt');

table_SW2=table_SW(ismember(table_SW.SubID,FullSubID),:);
writetable(table_SW2,'/Users/tand0009/Data/CTET_Dockree/CTET_SWdetection_thr90_allE_P2P_avDens_behav_vec_full_v2.txt');


%%
path_fieldtrip='/Users/tand0009/Work/local/fieldtrip/';
addpath(path_fieldtrip);
ft_defaults;

cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel=setdiff(hdr.label,{'Iz'});
layout=ft_prepare_layout(cfg);

Drugs={'PLA','CIT','MPH','ATM'};
limMax=16;
figure;
for nDrug=1:4
    subplot(1,4,nDrug);
    temp_topo=[];
    for nE=1:63
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
% for nDrug2=2:4
%     subplot(2,4,4+nDrug2);
%     temp_topo=(nanmean(all_slow_Waves(ismember(all_drug_types,Drugs{nDrug2}),4:67),1))-...
%         (nanmean(all_slow_Waves(ismember(all_drug_types,Drugs{1}),4:67),1));
%     simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
%     title([Drugs{nDrug2} '-' Drugs{1}]); h=colorbar;  ylabel(h, 'waves/min')
%     caxis([-1 1]*max(abs(temp_topo)))
%     h=colorbar;
%     format_fig;
% end

%%
cmap=cbrewer('seq','YlOrRd',64);

figure;
temp_topo=[];
for nE=1:64
    temp_topo(nE)=mean(table_SW.SWdens(match_str(table_SW.Elec,layout.label{nE})));
end
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title(Drugs{nDrug}); h=colorbar;  ylabel(h, 'waves/min')
caxis([0 1]*limMax)
h=colorbar;
colormap(cmap);
%     set(h,'Position',[0.85 0.7 0.04 0.2])
format_fig;

OrderEl={'AF3','AF4','AF7','AF8','AFz','C1','C2','C3','C4','C5','C6','CP1','CP2','CP3','CP4','CP5','CP6','CPz','Cz','F1','F2','F3','F4','F5','F6','F7','F8','FC1','FC2','FC3','FC4','FC5','FC6',...
    'FCz','Fp1','Fp2','Fpz','FT7','FT8','Fz','Iz','O1','O2','Oz','P1','P10','P2','P3','P4','P5','P6','P7','P8','P9','PO3','PO4','PO7','PO8','POz','Pz','T7','T8','TP7','TP8'};
Pos=[];
for nE=1:64
    Pos(:,nE)=layout.pos(match_str(layout.label,OrderEl{nE}),:)';
end

%%
myElecs=unique(table_SW.Elec);
for nE=1:length(myElecs)
    fprintf('... %g\n',nE)
    temp_table=table_SW(table_SW.Elec==myElecs(nE) & (table_SW.Drug=='PLA' | table_SW.Drug=='CIT'),:);
    temp_table.Drug=removecats(temp_table.Drug);
    mdl0=fitlme(temp_table,'Miss~1+(1|SubID)');
    mdl1=fitlme(temp_table,'Miss~1+SWdens+(1|SubID)');
    mdl2=fitlme(temp_table,'Miss~1+Drug+(1|SubID)');
    mdl3=fitlme(temp_table,'Miss~1+SWdens+Drug+(1|SubID)');
    res=compare(mdl0,mdl1);
    EffectSW_Miss(nE)=res.pValue(2);
    res=compare(mdl0,mdl2);
    EffectDrug_Miss(nE)=res.pValue(2);
    res=compare(mdl1,mdl2);
    Med_Miss(nE)=res.pValue(2);
    
    mdl0=fitlme(temp_table,'FA~1+(1|SubID)');
    mdl1=fitlme(temp_table,'FA~1+SWdens+(1|SubID)');
    mdl2=fitlme(temp_table,'FA~1+Drug+(1|SubID)');
    res=compare(mdl0,mdl1);
    EffectSW_FA(nE)=res.pValue(2);
    res=compare(mdl0,mdl2);
    EffectDrug_FA(nE)=res.pValue(2);
    res=compare(mdl2,mdl1);
    Med_FA(nE)=res.pValue(2);
    
    
    mdl0=fitlme(temp_table,'Hit_RT~1+(1|SubID)');
    mdl1=fitlme(temp_table,'Hit_RT~1+SWdens+(1|SubID)');
    mdl2=fitlme(temp_table,'Hit_RT~1+Drug+(1|SubID)');
    res=compare(mdl0,mdl1);
    EffectSW_Hit_RT(nE)=res.pValue(2);
    res=compare(mdl0,mdl2);
    EffectDrug_Hit_RT(nE)=res.pValue(2);
    res=compare(mdl2,mdl1);
    Med_Hit_RT(nE)=res.pValue(2);
end

%%
cmap=cbrewer('seq','YlOrRd',64);

figure;
subplot(1,3,1)
simpleTopoPlot_ft(EffectSW_FA', layout,'on',[],0,1);
colorbar; caxis([0 1])

subplot(1,3,2)
simpleTopoPlot_ft(EffectSW_Miss', layout,'on',[],0,1);
colorbar; caxis([0 1])

subplot(1,3,3)
simpleTopoPlot_ft(EffectSW_Hit_RT', layout,'on',[],0,1);
colorbar; caxis([0 1])

%%
clear X Y M
nc=0;
Subs=unique(table_SW.SubID);
sub_table_SW=table_SW((table_SW.Drug=='MPH' | table_SW.Drug=='PLA'),:);
sub_table_SW.Drug=removecats(sub_table_SW.Drug);
% for nS=1:length(Subs)
%     if length(unique(sub_table_SW.Drug(sub_table_SW.SubID==Subs(nS))))==2
%         nc=nc+1;
%         X{nc}={sub_table_SW.Miss(sub_table_SW.Elec=='Fz' & sub_table_SW.SubID==Subs(nS))};
%         Y{nc}={sub_table_SW.Drug(sub_table_SW.Elec=='Fz' & sub_table_SW.SubID==Subs(nS))=='MPH'};
%         M{nc}={sub_table_SW.SWdens(sub_table_SW.Elec=='Fz' & sub_table_SW.SubID==Subs(nS))};
        
        for nBl=1:10
            X(:,nBl)=sub_table_SW.Miss(sub_table_SW.Elec=='Fz'  & double(sub_table_SW.BlockN)==nBl);
            Y(:,nBl)=sub_table_SW.Drug(sub_table_SW.Elec=='Fz'  & double(sub_table_SW.BlockN)==nBl)=='MPH';
            M(:,nBl)=sub_table_SW.SWdens(sub_table_SW.Elec=='Fz'  & double(sub_table_SW.BlockN)==nBl);
        end
%     end
% end
[paths, stats2] = mediation(X', Y', M', 'plots',[],'verbose');