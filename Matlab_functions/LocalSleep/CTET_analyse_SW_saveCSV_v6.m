%%
clear all
close all

run ../localdef.m

addpath(path_fieldtrip);
addpath(path_localsleep);
ft_defaults;

addpath(genpath(path_LSCPtools));

files=dir([data_path filesep '*.bdf']);
filesPLA=dir([data_path filesep '*PLA.bdf']);

table=readtable([save_path filesep 'CTET_behav_res.txt']);
tableblock=readtable([save_path filesep 'CTET_behav_resblock.txt']);

%% Get the thresholds
sw_thr=[];
for nF=1:length(filesPLA)
    File_Name=filesPLA(nF).name;
    fprintf('... processing %s\n',File_Name);
    septag=findstr(File_Name,'_');
    SubN=str2num(File_Name(1:septag(1)-1));
    SessN=str2num(File_Name(septag(3)-1));
    DrugC=(File_Name(septag(3)+1:septag(3)+3));
    
    if exist([data_path filesep '..' filesep 'SW_detection' filesep 'CIcfe_blocks_ft_allSW_' File_Name(1:end-4) '.mat'])==0
        continue;
    end
    load([data_path filesep '..' filesep 'SW_detection' filesep 'CIcfe_blocks_ft_allSW_' File_Name(1:end-4)]); %,'all_Waves','hdr')
    Fs=256;
    
    %%% clean detection
    paramSW.prticle_Thr=90; % 80 or 90 or 95
    paramSW.LimFrqW=[1 4]; % [1 4] or [4 10]
    paramSW.AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=[];
    paramSW.art_ampl=150;
    paramSW.max_posampl=75;
    paramSW.max_Freq=4;
    paramSW.byElec=1;
    
    all_Waves=double(all_Waves);
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./Fs);
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
    
     if sum(sw_thr(:,1)==SubN)==0
        continue;
    end
    if exist([data_path filesep '..' filesep 'SW_detection' filesep 'CIcfe_blocks_ft_allSW_' File_Name(1:end-4) '.mat'])==0
        continue;
    end
    load([data_path filesep '..' filesep 'SW_detection' filesep 'CIcfe_blocks_ft_allSW_' File_Name(1:end-4)]); %,'all_Waves','hdr')
    load([data_path filesep 'CIcfe_blocks_ft_' File_Name(1:end-4)]); %,'all_Waves','hdr')
    cfg=[];
    cfg.reref      = 'yes';
    cfg.refchannel = {'M1','M2'};
    data_clean = ft_preprocessing(cfg,data_clean);
    
    
    all_Waves=double(all_Waves);
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./Fs);
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
        temp_slow_Waves=slow_Waves(slow_Waves(:,2)==nbl,:);
        nout=histc(temp_slow_Waves(:,3),1:64);
        temp_P2P=[];
        temp_NegSl=[];
        temp_PosSl=[];
        for nE=1:64
            temp_P2P(nE)=nanmean(temp_slow_Waves(temp_slow_Waves(:,3)==nE,4));
            temp_NegSl(nE)=nanmean(temp_slow_Waves(temp_slow_Waves(:,3)==nE,12));
            temp_PosSl(nE)=nanmean(temp_slow_Waves(temp_slow_Waves(:,3)==nE,13));
        end
        
        nout(match_str(data_clean.label,'Iz'))=NaN;
        temp_P2P(match_str(data_clean.label,'Iz'))=NaN;
        temp_NegSl(match_str(data_clean.label,'Iz'))=NaN;
        temp_PosSl(match_str(data_clean.label,'Iz'))=NaN;
        
        all_slow_Waves=[all_slow_Waves ; [nF SubN SessN nbl nout'/(size(data_clean.trial{nbl},2)/Fs/60)]];
        all_drug_types=[all_drug_types ; {DrugC}];
        
        temp_table2=table(table.SubID==SubN & table.SessN==SessN & table.BlockN==nbl,:);
        all_slow_Waves_vec=[all_slow_Waves_vec ; [repmat([nF SubN SessN nbl nanmean(table2array(temp_table2(:,7:11)),1)],64,1) (1:64)' nout/(size(data_clean.trial{nbl},2)/Fs/60) temp_P2P' temp_NegSl' temp_PosSl']];
        all_drug_types_vec=[all_drug_types_vec ; repmat({DrugC},64,1)];
        all_slow_Waves_vec2=[all_slow_Waves_vec2 ; [repmat([nF SubN SessN nbl nanmean(table2array(temp_table2(:,7:11)),1)],1,1) 0 nanmean(nout/(size(data_clean.trial{nbl},2)/Fs/60)) nanmean(temp_P2P) nanmean(temp_NegSl) nanmean(temp_PosSl)]];
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
all_slow_Waves_vec(:,6)=1-all_slow_Waves_vec(:,6);
all_slow_Waves_vec2(:,6)=1-all_slow_Waves_vec2(:,6);
%%
table_Thr=array2table(sw_thr,'VariableNames',{'FileN','Thr','Elec'});
table_Thr.Elec=categorical(table_Thr.Elec);
for nE=1:64
    table_Thr.Elec(table_Thr.Elec==num2str(nE))=data_clean.label{nE};
end
table_Thr(ismember(table_Thr.Elec,{'Iz'}),:)=[];
table_Thr.Elec=removecats(table_Thr.Elec);

table_SW=array2table(all_slow_Waves_vec,'VariableNames',{'FileN','SubID','SessN','BlockN','rawRT','FA','Hit','ITI','RT','Elec','SWdens','P2P','NegSlope','PosSlope'});
table_SW.SubID=categorical(table_SW.SubID);
table_SW.SessN=categorical(table_SW.SessN);
table_SW.Elec=categorical(table_SW.Elec);
for nE=1:64
    table_SW.Elec(table_SW.Elec==num2str(nE))=data_clean.label{nE};
end
table_SW.Drug=all_drug_types_vec;
table_SW(ismember(table_SW.Elec,{'Iz'}),:)=[];
table_SW.Elec=removecats(table_SW.Elec);
table_SW.Drug=categorical(table_SW.Drug);
table_SW.Drug=reordercats(table_SW.Drug,[4 1 2 3]);

%%% clean for extreme SWdens values
table_SW.SWdens(table_SW.SWdens>(mean(table_SW.SWdens)+5*std(table_SW.SWdens)))=NaN;

mdl0=fitlme(table_SW,'SWdens~1+(1|SubID)');
mdl1=fitlme(table_SW,'SWdens~1+Elec+(1|SubID)');
mdl2=fitlme(table_SW,'SWdens~1+Elec+BlockN+(1|SubID)');
mdl3=fitlme(table_SW,'SWdens~1+Elec+BlockN+Drug+(1|SubID)');
mdl4=fitlme(table_SW,'SWdens~1+Elec*Drug+BlockN*Drug+(1|SubID)');

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
    writetable(table_SW,[save_path filesep 'CTET_SWdetection_thr90_byE_P2P_behav_vec_v6.txt']);
    writetable(table_SW2,[save_path filesep 'CTET_SWdetection_thr90_byE_P2P_behav_vec_full_v6.txt']);
else
    writetable(table_SW,[save_path filesep 'CTET_SWdetection_thr90_allE_P2P_behav_vec_v6.txt']);
    writetable(table_SW2,[save_path filesep 'CTET_SWdetection_thr90_allE_P2P_behav_vec_full_v6.txt']);
end
table_allSW=table_SW;
%%
table_SW=array2table(all_slow_Waves_vec2,'VariableNames',{'FileN','SubID','SessN','BlockN','rawRT','corrNT','corrTG','ITI','RT','Elec','SWdens','P2P','NegSlope','PosSlope'});
table_SW.SubID=categorical(table_SW.SubID);
table_SW.SessN=categorical(table_SW.SessN);
table_SW.BlockN=ordinal(table_SW.BlockN);
table_SW.Elec=categorical(table_SW.Elec);
for nE=1:64
    table_SW.Elec(table_SW.Elec==num2str(nE))=data_clean.label{nE};
end
table_SW.Drug=all_drug_types_vec2;
table_SW(ismember(table_SW.Elec,{'Iz'}),:)=[];
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
    writetable(table_SW,[save_path filesep 'CTET_SWdetection_thr90_byE_P2P_avDens_behav_vec_v6.txt']);
    writetable(table_SW2,[save_path filesep 'CTET_SWdetection_thr90_byE_P2P_avDens_behav_vec_full_v6.txt']);
    
else
    writetable(table_SW,[save_path filesep 'CTET_SWdetection_thr90_allE_P2P_avDens_behav_vec_v6.txt']);
    writetable(table_SW2,[save_path filesep 'CTET_SWdetection_thr90_allE_P2P_avDens_behav_vec_full_v6.txt']);
end

%%
path_fieldtrip='/Users/tand0009/Work/local/fieldtrip/';
addpath(path_fieldtrip);
ft_defaults;

cfg = [];
cfg.layout = 'biosemi64.lay';
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
    caxis([2 15])
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