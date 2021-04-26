%%
clear all
% close all

run ../localdef.m

addpath(path_fieldtrip);
addpath(path_localsleep);
ft_defaults;

addpath(genpath(path_LSCPtools));

% data_path='/Volumes/shared/R-MNHS-SPP/Bellgrove-data/Jess Barnes EEG Backup Data/EEG_CTET/';
% data_path='/Volumes/shared/R-MNHS-SPP/Bellgrove-data/Jess Barnes EEG Backup Data/EEG_CTET/';
files=dir([data_path filesep '*.bdf']);
filesPLA=dir([data_path filesep '*PLA.bdf']);

% table=readtable('/Users/tand0009/Data/CTET_Dockree/CTET_behav_res.txt');
% tableblock=readtable('/Users/tand0009/Data/CTET_Dockree/CTET_behav_resblock.txt');

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
    paramSW.byElec=0;
    
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

% %%

% 
% figure;
% temp_topo=[];
% for nE=1:length(layout.label)-2
%     findElec=match_str(data_clean.label,layout.label{nE});
%     temp_topo(nE)=nanmean(sw_thr(sw_thr(:,3)==findElec,2));
% end
% simpleTopoPlot_ft(temp_topo', layout,'labels',[],0,1);
% colorbar;
%%
myERP_Elec={'Fz','FCz','Cz','CPz','Oz'};
%%
res_mat=[];
drug_cond=[];
all_slow_Waves=[];
all_drug_types=[];
mean_SW_ERP_byElec=[];
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
        
    %%% clean detection
    paramSW.prticle_Thr=90; % 80 or 90 or 95
    paramSW.LimFrqW=[1 4]; % [1 4] or [4 10]
    paramSW.AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=[];
    paramSW.art_ampl=150;
    paramSW.max_posampl=75;
    paramSW.max_Freq=4;
    
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
    temp_ERP=cell(1,length(myERP_Elec));
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
            
            if ismember(data_clean.label(nE),myERP_Elec)
                % get ERP for slow waves
                temp_SW=temp_slow_Waves(temp_slow_Waves(:,3)==nE,5);
                for m=1:length(temp_SW)
                    if temp_SW(m)-0.5*data_clean.fsample>0 && temp_SW(m)+1*data_clean.fsample<size(data_clean.trial{nbl},2)
                    vec=data_clean.trial{nbl}(nE,(temp_SW(m)-0.5*data_clean.fsample):(temp_SW(m)+1*data_clean.fsample));
                    vec=vec-mean(vec(1:0.5*data_clean.fsample));
                    temp_ERP{find(ismember(myERP_Elec,data_clean.label(nE)))}=[temp_ERP{find(ismember(myERP_Elec,data_clean.label(nE)))} ; vec];
                    end
                end
            end
        end
        nout(match_str(data_clean.label,'Iz'))=NaN;
        temp_P2P(match_str(data_clean.label,'Iz'))=NaN;
        temp_NegSl(match_str(data_clean.label,'Iz'))=NaN;
        temp_PosSl(match_str(data_clean.label,'Iz'))=NaN;
        
        all_slow_Waves=[all_slow_Waves ; [nF SubN SessN nbl nout'/(size(data_clean.trial{nbl},2)/Fs/60)]];
        all_drug_types=[all_drug_types ; {DrugC}];
%         all_slow_Waves_vec=[all_slow_Waves_vec ; [repmat([nF SubN SessN nbl table2array(temp_table2(1,4:9))],64,1) (1:64)' nout/(size(data_clean.trial{nbl},2)/Fs/60) temp_P2P' temp_NegSl' temp_PosSl']];
%         all_drug_types_vec=[all_drug_types_vec ; repmat({DrugC},64,1)];
%         all_slow_Waves_vec2=[all_slow_Waves_vec2 ; [repmat([nF SubN SessN nbl table2array(temp_table2(1,4:9))],1,1) 0 nanmean(nout/(size(data_clean.trial{nbl},2)/Fs/60)) nanmean(temp_P2P) nanmean(temp_NegSl) nanmean(temp_PosSl)]];
%         all_drug_types_vec2=[all_drug_types_vec2 ; repmat({DrugC},1,1)];
    end
    temp_SW_ERP_byElec=[];
    for j=1:length(myERP_Elec)
        temp_SW_ERP_byElec(j,:)=mean(temp_ERP{j});
    end
    mean_SW_ERP_byElec=cat(3,mean_SW_ERP_byElec,temp_SW_ERP_byElec);
end


%%
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel=data_clean.label;
cfg.channel(match_str(cfg.channel,'Iz'))=[];
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);


Drugs={'PLA','CIT','MPH','ATM'};
figure;
for nDrug=1:4
    subplot(2,4,nDrug);
    temp_topo=[];
    for nE=1:length(layout.label)-2
        findElec=match_str(data_clean.label,layout.label{nE});
        temp_topo(nE)=nanmean(all_slow_Waves(match_str(all_drug_types,Drugs{nDrug}),4+findElec),1);
    end
    %     temp_topo=(nanmean(all_slow_Waves(:,4:67),1));
    simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
    title(Drugs{nDrug}); h=colorbar;  ylabel(h, 'waves/min')
    caxis([4 7])
    h=colorbar;
    %     set(h,'Position',[0.85 0.7 0.04 0.2])
    format_fig;
    %     set(h,'FontSize',22);
end
for nDrug=2:4
    subplot(2,4,4+nDrug);
    temp_topo=[];
    for nE=1:length(layout.label)-2
        findElec=match_str(data_clean.label,layout.label{nE});
        temp_topo(nE)=nanmean(all_slow_Waves(match_str(all_drug_types,Drugs{nDrug}),4+findElec),1)-...
            nanmean(all_slow_Waves(match_str(all_drug_types,Drugs{1}),4+findElec),1);
    end
    %     temp_topo=(nanmean(all_slow_Waves(:,4:67),1));
    simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
    title(Drugs{nDrug}); h=colorbar;  ylabel(h, 'waves/min')
    caxis([-1 1]*1.25)
    h=colorbar;
    %     set(h,'Position',[0.85 0.7 0.04 0.2])
    format_fig;
    %     set(h,'FontSize',22);
end

%%
Drugs={'PLA','CIT','MPH','ATM'};
figure;
for nBl=1:10
    subplot(2,5,nBl);
    temp_topo=[];
    for nE=1:length(layout.label)-2
        findElec=match_str(data_clean.label,layout.label{nE});
        temp_topo(nE)=nanmean(all_slow_Waves(all_slow_Waves(:,4)==nBl,4+findElec),1);
    end
    %     temp_topo=(nanmean(all_slow_Waves(:,4:67),1));
    simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
    title(Drugs{nDrug}); h=colorbar;  ylabel(h, 'waves/min')
    caxis([4 6])
    h=colorbar;
    %     set(h,'Position',[0.85 0.7 0.04 0.2])
    format_fig;
    %     set(h,'FontSize',22);
end
%%
figure; 
plot(squeeze(mean(mean_SW_ERP_byElec,3))')
