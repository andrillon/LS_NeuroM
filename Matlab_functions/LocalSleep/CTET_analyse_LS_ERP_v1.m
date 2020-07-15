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
            %             thr_Wave=prctile(all_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
            thr_Wave=prctile(thisE_Waves(:,paramSW.AmpCriterionIdx),paramSW.prticle_Thr);
        end
        sw_thr=[sw_thr ; [SubN thr_Wave nE]];
        
    end
end

%%
nFc=0;
res_mat=[];
drug_cond=[];
EOI={'Fz','Cz','Pz','Oz'};
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
    data=ft_read_data([files(nF).folder filesep files(nF).name]);
    hdr=ft_read_header([files(nF).folder filesep files(nF).name]);
    if hdr.Fs~=1024
        continue;
    end
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
    nFc=nFc+1;
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
    
    ERP_SW=cell(1,length(EOI));
    for nChan=1:length(EOI)
        chanidx=find(ismember(hdr.label,EOI{nChan}));
        chanidx1=find(ismember(hdr.label,'M1'));
        chanidx2=find(ismember(hdr.label,'M2'));
        temp_data=data(chanidx,:)-0.5*data(chanidx1,:)-0.5*data(chanidx2,:);
        temp_data=bandpass(temp_data,hdr.Fs,0.1,30,3);
        
        for nbl=1:10
            temp_table=table(table.SubID==SubN & table.SessN==SessN & ismember(table.BlockN,num2str(nbl)),:);
            temp_table2=tableblock(tableblock.SubID==SubN & tableblock.SessN==SessN & tableblock.BlockN==nbl,:);
            if isempty(temp_table)
                continue;
            end
            min_sample=min(temp_table.Sample);
            max_sample=max(temp_table.Sample);
            
            temp_slow_Waves=slow_Waves(slow_Waves(:,3)==chanidx & slow_Waves(:,5)>min_sample & slow_Waves(:,5)<max_sample,:);
            for nW=1:size(temp_slow_Waves,1)
                if min((-0.3*hdr.Fs:hdr.Fs)+temp_slow_Waves(nW,5))>0 && max((-0.3*hdr.Fs:hdr.Fs)+temp_slow_Waves(nW,5))<length(temp_data)
                    temp=temp_data((-0.25*hdr.Fs:hdr.Fs)+temp_slow_Waves(nW,5));
                    temp=temp-mean(temp(1:0.25*hdr.Fs));
                    ERP_SW{nChan}=[ERP_SW{nChan} ; temp];
                end
            end
        end
        if size(ERP_SW{nChan},1)<30
            mean_ERP_SW(nFc,nChan,:)=nan(1,1281);
        else
            mean_ERP_SW(nFc,nChan,:)=nanmean(ERP_SW{nChan},1);
        end
        mean_ERP_Drugs{nFc}={DrugC};
        
    end
    
end

%%
xTime=-0.25:1/hdr.Fs:1;
figure; hold on;
for nCh=1:4
    %     subplot(2,2,nCh); format_fig;
    plot(xTime,squeeze(nanmean(mean_ERP_SW(:,nCh,:),1)));
    title(EOI{nCh});
    format_fig;
    xlim([-0.2 .8])
end
legend(EOI)