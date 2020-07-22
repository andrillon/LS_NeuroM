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

%%
res_mat=[];
drug_cond=[];
for nF=1:length(files)
    File_Name=files(nF).name;
    fprintf('... processing %s\n',File_Name);
    septag=findstr(File_Name,'_');
    SubN=str2num(File_Name(1:septag(1)-1));
    SessN=str2num(File_Name(septag(3)-1));
    DrugC=(File_Name(septag(3)+1:septag(3)+3));
    
    % data=ft_read_data('/Volumes/tLab_BackUp1/Monash/CTET_Dockree/EEG_CTET/01_ctet_session1_ATM.bdf');
    hdr=ft_read_header([data_path filesep File_Name]);
    events=ft_read_event([data_path filesep File_Name]);
    dat=ft_read_data([data_path filesep File_Name]);
    
    all_Waves=[];
    
    temp_data=dat(1:64,:);
    temp_data=temp_data-repmat(mean(temp_data(match_str(hdr.label,{'TP7','TP8'}),:),1),size(temp_data,1),1);
    temp_data=temp_data-repmat(mean(temp_data,2),1,size(temp_data,2));
    
    [twa_results]=twalldetectnew_TA_v2(temp_data,hdr.Fs,0);
    for nE=1:64
        all_Waves=[all_Waves ; [repmat([1 nF nE],length(abs(cell2mat(twa_results.channels(nE).maxnegpkamp))),1) abs(cell2mat(twa_results.channels(nE).maxnegpkamp))'+abs(cell2mat(twa_results.channels(nE).maxpospkamp))' ...
            cell2mat(twa_results.channels(nE).negzx)' ...
            cell2mat(twa_results.channels(nE).poszx)' ...
            cell2mat(twa_results.channels(nE).wvend)' ...
            cell2mat(twa_results.channels(nE).maxnegpk)' ...
            cell2mat(twa_results.channels(nE).maxnegpkamp)' ...
            cell2mat(twa_results.channels(nE).maxpospk)' ...
            cell2mat(twa_results.channels(nE).maxpospkamp)' ...
            cell2mat(twa_results.channels(nE).mxdnslp)' ...
            cell2mat(twa_results.channels(nE).mxupslp)' ...
            cell2mat(twa_results.channels(nE).maxampwn)' ...
            cell2mat(twa_results.channels(nE).minampwn)' ...
            ]];
    end
    fprintf('\n')
    save([save_path filesep 'SW_detection' filesep 'PH_CTET_allSW_' File_Name(1:end-4)],'all_Waves','hdr')
    
    %%% clean detection
    paramSW.prticle_Thr=90; % 80 or 90 or 95
    paramSW.LimFrqW=[1 4]; % [1 4] or [4 10]
    paramSW.AmpCriterionIdx=4; % 9 (MaxNegpkAmp) or 11 (MaxPosPeakAmp) or 4 (P2P)
    paramSW.fixThr=75/2;
    paramSW.art_ampl=150;
    paramSW.max_posampl=75;
    paramSW.max_Freq=7;
    
    all_Waves=double(all_Waves);
    all_freq=1./(abs((all_Waves(:,5)-all_Waves(:,7)))./hdr.Fs);
    fprintf('... ... %g %% waves discarded because of frequency\n',mean(all_freq>paramSW.max_Freq)*100)
    fprintf('... ... %g %% waves discarded because of max P2P ampl\n',mean(all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl)*100)
    fprintf('... ... %g %% waves discarded because of max pos ampl\n',mean(all_Waves(:,11)>paramSW.max_posampl | all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl)*100)
    all_Waves(all_freq>paramSW.max_Freq | all_Waves(:,paramSW.AmpCriterionIdx)>paramSW.art_ampl | all_Waves(:,11)>paramSW.max_posampl| all_Waves(:,14)>paramSW.art_ampl| abs(all_Waves(:,15))>paramSW.art_ampl,:)=[];
    
    thr_Wave=[];
    slow_Waves=[];
    for nE=1:64
        thisE_Waves=all_Waves(all_Waves(:,3)==nE,:);
        temp_p2p=thisE_Waves(:,paramSW.AmpCriterionIdx);
        
        if ~isempty(paramSW.fixThr)
            thr_Wave(nE)=paramSW.fixThr;
        else
            thr_Wave(nE)=prctile(thisE_Waves(:,AmpCriterionIdx),paramSW.prticle_Thr);
        end
        slow_Waves=[slow_Waves ; thisE_Waves(temp_p2p>thr_Wave(nE),:)];
    end
    save([save_path filesep 'SW_detection' filesep 'PH_CTET_SW_' File_Name(1:end-4)],'slow_Waves','hdr','paramSW')
    
end
