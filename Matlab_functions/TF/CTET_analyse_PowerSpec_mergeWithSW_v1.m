%%
clear all
close all
run ../localdef.m

path_fieldtrip='/Users/tand0009/Work/local/fieldtrip/';
addpath(path_fieldtrip);
ft_defaults;
fooof_path='/Users/tand0009/WorkGit/projects/ext/fooof_mat/';
addpath(genpath(fooof_path))
path_LSCPtools='/Users/tand0009/WorkGit/LSCPtools/';
addpath(genpath(path_LSCPtools));

% data_path='/Volumes/shared/R-MNHS-SPP/Bellgrove-data/Jess Barnes EEG Backup Data/EEG_CTET/';
files=dir([data_path filesep 'CIcfe_ft_*.mat']);

table_SW=readtable([save_path filesep 'CTET_SWdetection_thr90_byE_P2P_behav_vec_v6.txt']);
table_avSW=readtable([save_path filesep 'CTET_SWdetection_thr90_byE_P2P_avDens_behav_vec_v6.txt']);

%%
res_mat=[];
drug_cond=[];
myfreqs=1:0.1:40;
nFc=0;
redo=0;

all_pow=[];
all_SNR=[];
all_SNRtag=[];
table_avSW.SSVEP=nan(size(table_avSW,1),1);
table_SW.SSVEP=nan(size(table_SW,1),1);
table_avSW.Alpha=nan(size(table_avSW,1),1);
table_SW.Alpha=nan(size(table_SW,1),1);
for nF=1:length(ListSubjectsID)
    for nDrug=1:4
        thisF=find_trials({files.name},sprintf('CIcfe_ft_%02d_ctet_session._%s.mat',ListSubjectsID(nF),ColorsDlabels{nDrug}));
        if ~isempty(thisF)
            
            File_Name=files(thisF).name;
            File_Name=File_Name(10:end);
            septag=findstr(File_Name,'_');
            SubN=str2num(File_Name(1:septag(1)-1));
            SessN=str2num(File_Name(septag(3)-1));
            DrugC=(File_Name(septag(3)+1:septag(3)+3));
            
            if ~ismember(SubN,ListSubjectsID)
                fprintf('... %s not in subject list\n',File_Name);
                continue;
            end
            if nDrug==1
                nFc=nFc+1;
            end
            fprintf('... processing %s\n',File_Name);
            if redo==1 || exist([data_path filesep 'CTET_' File_Name(1:end-4) '_FFT_perBlock_byElec_ICAcleaned.mat'])==0
                subj_pow=[];
                subj_SNR=[];
                subj_SNRtag=[];
                
                load([data_path filesep 'CIcfe_blocks_ft_' File_Name(1:end-4)]);
                
                fprintf('%2.0f-%2.0f\n',0,0)
                for nBl=1:length(data_clean.trial)
                    
                    temp_data=data_clean.trial{nBl}(1:64,data_clean.time{nBl}>0);
                    temp_data=temp_data-repmat(mean(temp_data,2),1,size(temp_data,2));
                    
                    
                    for nEl=1:64
                        fprintf('\b\b\b\b\b\b%2.0f-%2.0f\n',nBl,nEl)
                        block_data=temp_data(nEl,:);
                        
                        %                 w_window=10*data_clean.fsample;
                        %                 w_overlap=w_window/2;
                        %                 df=0.1;
                        %                 freqV=2:0.1:40;
                        %                 [pow,faxis] = pwelch(block_data,w_window,w_overlap,freqV,data_clean.fsample,'psd');
                        
                        param=[];
                        param.method='welch';
                        param.w_window=10*data_clean.fsample;
                        param.w_overlap=param.w_window/2;
                        param.df=0.1;
                        param.freqV=2:0.1:40;
                        param.mindist=0.5;
                        
                        [logSNR, faxis, logpow]=get_logSNR(block_data,data_clean.fsample,param);
                        
                        subj_pow(nBl,nEl,:)=logpow(faxis<40);
                        subj_SNR(nBl,nEl,:)=logSNR(faxis<40);
                        [~,closestidx]=findclosest(faxis,25);
                        subj_SNRtag(nBl,nEl)=logSNR(closestidx);
                        faxis=faxis(faxis<40);
                    end
                end
                %         all_Drugs{nFc}=DrugC;
                %         all_SubInfo(nFc,:)=[SubN SessN];
                %     save([save_path filesep 'SW_detection' filesep 'PH_CTET_SW_' File_Name(1:end-4)],'slow_Waves','hdr','paramSW')
                save([data_path 'CTET_' File_Name(1:end-4) '_FFT_perBlock_byElec_ICAcleaned.mat'],'subj_pow','subj_SNRtag','subj_SNR','faxis');
            else
                if exist('data_clean')==0
                    load([data_path filesep 'CIcfe_blocks_ft_' File_Name(1:end-4)]);
                end
                load([data_path 'CTET_' File_Name(1:end-4) '_FFT_perBlock_byElec_ICAcleaned.mat']); %,'subj_pow','subj_SNRtag','subj_SNR','faxis');
            end
            %             all_pow(nFc,nDrug,:,:,:)=subj_pow;
            %             all_SNR(nFc,nDrug,:,:,:)=subj_SNR;
            %             all_SNRtag(nFc,nDrug,:,:)=subj_SNRtag;
      alphaFreqs=find((faxis>8 & faxis<11));
      
            table_avSW.SSVEP(table_avSW.SubID==SubN & table_avSW.SessN==SessN)=subj_SNRtag(:,match_str(data_clean.label,'Oz'));
            table_avSW.Alpha(table_avSW.SubID==SubN & table_avSW.SessN==SessN)=squeeze(mean(subj_pow(:,match_str(data_clean.label,'Oz'),alphaFreqs),3));

            
            for nE=1:size(subj_SNRtag,2)
                if sum(table_SW.SubID==SubN & table_SW.SessN==SessN & ismember(table_SW.Elec,data_clean.label{nE}))~=0
                    table_SW.SSVEP(table_SW.SubID==SubN & table_SW.SessN==SessN & ismember(table_SW.Elec,data_clean.label{nE}))=subj_SNRtag(:,nE);
                    table_SW.Alpha(table_SW.SubID==SubN & table_SW.SessN==SessN & ismember(table_SW.Elec,data_clean.label{nE}))=squeeze(mean(subj_pow(:,nE,alphaFreqs),3));
                end
            end
        elseif nFc~=0
            %             all_pow(nFc,nDrug,:,:,:)=nan(10,64,640);
%             all_SNR(nFc,nDrug,:,:,:)=nan(10,64,640);
%             all_SNRtag(nFc,nDrug,:,:)=nan(10,64);
        end
    end
end

%%
writetable(table_SW,[save_path filesep 'CTET_SWdetection_thr90_byE_P2P_behav_TF_vec_v6.txt']);

