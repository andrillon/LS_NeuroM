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
files=dir([data_path filesep '*.bdf']);

table=readtable('/Users/tand0009/Data/CTET_Dockree/CTET_behav_res.txt');
%%
res_mat=[];
drug_cond=[];
myfreqs=1:0.1:40;
nFc=0;
redo=1; complete=0;

for nF=1:length(files)
    File_Name=files(nF).name;
    fprintf('... processing %s\n',File_Name);
    septag=findstr(File_Name,'_');
    SubN=str2num(File_Name(1:septag(1)-1));
    SessN=str2num(File_Name(septag(3)-1));
    DrugC=(File_Name(septag(3)+1:septag(3)+3));
    hdr_ori=ft_read_header([data_path filesep File_Name]);
    if hdr_ori.Fs~=1024 || length(unique(table.BlockN(table.SubID==SubN & table.SessN==SessN)))~=10
        continue;
    end
    nFc=nFc+1;
    if redo==1 || exist(['/Users/tand0009/Data/CTET_Dockree/CTET_' File_Name(1:end-4) '_FOOOF_FFT_perBlock_byElec_avMast_ICAcleaned.mat'])==0
        all_peaks=[];
        all_bg=[];
        all_pow=[];
        all_spec=[];
        all_fspec=[];
        all_fit=[];
        
        load([data_path filesep 'CIcfe_blocks_ft_' File_Name(1:end-4)]);
        if length(data_clean.trial)~=10
            continue;
        end
        
        fprintf('%2.0f-%2.0f\n',0,0)
        for nBl=1:length(data_clean.trial)
            
            temp_data=data_clean.trial{nBl}(1:64,data_clean.time{nBl}>0);
            temp_data=temp_data-repmat(mean(temp_data(match_str(data_clean.label,{'TP7','TP8'}),:),1),size(temp_data,1),1);
            temp_data=temp_data-repmat(mean(temp_data,2),1,size(temp_data,2));
            
            
            for nEl=1:64
                fprintf('\b\b\b\b\b\b%2.0f-%2.0f\n',nBl,nEl)
                block_data=temp_data(nEl,:);
                
                w_window=10*data_clean.fsample;
                w_overlap=w_window/2;
                df=0.2;
                freqV=1:0.1:50;
                [pow,faxis] = pwelch(block_data,w_window,w_overlap,freqV,data_clean.fsample,'psd');
                
                f_range = [2, 50];
                settings = struct();  % Use defaults
                try
                    fooof_results = fooof(faxis, pow, f_range, settings,1);
                    
                    myfreqs=fooof_results.freqs;
                    
                    all_pow(nBl,nEl,:)=log10(pow);
                    all_peaks=[all_peaks ; [repmat([SubN SessN nFc nBl nEl],size(fooof_results.peak_params,1),1) fooof_results.peak_params]];
                    all_bg=[all_bg ; [SubN SessN nFc nBl nEl fooof_results.background_params]];
                    
                    all_spec(nBl,nEl,:)=fooof_results.power_spectrum;
                    all_fspec(nBl,nEl,:)=fooof_results.fooofed_spectrum;
                    all_fit(nBl,nEl,:)=fooof_results.ap_fit;
                catch
                    all_pow(nBl,nEl,:)=nan(1,length(faxis));
                    all_bg=[all_bg ; [SubN SessN nFc nBl nEl nan(1,2)]];
                    
                    all_spec(nBl,nEl,:)=nan(1,length(myfreqs));
                    all_fspec(nBl,nEl,:)=nan(1,length(myfreqs));
                    all_fit(nBl,nEl,:)=nan(1,length(myfreqs));
                end
            end
        end
        %         all_Drugs{nFc}=DrugC;
        %         all_SubInfo(nFc,:)=[SubN SessN];
        %     save([save_path filesep 'SW_detection' filesep 'PH_CTET_SW_' File_Name(1:end-4)],'slow_Waves','hdr','paramSW')
        save(['/Users/tand0009/Data/CTET_Dockree/CTET_' File_Name(1:end-4) '_FOOOF_FFT_perBlock_byElec_avMast_ICAcleaned.mat'],'all_pow','myfreqs','all_peaks','all_bg','DrugC',...
            'all_spec','all_fspec','all_fit');
    end
end



%%
% thiCh=match_str(hdr.label,'Pz');
% figure;
%     hold on;
% Drugs={'PLA','MPH','CIT','ATM'};
% for nDrug=1:4
% [~,hp(nDrug)]=simpleTplot(myfreqs,squeeze(mean(all_pow(ismember(all_Drugs,Drugs{nDrug}),:,thiCh,:),2)),0,Colors(nDrug,:),0,'-',.5,1,10,[],2);
% end
%  xlim([2 20])
% legend(hp,Drugs);
% %%
% table_bg=array2table(all_bg,'VariableNames',{'SubID','SessN','nFile','BlockN','ElecN','offset','slope'});
% table_bg.SubID=categorical(table_bg.SubID);
% table_bg.Drug=cell(size(table_bg,1),1);
% nFiles=unique(table_bg.nFile);
% for k=1:length(nFiles)
%     table_bg.Drug(table_bg.nFile==nFiles(k))=repmat(all_Drugs(nFiles(k)),sum(table_bg.nFile==nFiles(k)),1);
% end
% table_bg.Drug=categorical(table_bg.Drug);
%
% table_bg.Elec=cell(size(table_bg,1),1);
% for k=1:64
%     table_bg.Elec(table_bg.ElecN==k)=repmat(hdr.label(k),sum(table_bg.ElecN==k),1);
% end
% table_bg.Elec=categorical(table_bg.Elec);
% table_bg.Drug=reordercats(table_bg.Drug,[4 1 2 3]);
%
% writetable(table_bg,'/Users/tand0009/Data/CTET_Dockree/CTET_FOOOF_Background_perBlock_byElec_avMast.txt');
%
% %%21
% table_pk=array2table(all_peaks,'VariableNames',{'SubID','SessN','nFile','BlockN','ElecN','peak_freq','peak_amp','peak_width'});
% table_pk.SubID=categorical(table_pk.SubID);
% uS=unique(table_pk.SubID);
% table_pk.Drug=cell(size(table_pk,1),1);
% nFiles=unique(table_pk.nFile);
% for k=1:length(nFiles)
%     table_pk.Drug(table_pk.nFile==nFiles(k))=repmat(all_Drugs(nFiles(k)),sum(table_pk.nFile==nFiles(k)),1);
% end
% table_pk.Drug=categorical(table_pk.Drug);
%
% table_pk.Elec=cell(size(table_pk,1),1);
% for k=1:64
%     table_pk.Elec(table_pk.ElecN==k)=repmat(hdr.label(k),sum(table_pk.ElecN==k),1);
% end
% table_pk.Elec=categorical(table_pk.Elec);
% table_pk.Drug=reordercats(table_pk.Drug,[4 1 2 3]);
%
% writetable(table_pk,'/Users/tand0009/Data/CTET_Dockree/CTET_FOOOF_Peaks_perBlock_byElec_avMast.txt');
%
% %%
% mdl0= fitlme(table_bg,sprintf('slope~1+BlockN+Elec+(1|SubID)'));
% mdl1= fitlme(table_bg,sprintf('slope~1+BlockN+Elec+Drug+(1|SubID)'));
% mdl2= fitlme(table_bg,sprintf('slope~1+BlockN+Elec*Drug+(1|SubID)'));
% mdl3= fitlme(table_bg,sprintf('slope~1+BlockN*Elec*Drug+(1|SubID)'));
%
% %%
% table_pkap=table_pk(table_pk.peak_freq>8 & table_pk.peak_freq<12,:);
%
% mdl0= fitlme(table_pkap,sprintf('peak_freq~1+BlockN+Elec+(1|SubID)'));
% mdl1= fitlme(table_pkap,sprintf('peak_freq~1+BlockN+Elec+Drug+(1|SubID)'));
% mdl2= fitlme(table_pkap,sprintf('peak_freq~1+BlockN+Elec*Drug+(1|SubID)'));
% mdl3= fitlme(table_pkap,sprintf('peak_freq~1+BlockN*Elec*Drug+(1|SubID)'));

% %%
% filestomove=dir(['/Users/tand0009/Data/CTET_Dockree/CTET_*_FOOOF_FFT_perBlock_byElec_avMast.txt']);
% for nF=1:length(filestomove)
% movefile([filestomove(nF).folder filesep filestomove(nF).name],[filestomove(nF).folder filesep filestomove(nF).name(1:end-4) '.mat']);
% end