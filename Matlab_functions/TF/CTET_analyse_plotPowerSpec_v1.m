%%
clear all
close all
run ../localdef.m

addpath(path_fieldtrip);
ft_defaults;
ft_defaults;
addpath(genpath(path_LSCPtools));

% data_path='/Volumes/shared/R-MNHS-SPP/Bellgrove-data/Jess Barnes EEG Backup Data/EEG_CTET/';
files=dir([data_path filesep 'CIcfe_ft_*.mat']);

%%
res_mat=[];
drug_cond=[];
myfreqs=1:0.1:40;
nFc=0;
redo=0;

all_pow=[];
all_SNR=[];
all_SNRtag=[];

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
            all_pow(nFc,nDrug,:,:,:)=subj_pow;
            all_SNR(nFc,nDrug,:,:,:)=subj_SNR;
            all_SNRtag(nFc,nDrug,:,:)=subj_SNRtag;
        elseif nFc~=0
            all_pow(nFc,nDrug,:,:,:)=nan(10,64,640);
            all_SNR(nFc,nDrug,:,:,:)=nan(10,64,640);
            all_SNRtag(nFc,nDrug,:,:)=nan(10,64);
        end
    end
end

%% Topographies 25Hz tag
chLabels=data_clean.label(1:64);
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel=chLabels(~ismember(chLabels,{'Iz','M1','M2'}));
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);
correspCh=[];
for nCh=1:length(layout.label)-2
    correspCh(nCh)=match_str(chLabels,layout.label(nCh));
end

%%
figure; set(gcf,'Position',[1     1   320 320]); %575 300
[ha pos]=tight_subplot(2,2,0.05,0.05,0.05);
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
cmap(cmap<0)=0;
OrderPos=[1 3 2 4];
for nD=1:4
      hs=subplot(2,2,nD); format_fig;
            set(hs,'Position',pos{nD})
    temp_topo=squeeze(nanmean(nanmean(all_SNRtag(:,OrderPos(nD),:,correspCh),1),3));
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    colormap(cmap);
if nD==4
        hb=colorbar('Position',[0.9233    0.4133    0.0568    0.1567]);
end
    caxis([0 1]*3);
%     title(ColorsDlabels{OrderPos(nD)});
end
print('-dpng', '-r300', '../../Figures/Topo_FreqTag_v5.png')

%%
ColorRec=[255 255 0]/256;
cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);

figure; set(gcf,'Position',[213         173        320*2         320]);
% subplot(1,4,1);
jbfill([24.5 25.5],[-.7 -.7],[2 2],ColorRec,ColorRec,1,0.2);
format_fig;
hold on;
hp=[];
for nD=1:4
    [~,hp(nD)]=simpleTplot(faxis',squeeze(nanmean(all_SNR(:,OrderPos(nD),:,match_str(chLabels,'Cz'),:),3)),0,ColorsD{OrderPos(nD)},[0 0.05 0.0001 1000],'-',0.1,1,0,1,1);
end
xlim([2 30])
ylim([-.7 2])
xlabel('Frequency (Hz)')
ylabel('SNR')
set(gca,'FontSize',22)
set(gca,'LineWidth',2);
hl=legend(hp,ColorsDlabels(OrderPos),'box','off','position',[0.1750    0.6234    0.1359    0.3172]);
print('-dpng', '-r300', '../../Figures/Topo_FreqTag_Clusters_byFreq_v5.png')

%%
PosDrugs={[3 1],[2 1],[4 1];[],[2 3],[4 3];[],[],[4 2]};
cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);
figure; set(gcf,'Position',[1 1 320*3 320*1]);
[ha pos]=tight_subplot(1,3,0.02,0.05,0.05);
winTime=[0.05 0.3];
for nD=2:4
%     for nD2=1:size(PosDrugs,2)
%         if isempty(PosDrugs{nD,nD2})
%             hs=subplot(3,3,3*(nD-1)+(nD2));
%             set(hs,'Position',pos{3*(nD-1)+(nD2)})
%             set(gcf,'Color','w')
%             set(gca,'Xcolor','w','Ycolor','w')
%             continue;
%         end
 hs=subplot(1,3,nD-1); format_fig;
        set(hs,'Position',pos{nD-1})
        temp_topo=[];
    for nCh=1:length(layout.label)-2
        temp1=squeeze(nanmean(all_SNRtag(:,OrderPos(nD),:,match_str(chLabels,layout.label(nCh))),3));
        temp0=squeeze(nanmean(all_SNRtag(:,1,:,match_str(chLabels,layout.label(nCh))),3));
        [h, pV, ~ , stats]=ttest(temp1,temp0);
        temp_topo(nCh)=stats.tstat;%-...
        %             squeeze(nanmean(nanmean(nanmean(all_ERP_NT_offset(:,nD,match_str(chLabels,layout.label(nCh)),xTime_offset>winTime(1) & xTime_offset<winTime(2)),1),2),4));
    end
    temp_topo1=squeeze(nanmean(all_SNRtag(:,OrderPos(nD),:,correspCh),3));
    temp_topo2=squeeze(nanmean(all_SNRtag(:,1,:,correspCh),3));
    [stat] = compute_Topo_clusterPerm_v2(temp_topo1,temp_topo2,0,chLabels(correspCh),data_clean.fsample,0.05,0.05,1000,layout);
    
    
    simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
    colormap(cmap2);
    
    if isfield(stat,'posclusters') && ~isempty(stat.posclusters)
        sigClusters=find([stat.posclusters.prob]<0.05);
        for k=1:length(sigClusters)
            ft_plot_lay_me(layout, 'chanindx',find(stat.posclusterslabelmat==sigClusters(k)),'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
        end
    end
    if isfield(stat,'negclusters') && ~isempty(stat.negclusters)
        sigClusters=find([stat.negclusters.prob]<0.05);
        for k=1:length(sigClusters)
            ft_plot_lay_me(layout, 'chanindx',find(stat.negclusterslabelmat==sigClusters(k)),'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
        end
    end
    if nD==4
        hb=colorbar('Position',[0.975    0.33   0.02    0.33]);
    end
    caxis([-1 1]*5);
%     title(sprintf('%s vs %s',ColorsDlabels{PosDrugs{nD,nD2}(1)},ColorsDlabels{PosDrugs{nD,nD2}(2)}));
end

print('-dpng', '-r300', '../../Figures/Topo_FreqTag_Clusters_v5.png')


%% Alpha 8-11Hz
figure; set(gcf,'Position',[1     1   320 320]); %575 300
[ha pos]=tight_subplot(2,2,0.05,0.05,0.05);
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
cmap(cmap<0)=0;
for nD=1:4
    hs=subplot(2,2,nD); format_fig;
    set(hs,'Position',pos{nD})
    temp_topo=squeeze(nanmean(nanmean(nanmean(all_pow(:,OrderPos(nD),:,correspCh,faxis>8 & faxis<11),1),3),5));
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    colormap(cmap);
    if nD==4
        hb=colorbar('Position',[0.9233    0.4133    0.0568    0.1567]);
    end
    caxis([-2.5 -1]);
    %     title(ColorsDlabels{nD});
end
print('-dpng', '-r300', '../../Figures/Topo_Alpha_v5.png')

%%
% alphaFreqs=find((faxis>8 & faxis<8.6) | (faxis>8.95 & faxis<9.8) | (faxis>10.125 & faxis<11));
alphaFreqs=find((faxis>8 & faxis<11));
figure; set(gcf,'Position',[213         173        320*2         320]);
% subplot(1,4,1); 
format_fig;
jbfill([8 11],[-5 -5],[-.5 -.5],ColorRec,ColorRec,1,0.2);
hold on;
hp=[];
for nD=1:4
    [~,hp(nD)]=simpleTplot(faxis',squeeze(nanmean(all_pow(:,OrderPos(nD),:,match_str(chLabels,'Cz'),:),3)),0,ColorsD{OrderPos(nD)},[0 0.05 0.0001 1000],'-',0.1,1,0,1,1);
end
xlim([2 30])
ylim([[-5 -0.5]])
xlabel('Frequency (Hz)')
ylabel('PSD (log)')

set(gca,'FontSize',22)
set(gca,'LineWidth',2);
hl=legend(hp,ColorsDlabels(OrderPos),'box','off','position',[0.5016    0.5984    0.1359    0.3172]);
print('-dpng', '-r300', '../../Figures/Topo_Alpha_Clusters_byFreq_v5.png')

%%
figure; set(gcf,'Position',[1 1 320*3 320*1]);
[ha pos]=tight_subplot(1,3,0.02,0.05,0.05);
winTime=[0.05 0.3];
for nD=2:4 %size(PosDrugs,1)
%     for nD2=1:size(PosDrugs,2)
%         if isempty(PosDrugs{nD,nD2})
%             hs=subplot(3,3,3*(nD-1)+(nD2));
%             set(hs,'Position',pos{3*(nD-1)+(nD2)})
%             set(gcf,'Color','w')
%             set(gca,'Xcolor','w','Ycolor','w')
%             continue;
%         end
        hs=subplot(1,3,nD-1); format_fig;
        set(hs,'Position',pos{nD-1})
        
    temp_topo=[];
    for nCh=1:length(layout.label)-2
        temp1=squeeze(nanmean(nanmean(all_pow(:,OrderPos(nD),:,match_str(chLabels,layout.label(nCh)),alphaFreqs),3),5));
        temp0=squeeze(nanmean(nanmean(all_pow(:,1,:,match_str(chLabels,layout.label(nCh)),alphaFreqs),3),5));
        [h, pV, ~ , stats]=ttest(temp1,temp0);
        temp_topo(nCh)=stats.tstat;%-...
        %             squeeze(nanmean(nanmean(nanmean(all_ERP_NT_offset(:,nD,match_str(chLabels,layout.label(nCh)),xTime_offset>winTime(1) & xTime_offset<winTime(2)),1),2),4));
    end
    temp_topo1=squeeze(nanmean(nanmean(all_pow(:,OrderPos(nD),:,correspCh,alphaFreqs),3),5));
    temp_topo2=squeeze(nanmean(nanmean(all_pow(:,1,:,correspCh,alphaFreqs),3),5));
    [stat] = compute_Topo_clusterPerm_v2(temp_topo1,temp_topo2,0,chLabels(correspCh),data_clean.fsample,0.05,0.05,1000,layout);
    
    
    simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
    colormap(cmap2);
    
    if isfield(stat,'posclusters') && ~isempty(stat.posclusters)
        sigClusters=find([stat.posclusters.prob]<0.05);
        for k=1:length(sigClusters)
            ft_plot_lay_me(layout, 'chanindx',find(stat.posclusterslabelmat==sigClusters(k)),'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
        end
    end
    if isfield(stat,'negclusters') && ~isempty(stat.negclusters)
        sigClusters=find([stat.negclusters.prob]<0.05);
        for k=1:length(sigClusters)
            ft_plot_lay_me(layout, 'chanindx',find(stat.negclusterslabelmat==sigClusters(k)),'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
        end
    end
    if nD==4
        hb=colorbar('Position',[0.975    0.33   0.02    0.33]);
    end
    caxis([-1 1]*5);
%     title(sprintf('%s vs %s',ColorsDlabels{PosDrugs{nD,nD2}(1)},ColorsDlabels{PosDrugs{nD,nD2}(2)}));
%     end
end

print('-dpng', '-r300', '../../Figures/Topo_Alpha_Clusters_v5.png')


%%
[~,closestindex]=findclosest(faxis,26);
all_pow2=all_pow-repmat(squeeze(all_pow(:,:,:,:,closestindex)),[1 1 1 1 length(faxis)]);

thisCh='Fz';
figure;
% subplot(1,4,1); 
format_fig;
hold on;
hp=[];
for nD=1:4
    [~,hp(nD)]=simpleTplot(faxis',squeeze(nanmean(all_pow(:,nD,:,match_str(chLabels,thisCh),:),3)),0,Colors(nD,:),[0 0.05 0.0001 1000],'-',0.1,1,0,1,1);
end
xlim([2 40])
% ylim([[-5 -0.5]])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
legend(hp,ColorsDlabels)

% figure; set(gcf,'Position',[213         173        1027/4         805/3]);
% 
% print('-dpng', '-r300', '../../Figures/Cz_PowSpec_v5.png')
% 
% figure; set(gcf,'Position',[213         173        1027/4         805/3]);
% 
% % ylim([[-5 -0.5]])