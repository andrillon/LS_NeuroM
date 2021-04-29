%%
clear all;
% close all;
run ../localdef.m

addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_fooof))
addpath(genpath(path_LSCPtools));

files=dir([data_path filesep 'CIcfe_ft_*.mat']);

% table=readtable([save_path 'CTET_behav_res.txt']);
% load('../ICA_Artifacts.mat')

%% Layout
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

load ../EEG_Cleaning.mat
redo=1;
%%
if redo==1
    nFc=0;
    all_ERP_NT=[];
    all_ERP_TG=[];
    all_DrugC=[];
    for nF=1:length(files)
        File_Name=files(nF).name;
        File_Name=File_Name(10:end);
        septag=findstr(File_Name,'_');
        SubN=str2num(File_Name(1:septag(1)-1));
        SessN=str2num(File_Name(septag(3)-1));
        DrugC=(File_Name(septag(3)+1:septag(3)+3));
        
        table=readtable([save_path 'CTET_behav_' File_Name(1:end-4) '.txt']);
        NT_ITI(nF)=mode(table.ITI);
        TG_ITI(nF)=mode(table.ITI(table.ITI>mode(table.ITI)+0.1));
        
        if length(unique(table.BlockN))~=10
            continue;
        end
        nFc=nFc+1;
        fprintf('... processing %s\n',File_Name);
        load([data_path filesep 'CIcfe_ft_' File_Name(1:end-4)]);
        
        cfg=[];
        cfg.reref           = 'yes';
        cfg.refchannel      = {'M1','M2'};
        cfg.demean          = 'yes';
        cfg.baselinewindow  = [-0.2 0];
        data_clean = ft_preprocessing(cfg,data_clean);
        
        thisF=match_str(EEGCleaning.File,['fe_ft_' File_Name(1:end-4) '.mat']);
        eval(['badTrials=[' EEGCleaning.Trials{thisF} '];']);
        table(ismember(table.TrialN,badTrials),:)=[];
        table.StimType(find(diff(table.BlockN)==1 | diff(table.StimType)==1))=NaN;
        
        % COMPUTE BOTH ONSET AND OFFSET ERP
        cfgerp        = [];
        %     cfgerp.latency        = [-0.2 0.8];
        cfgerp.trials = find(table.StimType==0); %cfgerp.trials(cfgerp.trials>length(data_clean.trial))=[];
        av_data_NT = ft_timelockanalysis(cfgerp, data_clean);
        cfgerp.trials = find(table.StimType==1); %cfgerp.trials(cfgerp.trials>length(data_clean.trial))=[];
        av_data_TG = ft_timelockanalysis(cfgerp, data_clean);
        
        ERP_NT=av_data_NT.avg; %(:,av_data_NT.time>-0.2 & av_data_NT.time<0.8);
        %     [~,idx]=findclosest(av_data_TG.time,NT_ITI(nF)-0.2);
        %     ERP_NT=ERP_NT(:,idx:idx+176);
        %     ERP_NT=ERP_NT-repmat(mean(ERP_NT(:,1:50),2),[1 size(ERP_NT,2)]);
        all_ERP_NT(nFc,:,:)=ERP_NT;
        
        ERP_TG=av_data_TG.avg; %(:,av_data_TG.time>-0.2 & av_data_TG.time<0.8);
        %     [~,idx]=findclosest(av_data_TG.time,TG_ITI(nF)-0.2);
        %     ERP_TG=ERP_TG(:,idx:idx+176);
        %     ERP_TG=ERP_TG-repmat(mean(ERP_TG(:,1:50),2),[1 size(ERP_TG,2)]);
        all_ERP_TG(nFc,:,:)=ERP_TG;
        
        all_DrugC{nFc}=DrugC;
        %     plot(av_data_NT.time,av_data_NT.avg(match_str(av_data_NT.label,'Oz'),:)','r')
        %     plot(av_data_TG.time,av_data_TG.avg(match_str(av_data_TG.label,'Oz'),:)','k')
        
        %     cfgerp           = [];
        %     cfgerp.trials   = find(table.StimType==0);
        %     av_data_NT = ft_timelockanalysis(cfgerp, data_clean);
        
        cfgerp2        = [];
        cfgerp2.trials = find(table.StimType==0)+1; cfgerp2.trials(cfgerp2.trials>length(data_clean.trial))=[];
        av_data_NT_offset = ft_timelockanalysis(cfgerp2, data_clean);
        cfgerp2.trials = find(table.StimType==1)+1; cfgerp2.trials(cfgerp2.trials>length(data_clean.trial))=[];
        av_data_TG_offset = ft_timelockanalysis(cfgerp2, data_clean);
        
        ERP_NT_offset=av_data_NT_offset.avg(:,av_data_NT_offset.time>-0.2 & av_data_NT_offset.time<0.8)-repmat(mean(av_data_NT_offset.avg(:,av_data_NT_offset.time>-0.2 & av_data_NT_offset.time<0),2),1,sum(av_data_NT_offset.time>-0.2 & av_data_NT_offset.time<0.8));
        all_ERP_NT_offset(nFc,:,:)=ERP_NT_offset;
        
        ERP_TG_offset=av_data_TG_offset.avg(:,av_data_TG_offset.time>-0.2 & av_data_TG_offset.time<0.8)-repmat(mean(av_data_TG_offset.avg(:,av_data_TG_offset.time>-0.2 & av_data_TG_offset.time<0),2),1,sum(av_data_TG_offset.time>-0.2 & av_data_TG_offset.time<0.8));
        all_ERP_TG_offset(nFc,:,:)=ERP_TG_offset;
    end
    
    xTime=av_data_NT.time;
    chLabels=av_data_NT.label;
    xTime_offset=av_data_TG_offset.time(av_data_TG_offset.time>-0.2 & av_data_TG_offset.time<0.8);
    save([data_path filesep 'CIcfe_ERPav_TG_NT'],'all_DrugC','all_ERP_TG','all_ERP_NT','xTime','NT_ITI','TG_ITI','chLabels',...
        'all_ERP_NT_offset','all_ERP_NT_offset','xTime_offset')
    
else
    load([data_path filesep 'CIcfe_ERPav_TG_NT'])
    
end

%%

thisCh=match_str(chLabels,'Oz');

figure;
plot(xTime,squeeze(nanmean(all_ERP_NT(:,thisCh,:),1)),'k')
hold on;
plot(xTime,squeeze(nanmean(all_ERP_TG(:,thisCh,:),1)),'r')

all_DrugC2=all_DrugC;
all_DrugC2(NT_ITI>1)={'wrongITI'};

figure;
set(gcf,'Position',[1         378        1396         420]);
% ColorsD={[1 1 1]*0.5,[84 81 153]/255,[170 75 21]/255,[75 128 90]/255};
% ColorsDlabels={'PLA','ATM','MPH','CIT'};
for nDrug=1:4
    subplot(1,4,nDrug);
    hp=[];
    [~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_NT(match_str(all_DrugC2,ColorsDlabels{nDrug}),thisCh,:)),0,'k',0,'-',0.5,1,5,1,3);
    hold on;
    [~,hp(2)]=simpleTplot(xTime,squeeze(all_ERP_TG(match_str(all_DrugC2,ColorsDlabels{nDrug}),thisCh,:)),0,'r',0,'-',0.5,1,5,1,3);
    title(ColorsDlabels{nDrug})
    ylim([-2 10])
    xlim([-0.2 1.8])
    format_fig;
    xlabel('Time from Onset')
    ylabel('ERP (\muV)')
end

%%
figure;
set(gcf,'Position',[1         378        1396         420]);
% ColorsD={[1 1 1]*0.5,[84 81 153]/255,[170 75 21]/255,[75 128 90]/255};
% ColorsDlabels={'PLA','ATM','MPH','CIT'};
    hp=[];
for nDrug=1:4
    [~,hp(nDrug)]=simpleTplot(xTime,squeeze(all_ERP_TG(match_str(all_DrugC2,ColorsDlabels{nDrug}),thisCh,:)-...
        all_ERP_NT(match_str(all_DrugC2,ColorsDlabels{nDrug}),thisCh,:)),0,Colors(nDrug,:),[2 0.05 0.0001 1000],'-',0.5,1,5,1,3);
    title(ColorsDlabels{nDrug})
%     ylim([-2 10])
    xlim([-0.2 1.8])
    format_fig;
    xlabel('Time from Onset')
    ylabel('ERP (\muV)')
end
legend(hp,ColorsDlabels)
%%
winTime1=[0.8 1.15];
winTime2=[1.25 1.5];

cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel=chLabels;
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

figure;
subplot(1,2,1); format_fig;
topo_TG_win1=squeeze(nanmean(nanmean(all_ERP_TG(:,:,xTime>winTime1(1) & xTime<winTime1(2)),1),3));
% simpleTopoPlot_ft(topo_TG_win1(1:64)', layout,'labels',[],0,1);
% colorbar;

% subplot(2,2,1); format_fig;
topo_NT_win1=squeeze(nanmean(nanmean(all_ERP_NT(:,:,xTime>winTime1(1) & xTime<winTime1(2)),1),3));
simpleTopoPlot_ft(topo_TG_win1(1:64)'-topo_NT_win1(1:64)', layout,'labels',[],0,1);
colorbar;
caxis([-1 1]*3)

subplot(1,2,2); format_fig;
topo_TG_win2=squeeze(nanmean(nanmean(all_ERP_TG(:,:,xTime>winTime2(1) & xTime<winTime2(2)),1),3));
% simpleTopoPlot_ft(topo_TG_win2(1:64)', layout,'labels',[],0,1);
% colorbar;
% 
% subplot(2,2,3); format_fig;
topo_NT_win2=squeeze(nanmean(nanmean(all_ERP_NT(:,:,xTime>winTime2(1) & xTime<winTime2(2)),1),3));
simpleTopoPlot_ft(topo_TG_win2(1:64)'-topo_NT_win2(1:64)', layout,'labels',[],0,1);
colorbar;
caxis([-1 1]*5)

%%
figure;
for nDrug=1:4
    subplot(4,2,1+2*(nDrug-1)); format_fig;
    topo_TG_win1=squeeze(nanmean(nanmean(all_ERP_TG(match_str(all_DrugC2,ColorsDlabels{nDrug}),:,xTime>winTime1(1) & xTime<winTime1(2)),1),3));
    topo_NT_win1=squeeze(nanmean(nanmean(all_ERP_NT(match_str(all_DrugC2,ColorsDlabels{nDrug}),:,xTime>winTime1(1) & xTime<winTime1(2)),1),3));
    simpleTopoPlot_ft(topo_TG_win1(1:64)'-topo_NT_win1(1:64)', layout,'on',[],0,1);
    colorbar;
 caxis([-1 1]*3)
   
    subplot(4,2,2+2*(nDrug-1)); format_fig;
    topo_TG_win2=squeeze(nanmean(nanmean(all_ERP_TG(match_str(all_DrugC2,ColorsDlabels{nDrug}),:,xTime>winTime2(1) & xTime<winTime2(2)),1),3));
    topo_NT_win2=squeeze(nanmean(nanmean(all_ERP_NT(match_str(all_DrugC2,ColorsDlabels{nDrug}),:,xTime>winTime2(1) & xTime<winTime2(2)),1),3));
    simpleTopoPlot_ft(topo_TG_win2(1:64)'-topo_NT_win2(1:64)', layout,'on',[],0,1);
    colorbar;
caxis([-1 1]*5)
end