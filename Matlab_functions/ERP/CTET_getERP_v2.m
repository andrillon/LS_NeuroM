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
redo=0;
%%

if redo==1
    nFc=0;
    all_ERP_NT=[];
    all_ERP_TG=[];
    all_DrugC=[];
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
                
                table=readtable([save_path 'CTET_behav_' File_Name(1:end-4) '.txt']);
                
                if nDrug==1
                    nFc=nFc+1;
                end
                fprintf('... processing %s\n',File_Name);
                load([data_path filesep 'CIcfe_ft_' File_Name(1:end-4)]);
                
                %         cfg=[];
                %         cfg.reref           = 'yes';
                %         cfg.refchannel      = 'all';
                %         cfg.demean          = 'yes';
                %         cfg.baselinewindow  = [-0.2 0];
                %         data_clean = ft_preprocessing(cfg,data_clean);
                
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
                all_ERP_NT(nFc,nDrug,:,:)=ERP_NT;
                
                ERP_TG=av_data_TG.avg; %(:,av_data_TG.time>-0.2 & av_data_TG.time<0.8);
                %     [~,idx]=findclosest(av_data_TG.time,TG_ITI(nF)-0.2);
                %     ERP_TG=ERP_TG(:,idx:idx+176);
                %     ERP_TG=ERP_TG-repmat(mean(ERP_TG(:,1:50),2),[1 size(ERP_TG,2)]);
                all_ERP_TG(nFc,nDrug,:,:)=ERP_TG;
                
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
                all_ERP_NT_offset(nFc,nDrug,:,:)=ERP_NT_offset;
                
                ERP_TG_offset=av_data_TG_offset.avg(:,av_data_TG_offset.time>-0.2 & av_data_TG_offset.time<0.8)-repmat(mean(av_data_TG_offset.avg(:,av_data_TG_offset.time>-0.2 & av_data_TG_offset.time<0),2),1,sum(av_data_TG_offset.time>-0.2 & av_data_TG_offset.time<0.8));
                all_ERP_TG_offset(nFc,nDrug,:,:)=ERP_TG_offset;
            else
                all_ERP_NT(nFc,nDrug,:,:)=nan(66,589);
                all_ERP_TG(nFc,nDrug,:,:)=nan(66,589);
                all_ERP_NT_offset(nFc,nDrug,:,:)=nan(66,256);
                all_ERP_TG_offset(nFc,nDrug,:,:)=nan(66,256);
            end
        end
    end
    
    xTime=av_data_NT.time;
    chLabels=av_data_NT.label;
    xTime_offset=av_data_TG_offset.time(av_data_TG_offset.time>-0.2 & av_data_TG_offset.time<0.8);
    save([data_path filesep 'CIcfe_ERPav_TG_NT_paired'],'all_ERP_TG','all_ERP_NT','xTime','chLabels',...
        'all_ERP_NT_offset','all_ERP_TG_offset','xTime_offset')
    
else
    load([data_path filesep 'CIcfe_ERPav_TG_NT_paired'])
    
end

%%
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

thisCh=match_str(chLabels,'CPz');

figure;
set(gcf,'Position',[1         378        1396         420]);
% ColorsD={[1 1 1]*0.5,[84 81 153]/255,[170 75 21]/255,[75 128 90]/255};
% ColorsDlabels={'PLA','ATM','MPH','CIT'};
for nDrug=1:4
    subplot(1,4,nDrug);
    hp=[];
    [~,hp(1)]=simpleTplot(xTime,squeeze(all_ERP_NT(:,nDrug,thisCh,:)),0,'k',0,'-',0.1,1,5,1,3);
    hold on;
    [~,hp(2)]=simpleTplot(xTime,squeeze(all_ERP_TG(:,nDrug,thisCh,:)),0,'r',0,'-',0.1,1,5,1,3);
    title(ColorsDlabels{nDrug})
    ylim([-2 5.5])
    xlim([-0.2 1.8])
    format_fig;
    xlabel('Time from Onset')
    ylabel('ERP (\muV)')
end


%%
diffERP=squeeze(all_ERP_TG(:,:,thisCh,:)-all_ERP_NT(:,:,thisCh,:));
diffERP_Mat=[]; diffERP_Group=[];
for nD=1:4
diffERP_Mat=[diffERP_Mat ; squeeze(diffERP(:,nD,:))];
diffERP_Group=[diffERP_Group ; nD*ones(size(diffERP(:,nD,:),1),1)];
end
[realpos_lin ~]=get_cluster_permutation_aov(diffERP_Mat,diffERP_Group,...
    0.05,0.1,100,xTime);
    
%%
thisCh=match_str(chLabels,'C4');

figure;
set(gcf,'Position',[1         378        1396         420]);
% ColorsD={[1 1 1]*0.5,[84 81 153]/255,[170 75 21]/255,[75 128 90]/255};
% ColorsDlabels={'PLA','ATM','MPH','CIT'};
    hp=[];
for nDrug=1:4
    [~,hp(nDrug)]=simpleTplot(xTime,squeeze(all_ERP_TG(:,nDrug,thisCh,:)),0,Colors(nDrug,:),[0 0.05 0.0001 1000],'-',0.1,1,5,1,3);
    simpleTplot(xTime,squeeze(all_ERP_NT(:,nDrug,thisCh,:)),0,'k',[0 0.05 0.0001 1000],'-',0.1,1,5,1,3);
    title(ColorsDlabels{nDrug})
%     ylim([-2 10])
    xlim([-0.2 1.8])
    format_fig;
    xlabel('Time from Onset')
    ylabel('ERP (\muV)')
end
legend(hp,ColorsDlabels)


%%

figure;
set(gcf,'Position',[1         134         560         800]);
subplot(2,1,2);
thisCh=match_str(chLabels,'Pz');
jbfill([0.05 0.3],[-1.5 -1.5],[6 6],[50,205,50]/256,[50,205,50]/256,1,0.2);
    hp=[];
for nDrug=1:4
%     [~,hp(nDrug)]=simpleTplot(xTime_offset,squeeze(all_ERP_TG_offset(:,nDrug,thisCh,:)-...
%         all_ERP_NT_offset(:,nDrug,thisCh,:)),0,Colors(nDrug,:),[0 0.05 0.0001 1000],'-',0.1,1,5,1,3);
    [~,hp(nDrug)]=simpleTplot(xTime_offset,squeeze(all_ERP_TG_offset(:,nDrug,thisCh,:)),0,Colors(nDrug,:),[0 0.05 0.0001 1000],'-',0.1,1,5,1,3);
    [~,hp(nDrug)]=simpleTplot(xTime_offset,squeeze(all_ERP_NT_offset(:,nDrug,thisCh,:)),0,Colors(nDrug,:),[0 0.05 0.0001 1000],':',0.1,1,5,1,3);
    title('Pz')
    ylim([-1.5 6])
    xlim([-0.2 0.6])
    format_fig;
    xlabel('Time from Offset')
    ylabel('ERP (\muV)')
end
legend(hp,ColorsDlabels)

subplot(2,1,1);
thisCh=match_str(chLabels,'Fz');
jbfill([0.05 0.3],[-4.5 -4.5],[3 3],[50,205,50]/256,[50,205,50]/256,1,0.2);
    hp=[];
for nDrug=1:4
%     [~,hp(nDrug)]=simpleTplot(xTime_offset,squeeze(all_ERP_TG_offset(:,nDrug,thisCh,:)-...
%         all_ERP_NT_offset(:,nDrug,thisCh,:)),0,Colors(nDrug,:),[0 0.05 0.0001 1000],'-',0.1,1,5,1,3);
    [~,hp(nDrug)]=simpleTplot(xTime_offset,squeeze(all_ERP_TG_offset(:,nDrug,thisCh,:)),0,Colors(nDrug,:),[0 0.05 0.0001 1000],'-',0.1,1,5,1,3);
    [~,hp(nDrug)]=simpleTplot(xTime_offset,squeeze(all_ERP_NT_offset(:,nDrug,thisCh,:)),0,Colors(nDrug,:),[0 0.05 0.0001 1000],':',0.1,1,5,1,3);
    title('Fz')
    ylim([-4.5 3])
    xlim([-0.2 0.6])
    format_fig;
    xlabel('Time from Offset')
    ylabel('ERP (\muV)')
end
% legend(hp,ColorsDlabels)

print('-dpng', '-r300', '../../Figures/TimePlot_ERP_Offset.png')

%% Topographies on 0.9-1.1s post onset
% 
% cfg = [];
% cfg.layout = 'biosemi64.lay';
% cfg.channel=chLabels(~ismember(chLabels,{'Iz','M1','M2'}));
% cfg.center      = 'yes';
% layout=ft_prepare_layout(cfg);
% 
% figure;
% winTime=[0.9 1.1];
% for nD=1:4
%     subplot(1,4,nD); format_fig;
%     temp_topo=[];
%     for nCh=1:length(layout.label)-2
%         temp_topo(nCh)=squeeze(nanmean(nanmean(nanmean(all_ERP_TG(:,nD,match_str(chLabels,layout.label(nCh)),xTime>winTime(1) & xTime<winTime(2)),1),2),4))-...
%              squeeze(nanmean(nanmean(nanmean(all_ERP_TG(:,nD,match_str(chLabels,layout.label(nCh)),xTime>0.56 & xTime<0.64),1),2),4));
%     end
%     simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
%     colorbar;
% %     caxis([-1 1]*4);
%     title(ColorsDlabels{nD});
% end

%% Topographies on 0-0.3s post offset
figure; set(gcf,'Position',[1  1  509         878]);
winTime=[0.05 0.3];
for nD=1:4
    subplot(4,1,nD); format_fig;
    temp_topo=[];
    for nCh=1:length(layout.label)-2
        temp_topo(nCh)=squeeze(nanmean(nanmean(nanmean(all_ERP_TG_offset(:,nD,match_str(chLabels,layout.label(nCh)),xTime_offset>winTime(1) & xTime_offset<winTime(2)),1),2),4));%-...
%             squeeze(nanmean(nanmean(nanmean(all_ERP_NT_offset(:,nD,match_str(chLabels,layout.label(nCh)),xTime_offset>winTime(1) & xTime_offset<winTime(2)),1),2),4));
    end
    simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
if nD==4
        hb=colorbar('Position',[0.85    0.6373    0.025    0.1]);
end
caxis([-1 1]*4);
    title(ColorsDlabels{nD});
end

%%
%% Topographies on 0-0.3s post offset
figure; %set(gcf,'Position',[1  1  509         878]);
winTime=[0.05 0.3];
% for nD=1:4
%     subplot(4,1,nD); 
    format_fig;
    temp_topo=[];
    for nCh=1:length(layout.label)-2
        temp_topo(nCh)=squeeze(nanmean(nanmean(nanmean(all_ERP_TG_offset(:,:,match_str(chLabels,layout.label(nCh)),xTime_offset>winTime(1) & xTime_offset<winTime(2)),1),2),4));%-...
%             squeeze(nanmean(nanmean(nanmean(all_ERP_NT_offset(:,nD,match_str(chLabels,layout.label(nCh)),xTime_offset>winTime(1) & xTime_offset<winTime(2)),1),2),4));
    end
    simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
% if nD==4
        hb=colorbar('Position',[0.8500    0.6373    0.0411    0.2032]);
% end
caxis([-1 1]*4);
    title('ALL');
% end
print('-dpng', '-r300', '../../Figures/Topo_ERP_Offset_ALL.eps')

%%
PosDrugs={[3 1],[2 1],[4 1];[],[2 3],[4 3];[],[],[4 2]};
cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);
figure; set(gcf,'Position',[1 1 880 880]);
winTime=[0.05 0.3];
for nD=1:size(PosDrugs,1)
    for nD2=1:size(PosDrugs,2)
        if isempty(PosDrugs{nD,nD2})
            continue;
        end
    subplot(3,3,3*(nD-1)+(nD2)); format_fig;
    temp_topo=[];
    for nCh=1:length(layout.label)-2
        temp1=squeeze(nanmean(nanmean(all_ERP_TG_offset(:,PosDrugs{nD,nD2}(1),match_str(chLabels,layout.label(nCh)),xTime_offset>winTime(1) & xTime_offset<winTime(2)),2),4));
        temp0=squeeze(nanmean(nanmean(all_ERP_TG_offset(:,PosDrugs{nD,nD2}(2),match_str(chLabels,layout.label(nCh)),xTime_offset>winTime(1) & xTime_offset<winTime(2)),2),4));
        [h, pV, ~ , stats]=ttest(temp1,temp0);
        temp_topo(nCh)=stats.tstat;%-...
        %             squeeze(nanmean(nanmean(nanmean(all_ERP_NT_offset(:,nD,match_str(chLabels,layout.label(nCh)),xTime_offset>winTime(1) & xTime_offset<winTime(2)),1),2),4));
    end
    temp_topo1=squeeze(nanmean(nanmean(all_ERP_TG_offset(:,PosDrugs{nD,nD2}(2),correspCh,xTime_offset>winTime(1) & xTime_offset<winTime(2)),2),4));
    temp_topo2=squeeze(nanmean(nanmean(all_ERP_TG_offset(:,PosDrugs{nD,nD2}(1),correspCh,xTime_offset>winTime(1) & xTime_offset<winTime(2)),2),4));
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
        hb=colorbar('Position',[0.9195    0.6373    0.0143    0.2881]);
    end
    caxis([-1 1]*5);
    title(sprintf('%s vs %s',ColorsDlabels{PosDrugs{nD,nD2}(1)},ColorsDlabels{PosDrugs{nD,nD2}(2)}));
    end
end

print('-dpng', '-r300', '../../Figures/Topo_ERP_Offset_Clusters.eps')

%%
diffERP=squeeze(all_ERP_TG_offset(:,:,thisCh,:)-all_ERP_NT_offset(:,:,thisCh,:));
diffERP_Mat=[]; diffERP_Group=[];
for nD=1:4
diffERP_Mat=[diffERP_Mat ; squeeze(diffERP(:,nD,:))];
diffERP_Group=[diffERP_Group ; nD*ones(size(diffERP(:,nD,:),1),1)];
end
[realpos_lin]=get_cluster_permutation_aov(diffERP_Mat,diffERP_Group,...
    0.05,0.05,1000,xTime_offset);

allF=[];
for nT=1:size(diffERP_Mat,2)
    [p,anovatab,stats]=anova1(diffERP_Mat(:,nT),diffERP_Group,'off');
    allF(nT,1)=anovatab{2,5};
    allF(nT,2)=anovatab{2,6};
end