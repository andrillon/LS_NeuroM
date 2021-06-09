%%
clear all
close all

run ../localdef.m
addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_LSCPtools));
addpath(genpath([pwd filesep '..']));

% table=readtable('/Users/tand0009/Data/CTET_Dockree/CTET_behav_res.txt');
table_SW=readtable([save_path filesep 'CTET_SWdetection_thr90_byE_P2P_behav_vec_v6.txt']);
table_avSW=readtable([save_path filesep 'CTET_SWdetection_thr90_byE_P2P_avDens_behav_vec_v6.txt']);
table_behav=readtable([save_path filesep 'CTET_behav_resblock.txt']);
% prctile(table_SW.SWdens,99)
%%
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel=unique(table_SW.Elec);
cfg.channel(match_str(cfg.channel,'Iz'))=[];
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

cmap=cbrewer('seq','YlOrRd',64,'nearest'); % select a sequential colorscale from yellow to red (64)
Drugs={'PLA','ATM','CIT','MPH'};

limMax=[3.5 5.5];

figure;
temp_topo=[];
for nE=1:length(layout.label)-2
    temp_topo(nE)=nanmean(table_SW.SWdens(find_trials(table_SW.Elec,layout.label{nE})));
end
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title('All');
% caxis(limMax)
h=colorbar;
colormap(cmap);
ylabel(h, 'waves/min')
%     set(h,'Position',[0.85 0.7 0.04 0.2])
format_fig;
% print('-dpng', '-r300', '../../Figures/Topo_LS_SWdens.png')

%% By Block
figure; set(gcf,'Position',[228         128        1363         828]);
limMax=[3.5 6];
for nBl=1:10
    subplot(2,5,nBl)
    sub_table_SW=table_SW(table_SW.BlockN==nBl,:);
    temp_topo=[];
    for nE=1:length(layout.label)-2
        temp_topo(nE)=nanmean(sub_table_SW.SWdens(find_trials(sub_table_SW.Elec,layout.label{nE})));
    end
    simpleTopoPlot_ft(temp_topo', layout,'off',[],0,1);
    title(sprintf('Block %g',nBl));
    caxis(limMax)
    %     h=colorbar;
    colormap(cmap);
    %     ylabel(h, 'waves/min')
    %     set(h,'Position',[0.85 0.7 0.04 0.2])
    if nBl==10
        h=colorbar;
        set(h,'Position',[0.9 0.4 0.02 0.2])
        ylabel(h, 'waves/min')
    end
    format_fig;
end
% print('-dpng', '-r300', '../../Figures/Topo_LS_SWdens_byBlock.png')

figure;
temp_plot=[];
for nBl=1:10
    sub_table_SW=table_SW(table_SW.BlockN==nBl,:);
    temp=grpstats(sub_table_SW.SWdens,sub_table_SW.SubID);
    temp_plot(1,nBl)=nanmean(temp);
    temp_plot(2,nBl)=sem(temp);
end
plot(1:10,temp_plot(1,:),'Color','k','Marker','o','LineWidth',3)
hold on;
errorbar(1:10,temp_plot(1,:),temp_plot(2,:),'Color','k','LineWidth',2)
format_fig;
xlabel('Block')
xlim([0.5 10.5])
ylabel('waves/min')
% print('-dpng', '-r300', '../../Figures/Line_LS_SWdens_byBlock.png')

%% By Drug
limMax=[0 18];
cmap=colormap('hot'); cmap=flipud(cmap);

figure; hold on; %set(gcf,'Position',[201         428        1135         500]);
PosPlots=[1 3 2 4];
for nD=1:4
    subplot(2,2,PosPlots(nD));
    sub_table_SW=table_SW(ismember(table_SW.Drug,Drugs{nD}),:);
    temp_topo=[];
    for nE=1:length(layout.label)-2
        temp_topo(nE)=nanmean(sub_table_SW.SWdens(~cellfun(@isempty,regexp(sub_table_SW.Elec,layout.label{nE}))));
    end
    simpleTopoPlot_ft(temp_topo', layout,'off',[],0,1);
    format_fig;
    caxis(limMax)
    if nD==4
        h=colorbar;
        set(h,'Position',[0.93 0.4 0.02 0.2])
    end
    colormap(cmap);
    title(Drugs{nD})
end
% print('-dpng', '-r300', '../../Figures/Topo_LS_SWdens_byDrug.png')

figure;
for nD=1:4
    temp_block=[];
    for nBl=1:10
        sub_table_SW=table_SW(ismember(table_SW.Drug,Drugs{nD}) & table_SW.BlockN==nBl,:);
        temp=grpstats(sub_table_SW.SWdens,sub_table_SW.SubID);
        temp_block(1,nBl)=nanmean(temp);
        temp_block(2,nBl)=sem(temp);
    end
    plot((1:10)+(nD-2.5)/10,temp_block(1,:),'Color',Colors(nD,:),'Marker','o','LineWidth',3)
    hold on;
    errorbar((1:10)+(nD-2.5)/10,temp_block(1,:),temp_block(2,:),'Color',Colors(nD,:),'LineWidth',2)
end
format_fig;
xlabel('Block')
xlim([0.5 10.5])
ylabel('waves/min')
% print('-dpng', '-r300', '../../Figures/Line_LS_SWdens_byDrugAndBlock.png')


%% Models
table_SW2=table_SW;
table_SW2.SubID=categorical(table_SW2.SubID);
table_SW2.SessN=categorical(table_SW2.SessN);
table_SW2.Elec=categorical(table_SW2.Elec);
table_SW2.Drug=categorical(table_SW2.Drug);
table_SW2.Drug=reordercats(table_SW2.Drug,[4 1 2 3]);

mdl0=fitlme(table_SW2,'SWdens~1+(1|SubID)');
mdl1=fitlme(table_SW2,'SWdens~1+Elec+(1|SubID)');
mdl2=fitlme(table_SW2,'SWdens~1+Elec+BlockN+(1|SubID)');
mdl3=fitlme(table_SW2,'SWdens~1+Elec+BlockN+Drug+(1|SubID)');
mdl4=fitlme(table_SW2,'SWdens~1+Elec*(BlockN+Drug)+(1|SubID)');
mdl5=fitlme(table_SW2,'SWdens~1+Elec*BlockN*Drug+(1|SubID)');

compare(mdl4,mdl5)

%%
redo=1;
totperm=1000;
if redo==1
    temp_topo_tval=[];
    temp_topo_pval=[];
    % fprintf('%2.0f/%2.0f\n',0,64)
    SWdens_est=cell(1,2);
    for nE=1:length(layout.label)-2
        fprintf('%2.0f/%2.0f\n',nE,64)
        sub_table_SW2=table_SW2(match_str(table_SW.Elec,layout.label{nE}),:);
        mdl_byEle{nE}=fitlme(sub_table_SW2,'SWdens~1+BlockN+Drug+(1|SubID)');
        temp_topo_tval(nE,:)=double(mdl_byEle{nE}.Coefficients(2:end,4));
        temp_topo_pval(nE,:)=double(mdl_byEle{nE}.Coefficients(2:end,6));
        
        if nE==1
            [real_out, perm_out, out_pred_perm]=lme_perm_lsneurom(sub_table_SW2,'Drug','SWdens~1+pred+(1|SubID)',totperm);
        else
            [real_out, perm_out]=lme_perm_lsneurom(sub_table_SW2,'Drug','SWdens~1+pred+(1|SubID)',totperm, out_pred_perm);
        end
        SWdens_est{1}=[SWdens_est{1} ; [nE*ones(3,1) real_out]];
        for nDrug=1:3
            SWdens_est{2}=[SWdens_est{2} ; [nE*ones(totperm,1) perm_out{nDrug} nDrug*ones(totperm,1)]];
        end
    end
    save('../../Tables/model_SWdens_est_v6_byE','SWdens_est');
else
    load('../../Tables/model_SWdens_est_v6_byE');
end
%% Filter clusters
clus_alpha=0.05;
montecarlo_alpha=0.05/3;

cfg_neighb=[];
cfg_neighb.method = 'template';
cfg_neighb.layout='biosemi64.lay';
cfg.channel=unique(table_SW.Elec);
neighbours = ft_prepare_neighbours(cfg_neighb);
neighbours(~ismember({neighbours.label},unique(table_SW.Elec)))=[];
[SWdens_clus]=get_clusterperm_lme_lsneurom(SWdens_est,clus_alpha,montecarlo_alpha,totperm,neighbours,1);
%%
cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);
limNumClus=2;
limMax=8;%max(max(abs(temp_topo_tval)));
figure; set(gcf,'Position',[213         173        1027         805/3]);
ClustersByDrugs=cell(2,3);
for nDrug=1:3
    subplot(1,3,nDrug)
    
    temp_topo=SWdens_est{1}(SWdens_est{1}(:,5)==nDrug,3);
    temp_topo2=zeros(size(temp_topo));
    temp_topo3=zeros(size(temp_topo));
    temp_clus=SWdens_clus{nDrug};
    %     temp_topo(temp_pV>= fdr(temp_pV,0.05))=0;
    
    
    if ~isempty(temp_clus)
        for nclus=1:length(temp_clus)
            if length(match_str(layout.label,temp_clus{nclus}{2}))<limNumClus
                continue;
            end
            %             ft_plot_lay_me(layout, 'chanindx',match_str(layout.label,temp_clus{nclus}{2}),'pointsymbol','o','pointcolor','r','pointsize',64,'box','no','label','yes')
            temp_topo2(match_str(layout.label,temp_clus{nclus}{2}))=temp_topo(match_str(layout.label,temp_clus{nclus}{2}));
            temp_topo3(match_str(layout.label,temp_clus{nclus}{2}))=1;
            fprintf('... ... found %s cluster (%g) of %g electrodes (tval cluster=%g, Pmc=%g)\n',temp_clus{nclus}{1},nclus,length(temp_clus{nclus}{2}),temp_clus{nclus}{3},temp_clus{nclus}{4})
            if strcmp(temp_clus{nclus}{1},'pos')
                ClustersByDrugs{2,nDrug}=[ClustersByDrugs{2,nDrug} ; temp_clus{nclus}{2}];
            elseif strcmp(temp_clus{nclus}{1},'neg')
                ClustersByDrugs{1,nDrug}=[ClustersByDrugs{1,nDrug} ; temp_clus{nclus}{2}];
            end
        end
    end
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    ft_plot_lay_me(layout, 'chanindx',1:length(layout.label)-2,'pointsymbol','o','pointcolor',[1 1 1]*0.7,'pointsize',6,'box','no','label','no')
    format_fig;
    caxis([-1 1]*limMax)
    if ~isempty(temp_clus)
        for nclus=1:length(temp_clus)
            if length(match_str(layout.label,temp_clus{nclus}{2}))<limNumClus
                continue;
            end
            ft_plot_lay_me(layout, 'chanindx',match_str(layout.label,temp_clus{nclus}{2}),'pointsymbol','o','pointcolor','k','pointsize',36,'box','no','label','no')
        end
    end
    if nDrug==3
        h=colorbar;
        set(h,'Position',[0.93 0.4 0.02 0.2])
    end
    colormap(cmap2);
    
    title(Drugs{nDrug+1})
end
print('-dpng', '-r300', '../../Figures/Topo_LS_LME_DrugEffect_v6.png')

%% Make new table
% table_behav.C_Pos_ATM=nan(size(table_avSW,1),1);
% table_behav.C_Neg_ATM=nan(size(table_avSW,1),1);
table_behav.C_Pos_CIT=nan(size(table_behav,1),1);
% table_avSW.C_Neg_CIT=nan(size(table_avSW,1),1);
% table_avSW.C_Pos_MPH=nan(size(table_avSW,1),1);
% table_behav.C_Neg_MPH=nan(size(table_behav,1),1);

table_behav.Drug=table_behav.Treatment;
uS=unique(table_behav.SubID);
for nS=1:length(uS)
    for nB=1:10
%         table_avSW.C_Neg_ATM(ismember(table_avSW.SubID,uS(nS)) & table_avSW.BlockN==nB & ismember(table_avSW.Drug,'ATM'))=...
%             nanmean(table_SW.SWdens(ismember(table_SW.SubID,uS(nS)) & table_SW.BlockN==nB & ismember(table_SW.Drug,'ATM') & ismember(table_SW.Elec,ClustersByDrugs{1,1})));
%         table_avSW.C_Neg_ATM(ismember(table_avSW.SubID,uS(nS)) & table_avSW.BlockN==nB & ismember(table_avSW.Drug,'PLA'))=...
%             nanmean(table_SW.SWdens(ismember(table_SW.SubID,uS(nS)) & table_SW.BlockN==nB & ismember(table_SW.Drug,'PLA') & ismember(table_SW.Elec,ClustersByDrugs{1,1})));
%         
%         table_avSW.C_Pos_ATM(ismember(table_avSW.SubID,uS(nS)) & table_avSW.BlockN==nB & ismember(table_avSW.Drug,'ATM'))=...
%             nanmean(table_SW.SWdens(ismember(table_SW.SubID,uS(nS)) & table_SW.BlockN==nB & ismember(table_SW.Drug,'ATM') & ismember(table_SW.Elec,ClustersByDrugs{2,1})));
%         table_avSW.C_Pos_ATM(ismember(table_avSW.SubID,uS(nS)) & table_avSW.BlockN==nB & ismember(table_avSW.Drug,'PLA'))=...
%             nanmean(table_SW.SWdens(ismember(table_SW.SubID,uS(nS)) & table_SW.BlockN==nB & ismember(table_SW.Drug,'PLA') & ismember(table_SW.Elec,ClustersByDrugs{2,1})));
%         
%         table_avSW.C_Neg_CIT(ismember(table_avSW.SubID,uS(nS)) & table_avSW.BlockN==nB & ismember(table_avSW.Drug,'CIT'))=...
%             nanmean(table_SW.SWdens(ismember(table_SW.SubID,uS(nS)) & table_SW.BlockN==nB & ismember(table_SW.Drug,'CIT') & ismember(table_SW.Elec,ClustersByDrugs{1,2})));
%         table_avSW.C_Neg_CIT(ismember(table_avSW.SubID,uS(nS)) & table_avSW.BlockN==nB & ismember(table_avSW.Drug,'PLA'))=...
%             nanmean(table_SW.SWdens(ismember(table_SW.SubID,uS(nS)) & table_SW.BlockN==nB & ismember(table_SW.Drug,'PLA') & ismember(table_SW.Elec,ClustersByDrugs{1,2})));
        
        table_behav.C_Pos_CIT(ismember(table_behav.SubID,uS(nS)) & table_behav.BlockN==nB & ismember(table_behav.Drug,'CIT'))=...
            nanmean(table_SW.SWdens(ismember(table_SW.SubID,uS(nS)) & table_SW.BlockN==nB & ismember(table_SW.Drug,'CIT') & ismember(table_SW.Elec,ClustersByDrugs{2,2})));
        table_behav.C_Pos_CIT(ismember(table_behav.SubID,uS(nS)) & table_behav.BlockN==nB & ismember(table_behav.Drug,'PLA'))=...
            nanmean(table_SW.SWdens(ismember(table_SW.SubID,uS(nS)) & table_SW.BlockN==nB & ismember(table_SW.Drug,'PLA') & ismember(table_SW.Elec,ClustersByDrugs{2,2})));
        
%         table_avSW.C_Neg_MPH(ismember(table_avSW.SubID,uS(nS)) & table_avSW.BlockN==nB & ismember(table_avSW.Drug,'MPH'))=...
%             nanmean(table_SW.SWdens(ismember(table_SW.SubID,uS(nS)) & table_SW.BlockN==nB & ismember(table_SW.Drug,'MPH') & ismember(table_SW.Elec,ClustersByDrugs{1,3})));
%         table_avSW.C_Neg_MPH(ismember(table_avSW.SubID,uS(nS)) & table_avSW.BlockN==nB & ismember(table_avSW.Drug,'PLA'))=...
%             nanmean(table_SW.SWdens(ismember(table_SW.SubID,uS(nS)) & table_SW.BlockN==nB & ismember(table_SW.Drug,'PLA') & ismember(table_SW.Elec,ClustersByDrugs{1,3})));
        
%         table_avSW.C_Pos_MPH(ismember(table_avSW.SubID,uS(nS)) & table_avSW.BlockN==nB & ismember(table_avSW.Drug,'MPH'))=...
%             nanmean(table_SW.SWdens(ismember(table_SW.SubID,uS(nS)) & table_SW.BlockN==nB & ismember(table_SW.Drug,'MPH') & ismember(table_SW.Elec,ClustersByDrugs{2,3})));
%         table_avSW.C_Pos_MPH(ismember(table_avSW.SubID,uS(nS)) & table_avSW.BlockN==nB & ismember(table_avSW.Drug,'PLA'))=...
%             nanmean(table_SW.SWdens(ismember(table_SW.SubID,uS(nS)) & table_SW.BlockN==nB & ismember(table_SW.Drug,'PLA') & ismember(table_SW.Elec,ClustersByDrugs{2,3})));
    end
end
writetable(table_behav,[save_path filesep 'CTET_SWdetection_thr90_byE_P2P_avDens_behav_vec_v6_CLUSTERS.txt']);

% figure;
% for nDrug=1:3
%     subplot(1,3,nDrug)
%     temp_topo=temp_topo_tval(:,nDrug+4);
%     temp_pV=temp_topo_pval(:,nDrug+4);
% %     temp_topo(temp_pV>= fdr(temp_pV,0.05))=0;
%     simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
%     format_fig;
%     caxis([-1 1]*limMax)
%     if nDrug==3
%         h=colorbar;
%         set(h,'Position',[0.93 0.4 0.02 0.2])
%     end
%     colormap(cmap2);
%
%     title(sprintf('BlockN*%s',Drugs{nDrug+1}))
% end