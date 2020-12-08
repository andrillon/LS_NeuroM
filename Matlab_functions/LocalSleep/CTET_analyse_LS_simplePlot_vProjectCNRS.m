%%
clear all
close all

run ../localdef.m
addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_LSCPtools));
addpath(genpath([pwd filesep '..']));

% table=readtable('/Users/tand0009/Data/CTET_Dockree/CTET_behav_res.txt');
table_SW=readtable([save_path filesep 'CTET_SWdetection_thr90_byE_P2P_behav_vec_full_v3.txt']);
table_avSW=readtable([save_path filesep 'CTET_SWdetection_thr90_byE_P2P_avDens_behav_vec_full_v3.txt']);

% prctile(table_SW.SWdens,99)
%%
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel=unique(table_SW.Elec);
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)


%% By Drug
Drugs={'CIT','PLA','MPH'};
limMax=[3.5 6.5];
figure; hold on; set(gcf,'Position',[201         428        1135         500]);
tight_subplot(1,3);
for nD=1:length(Drugs)
    subplot(1,length(Drugs),nD);
    sub_table_SW=table_SW(ismember(table_SW.Drug,Drugs{nD}),:);
    temp_topo=[];
    for nE=1:length(layout.label)-2
        temp_topo(nE)=nanmean(sub_table_SW.SWdens(~cellfun(@isempty,regexp(sub_table_SW.Elec,layout.label{nE}))));
    end
    simpleTopoPlot_ft(temp_topo', layout,'off',[],0,1);
    format_fig;
    caxis(limMax)
    if nD==3
        h=colorbar;
        set(h,'Position',[0.93 0.4 0.02 0.2])
    end
    colormap(cmap);
    title(Drugs{nD})
end
print('-dpng', '-r300', '../../Figures/Topo_LS_SWdens_byDrug_CNRS.png')

%%
figure;
Drugs={'PLA','ATM','CIT','MPH'};
for nD=1:length(Drugs)
    if nD==2
        continue;
    end
    temp_block=[];
    for nBl=1:10
        sub_table_SW=table_SW(ismember(table_SW.Drug,Drugs{nD}) & table_SW.BlockN==nBl,:);
        temp=grpstats(sub_table_SW.SWdens,sub_table_SW.SubID);
        temp_block(1,nBl)=nanmean(temp);
        temp_block(2,nBl)=sem(temp);
    end
    plot((1:10)+(nD-2.5)/10,temp_block(1,:),'Color',Colors(nD,:),'LineWidth',3);
    hold on;
    errorbar((1:10)+(nD-2.5)/10,temp_block(1,:),temp_block(2,:),'Color',Colors(nD,:),'LineWidth',2)
    scatter((1:10)+(nD-2.5)/10,temp_block(1,:),'SizeData',144,'MarkerEdgeColor',Colors(nD,:),'MarkerFaceColor',Colors(nD,:),'MarkerFaceAlpha',0.7)
end
format_fig;
xlabel('Block')
xlim([0.5 10.5])
ylabel('waves/min')
% print('-dpng', '-r300', '../../Figures/Line_LS_SWdens_byDrugAndBlock_CNRS.png')

%%
Drugs={'PLA','ATM','CIT','MPH'};
figure;
set(gcf,'Position',[680   674   307   304]);
nc=0;
for nD=1:4
     if nD==2
        continue;
     end
    nc=nc+1;
    temp_plot=[];
    sub_table_SW=table_SW(ismember(table_SW.Drug,Drugs{nD}),:);
    temp_plot=grpstats(sub_table_SW.SWdens,sub_table_SW.SubID);
    %     simpleDotPlot(nD,temp_plot,144,Colors(nD,:),1,'k','o',[],3,1,0);
    %     h1 = raincloud_plot(temp_plot, 'box_on', 1, 'color', Colors(nD,:), 'alpha', 0.5,...
    %         'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
    %         'box_col_match', 1);
    simpleBarPlot(nc,temp_plot,ColorsD{match_str(ColorsDlabels,Drugs{nD})},0.85,'k',[],5);
end
xlim([0.2 3.8])
ylabel('waves/min')
set(gca,'XTick',1:3,'XTickLabel',{'PLA','CIT','MPH'});
format_fig;
% print('-dpng', '-r300', '../../Figures/Bar_LS_SWdens_byDrugAndBlock_CNRS.png')
export_fig(['../../Figures/Bar_LS_SWdens_byDrugAndBlock_CNRS.eps'],'-r 300')

%%
figure;
set(gcf,'Position',[680   674   307   304]);
nc=0;
for nD=1:4
     if nD==2
        continue;
     end
    nc=nc+1;
    temp_plot=[];
    sub_table_SW=table_SW(ismember(table_SW.Drug,Drugs{nD}),:);
    temp_plot=grpstats(sub_table_SW.Miss,sub_table_SW.SubID);
    %     simpleDotPlot(nD,temp_plot,144,Colors(nD,:),1,'k','o',[],3,1,0);
    %     h1 = raincloud_plot(temp_plot, 'box_on', 1, 'color', Colors(nD,:), 'alpha', 0.5,...
    %         'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
    %         'box_col_match', 1);
    simpleBarPlot(nc,100*temp_plot,[1 1 1;ColorsD{match_str(ColorsDlabels,Drugs{nD})}],0.85,'k',[],5);
end
xlim([0.2 3.8])
ylabel('misses (%)')
set(gca,'XTick',1:3,'XTickLabel',{'PLA','CIT','MPH'});

format_fig;
% print('-dpng', '-r300', '../../Figures/Bar_LS_Misses_byDrugAndBlock_CNRS.png')
export_fig(['../../Figures/Bar_LS_Misses_byDrugAndBlock_CNRS.eps'],'-r 300')
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
load('../../Tables/model_SWdens_est_v3');

% Filter clusters
clus_alpha=0.05;
montecarlo_alpha=0.05/3;
totperm=1000;
cfg_neighb=[];
cfg_neighb.method = 'template';
cfg_neighb.layout='biosemi64.lay';
cfg.channel=unique(table_SW.Elec);
neighbours = ft_prepare_neighbours(cfg_neighb);
neighbours(~ismember({neighbours.label},unique(table_SW.Elec)))=[];
[SWdens_clus]=get_clusterperm_lme_lsneurom(SWdens_est,clus_alpha,montecarlo_alpha,totperm,neighbours,1);
%
cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);
limNumClus=0;
limMax=10;%max(max(abs(temp_topo_tval)));
figure; set(gcf,'Position',[213         173        1027         805]);
ClustersByDrugs=cell(2,3);
for nDrug=2:3
    subplot(1,2,nDrug-1)
    
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
    simpleTopoPlot_ft(temp_topo2, layout,'on',[],0,1);
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
export_fig(['../../Figures/Topo_LS_LME_CNRS.eps'],'-r 300')

