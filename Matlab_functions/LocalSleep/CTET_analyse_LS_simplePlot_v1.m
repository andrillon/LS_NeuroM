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

% prctile(table_SW.SWdens,99)
%%
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel=unique(table_SW.Elec);
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
Drugs={'PLA','ATM','CIT','MPH'};

limMax=[3.5 5.5];

figure;
temp_topo=[];
for nE=1:length(layout.label)-2
    temp_topo(nE)=nanmean(table_SW.SWdens(find_trials(table_SW.Elec,layout.label{nE})));
end
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title('All');
caxis(limMax)
h=colorbar;
colormap(cmap);
ylabel(h, 'waves/min')
%     set(h,'Position',[0.85 0.7 0.04 0.2])
format_fig;
print('-dpng', '-r300', '../../Figures/Topo_LS_SWdens.png')

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
print('-dpng', '-r300', '../../Figures/Topo_LS_SWdens_byBlock.png')

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
print('-dpng', '-r300', '../../Figures/Line_LS_SWdens_byBlock.png')

%% By Drug
figure; hold on; set(gcf,'Position',[201         428        1135         500]);
for nD=1:4
    subplot(1,4,nD);
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
print('-dpng', '-r300', '../../Figures/Topo_LS_SWdens_byDrug.png')

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
print('-dpng', '-r300', '../../Figures/Line_LS_SWdens_byDrugAndBlock.png')


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
if redo==1
    totperm=1000;
    temp_topo_tval=[];
    temp_topo_pval=[];
    % fprintf('%2.0f/%2.0f\n',0,64)
    SWdens_est=cell(1,2);
    for nE=1:length(layout.label)-2
        fprintf('%2.0f/%2.0f\n',nE,64)
        sub_table_SW2=table_SW2(find_trials(table_SW.Elec,layout.label{nE}),:);
        mdl_byEle{nE}=fitlme(sub_table_SW2,'SWdens~1+BlockN+Drug+(1|SubID)');
        temp_topo_tval(nE,:)=double(mdl_byEle{nE}.Coefficients(2:end,4));
        temp_topo_pval(nE,:)=double(mdl_byEle{nE}.Coefficients(2:end,6));
        
        [real_out, perm_out]=lme_perm_lsneurom(sub_table_SW2,'Drug','SWdens~1+pred+BlockN+(1|SubID)',totperm);
        SWdens_est{1}=[SWdens_est{1} ; [nE*ones(3,1) real_out]];
        for nDrug=1:3
            SWdens_est{2}=[SWdens_est{2} ; [nE*ones(totperm,1) perm_out{nDrug} nDrug*ones(totperm,1)]];
        end
    end
    save('../../Tables/model_SWdens_est_v3','SWdens_est');
else
    load('../../Tables/model_SWdens_est_v3');
end
%% Filter clusters
clus_alpha=0.01;
montecarlo_alpha=0.05/3;

cfg_neighb=[];
cfg_neighb.method = 'template';
cfg_neighb.layout='biosemi64.lay';
cfg.channel=unique(table_SW.Elec);
neighbours = ft_prepare_neighbours(cfg_neighb);
neighbours(~ismember({neighbours.label},unique(table_SW.Elec)))=[];
[SWdens_clus]=get_clusterperm_lme_lsneurom(SWdens_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
%%
cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);
limNumClus=0;
limMax=10;%max(max(abs(temp_topo_tval)));
figure; set(gcf,'Position',[213         173        1027         805]);
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
            ft_plot_lay_me(layout, 'chanindx',match_str(layout.label,temp_clus{nclus}{2}),'pointsymbol','o','pointcolor','k','pointsize',12,'box','no','label','no')
        end
    end
    if nDrug==3
        h=colorbar;
        set(h,'Position',[0.93 0.4 0.02 0.2])
    end
    colormap(cmap2);
    
    title(Drugs{nDrug+1})
end
print('-dpng', '-r300', '../../Figures/Topo_LS_LME_DrugEffect_v3.png')
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