%%
clear all
close all

run ../localdef.m
addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_LSCPtools));
addpath(genpath([pwd filesep '..']));

% table=readtable('/Users/tand0009/Data/CTET_Dockree/CTET_behav_res.txt');
table_SW=readtable([save_path filesep 'CTET_SWdetection_thr90_byE_P2P_behav_TF_vec_v6.txt']);

%%
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel=unique(table_SW.Elec);
cfg.channel(match_str(cfg.channel,'Iz'))=[];
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
Drugs={'PLA','ATM','CIT','MPH'};

%% Models
table_SW2=table_SW;
table_SW2.SubID=categorical(table_SW2.SubID);
table_SW2.SessN=categorical(table_SW2.SessN);
table_SW2.Elec=categorical(table_SW2.Elec);
table_SW2.Drug=categorical(table_SW2.Drug);
table_SW2.Drug=reordercats(table_SW2.Drug,[4 1 2 3]);

table_SW2.Miss=1-table_SW2.Hit;
%%
redo=1;
totperm=100;
if redo==1
    % fprintf('%2.0f/%2.0f\n',0,64)
    SSVEP_est=cell(1,2);
    Alpha_est=cell(1,2);
    for nE=1:size(layout.label,1)-2
        fprintf('%2.0f/%2.0f\n',nE,size(layout.label,1)-2)
        sub_table_SW2=table_SW2(match_str(table_SW.Elec,layout.label{nE}),:);
        
        %%%% Alpha
        if nE==1
            [real_out, perm_out, out_perm_FA]=lme_perm_lsneurom(sub_table_SW2,'SWdens','Alpha~1+pred+Drug+(1|SubID)',totperm);
        else
            [real_out, perm_out]=lme_perm_lsneurom(sub_table_SW2,'SWdens','Alpha~1+pred+Drug+(1|SubID)',totperm,out_perm_FA);
        end
        Alpha_est{1}=[Alpha_est{1} ; [nE real_out]];
        Alpha_est{2}=[Alpha_est{2} ; [nE*ones(totperm,1) perm_out]];
        
        %%%% MISS
        if nE==1
            [real_out, perm_out, out_perm_Miss]=lme_perm_lsneurom(sub_table_SW2,'SWdens','SSVEP~1+pred+Drug+(1|SubID)',totperm);
        else
            [real_out, perm_out]=lme_perm_lsneurom(sub_table_SW2,'SWdens','SSVEP~1+pred+Drug+(1|SubID)',totperm,out_perm_Miss);
        end
        SSVEP_est{1}=[SSVEP_est{1} ; [nE real_out]];
        SSVEP_est{2}=[SSVEP_est{2} ; [nE*ones(totperm,1) perm_out]];
        
    end
    save('../../Tables/model_TF_Drug_est_v6b','Alpha_est','SSVEP_est');
else
    load('../../Tables/model_TF_Drug_est_v6b');
end
%% Filter clusters - SW in addition to Drug
clus_alpha=0.05;
montecarlo_alpha=0.05;

cfg_neighb=[];
cfg_neighb.method = 'template';
cfg_neighb.layout='biosemi64.lay';
cfg_neighb.channel=layout.label;
neighbours = ft_prepare_neighbours(cfg_neighb);

[Alpha_clus]=get_clusterperm_lme_lsneurom(Alpha_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
[SSVEP_clus]=get_clusterperm_lme_lsneurom(SSVEP_est,clus_alpha,montecarlo_alpha,totperm,neighbours);


cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);
limNumClus=2;
limMax=12;%max(max(abs(temp_topo_tval)));
figure; set(gcf,'Position',[213         173        1027         805/3]);
PlotTitles={'Alpha','SSVEP'};
for nPlot=1:2
    
    subplot(1,2,nPlot)
    switch nPlot
        case 1
            temp_topo=Alpha_est{1}(:,3);
            temp_clus=Alpha_clus;
        case 2
            temp_topo=SSVEP_est{1}(:,3);
            temp_clus=SSVEP_clus;
    end
    temp_topo2=zeros(size(temp_topo));
    temp_topo3=zeros(size(temp_topo));
    %         temp_topo(temp_pV>= fdr(temp_pV,0.05))=0;
    
    
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
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    ft_plot_lay_me(layout, 'chanindx',1:64,'pointsymbol','o','pointcolor',[1 1 1]*0.7,'pointsize',6,'box','no','label','no')
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
    if nPlot==2
        h=colorbar;
        set(h,'Position',[0.93 0.4 0.02 0.2])
    end
    colormap(cmap2);
    
    title(PlotTitles{nPlot})
end
print('-dpng', '-r300', '../../Figures/Topo_LME_TFEffect_byBlock_onTopOfDrug.png')


%% Filter clusters - SW instead of Drug
load('../../Tables/model_Behav_est_v6b');
clus_alpha=0.05;
montecarlo_alpha=0.05;

cfg_neighb=[];
cfg_neighb.method = 'template';
cfg_neighb.layout='biosemi64.lay';
cfg_neighb.channel=layout.label;
neighbours = ft_prepare_neighbours(cfg_neighb);

[Alpha_clus]=get_clusterperm_lme_lsneurom(Alpha_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
[SSVEP_clus]=get_clusterperm_lme_lsneurom(SSVEP_est,clus_alpha,montecarlo_alpha,totperm,neighbours);


cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);
limNumClus=2;
limMax=6;%max(max(abs(temp_topo_tval)));
figure; set(gcf,'Position',[213         173        1027         805/3]);
PlotTitles={'Alpha','SSVEP'};
for nPlot=1:2
    
    subplot(1,2,nPlot)
    switch nPlot
        case 1
            temp_topo=Alpha_est{1}(:,3);
            temp_clus=Alpha_clus;
        case 2
            temp_topo=SSVEP_est{1}(:,3);
            temp_clus=SSVEP_clus;
    end
    temp_topo2=zeros(size(temp_topo));
    temp_topo3=zeros(size(temp_topo));
    %         temp_topo(temp_pV>= fdr(temp_pV,0.05))=0;
    
    
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
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    ft_plot_lay_me(layout, 'chanindx',1:64,'pointsymbol','o','pointcolor',[1 1 1]*0.7,'pointsize',6,'box','no','label','no')
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
    if nPlot==2
        h=colorbar;
        set(h,'Position',[0.93 0.4 0.02 0.2])
    end
    colormap(cmap2);
    
    title(PlotTitles{nPlot})
end
print('-dpng', '-r300', '../../Figures/Topo_LME_TFEffect_byBlock_insteadOfDrug.png')

