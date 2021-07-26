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

table_SW2.dprime=nan(size(table_SW2,1),1);
table_SW2.crit=nan(size(table_SW2,1),1);
for k=1:size(table_SW2,1)
    [dp,crit]=calc_dprime(table_SW2.Hit(k),table_SW2.FA(k));
    table_SW2.dprime(k)=dp;
    table_SW2.crit(k)=crit;
end
%%
redo=0;
totperm=1000;
if redo==1
    % fprintf('%2.0f/%2.0f\n',0,64)
    Miss_est=cell(1,2);
    FA_est=cell(1,2);
    Hit_RT_est=cell(1,2);
    dprime_est=cell(1,2);
    crit_est=cell(1,2);
    for nE=1:size(layout.label,1)-2
        fprintf('%2.0f/%2.0f\n',nE,size(layout.label,1)-2)
        sub_table_SW2=table_SW2(match_str(table_SW.Elec,layout.label{nE}),:);
        
        %%%% FA
        if nE==1
            [real_out, perm_out, out_perm_FA]=lme_perm_lsneurom(sub_table_SW2,'SWdens','FA~1+pred+(1|SubID)',totperm);
        else
            [real_out, perm_out]=lme_perm_lsneurom(sub_table_SW2,'SWdens','FA~1+pred+(1|SubID)',totperm,out_perm_FA);
        end
        FA_est{1}=[FA_est{1} ; [nE real_out]];
        FA_est{2}=[FA_est{2} ; [nE*ones(totperm,1) perm_out]];
        
        %%%% MISS
        if nE==1
            [real_out, perm_out, out_perm_Miss]=lme_perm_lsneurom(sub_table_SW2,'SWdens','Miss~1+pred+(1|SubID)',totperm);
        else
            [real_out, perm_out]=lme_perm_lsneurom(sub_table_SW2,'SWdens','Miss~1+pred+(1|SubID)',totperm,out_perm_Miss);
        end
        Miss_est{1}=[Miss_est{1} ; [nE real_out]];
        Miss_est{2}=[Miss_est{2} ; [nE*ones(totperm,1) perm_out]];
        
        %%%% RT
        if nE==1
            [real_out, perm_out, out_perm_RT]=lme_perm_lsneurom(sub_table_SW2,'SWdens','RT~1+pred+(1|SubID)',totperm);
        else
            [real_out, perm_out]=lme_perm_lsneurom(sub_table_SW2,'SWdens','RT~1+pred+(1|SubID)',totperm,out_perm_RT);
        end
        Hit_RT_est{1}=[Hit_RT_est{1} ; [nE real_out]];
        Hit_RT_est{2}=[Hit_RT_est{2} ; [nE*ones(totperm,1) perm_out]];
        
        %%%% dprime
        if nE==1
            [real_out, perm_out, out_perm_RT]=lme_perm_lsneurom(sub_table_SW2,'SWdens','dprime~1+pred+(1|SubID)',totperm);
        else
            [real_out, perm_out]=lme_perm_lsneurom(sub_table_SW2,'SWdens','dprime~1+pred+(1|SubID)',totperm,out_perm_RT);
        end
        dprime_est{1}=[dprime_est{1} ; [nE real_out]];
        dprime_est{2}=[dprime_est{2} ; [nE*ones(totperm,1) perm_out]];
        
        
        %%%% crit
        if nE==1
            [real_out, perm_out, out_perm_RT]=lme_perm_lsneurom(sub_table_SW2,'SWdens','crit~1+pred+(1|SubID)',totperm);
        else
            [real_out, perm_out]=lme_perm_lsneurom(sub_table_SW2,'SWdens','crit~1+pred+(1|SubID)',totperm,out_perm_RT);
        end
        crit_est{1}=[crit_est{1} ; [nE real_out]];
        crit_est{2}=[crit_est{2} ; [nE*ones(totperm,1) perm_out]];
    end
    save('../../Tables/model_Behav_SW_est_v6','Miss_est','Hit_RT_est','FA_est','dprime_est','crit_est');
else
    load('../../Tables/model_Behav_SW_est_v6');
end

if redo==1
    % fprintf('%2.0f/%2.0f\n',0,64)
    Miss_est=cell(1,2);
    FA_est=cell(1,2);
    Hit_RT_est=cell(1,2);
    dprime_est=cell(1,2);
    crit_est=cell(1,2);
    for nE=1:size(layout.label,1)-2
        fprintf('%2.0f/%2.0f\n',nE,size(layout.label,1)-2)
        sub_table_SW2=table_SW2(match_str(table_SW.Elec,layout.label{nE}),:);
        
        %%%% FA
        if nE==1
            [real_out, perm_out, out_perm_FA]=lme_perm_lsneurom(sub_table_SW2,'Alpha','FA~1+pred+(1|SubID)',totperm);
        else
            [real_out, perm_out]=lme_perm_lsneurom(sub_table_SW2,'Alpha','FA~1+pred+(1|SubID)',totperm,out_perm_FA);
        end
        FA_est{1}=[FA_est{1} ; [nE real_out]];
        FA_est{2}=[FA_est{2} ; [nE*ones(totperm,1) perm_out]];
        
        %%%% MISS
        if nE==1
            [real_out, perm_out, out_perm_Miss]=lme_perm_lsneurom(sub_table_SW2,'Alpha','Miss~1+pred+(1|SubID)',totperm);
        else
            [real_out, perm_out]=lme_perm_lsneurom(sub_table_SW2,'Alpha','Miss~1+pred+(1|SubID)',totperm,out_perm_Miss);
        end
        Miss_est{1}=[Miss_est{1} ; [nE real_out]];
        Miss_est{2}=[Miss_est{2} ; [nE*ones(totperm,1) perm_out]];
        
        %%%% RT
        if nE==1
            [real_out, perm_out, out_perm_RT]=lme_perm_lsneurom(sub_table_SW2,'Alpha','RT~1+pred+(1|SubID)',totperm);
        else
            [real_out, perm_out]=lme_perm_lsneurom(sub_table_SW2,'Alpha','RT~1+pred+(1|SubID)',totperm,out_perm_RT);
        end
        Hit_RT_est{1}=[Hit_RT_est{1} ; [nE real_out]];
        Hit_RT_est{2}=[Hit_RT_est{2} ; [nE*ones(totperm,1) perm_out]];
        
        %%%% dprime
        if nE==1
            [real_out, perm_out, out_perm_RT]=lme_perm_lsneurom(sub_table_SW2,'Alpha','dprime~1+pred+(1|SubID)',totperm);
        else
            [real_out, perm_out]=lme_perm_lsneurom(sub_table_SW2,'Alpha','dprime~1+pred+(1|SubID)',totperm,out_perm_RT);
        end
        dprime_est{1}=[dprime_est{1} ; [nE real_out]];
        dprime_est{2}=[dprime_est{2} ; [nE*ones(totperm,1) perm_out]];
        
        
        %%%% crit
        if nE==1
            [real_out, perm_out, out_perm_RT]=lme_perm_lsneurom(sub_table_SW2,'Alpha','crit~1+pred+(1|SubID)',totperm);
        else
            [real_out, perm_out]=lme_perm_lsneurom(sub_table_SW2,'Alpha','crit~1+pred+(1|SubID)',totperm,out_perm_RT);
        end
        crit_est{1}=[crit_est{1} ; [nE real_out]];
        crit_est{2}=[crit_est{2} ; [nE*ones(totperm,1) perm_out]];
    end
    save('../../Tables/model_Behav_Alpha_est_v6','Miss_est','Hit_RT_est','FA_est','dprime_est','crit_est');
else
    load('../../Tables/model_Behav_Alpha_est_v6');
end
%% Filter clusters - SW in addition to Drug
load('../../Tables/model_Behav_SW_est_v6');
clus_alpha=0.05;
montecarlo_alpha=0.05;

cfg_neighb=[];
cfg_neighb.method = 'template';
cfg_neighb.layout='biosemi64.lay';
cfg_neighb.channel=layout.label;
neighbours = ft_prepare_neighbours(cfg_neighb);

[FA_clus]=get_clusterperm_lme_lsneurom(FA_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
[Miss_clus]=get_clusterperm_lme_lsneurom(Miss_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
[Hit_RT_clus]=get_clusterperm_lme_lsneurom(Hit_RT_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
[dprime_clus]=get_clusterperm_lme_lsneurom(dprime_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
[crit_clus]=get_clusterperm_lme_lsneurom(crit_est,clus_alpha,montecarlo_alpha,totperm,neighbours);



cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);
limNumClus=2;
limMax=6;%max(max(abs(temp_topo_tval)));
figure; set(gcf,'Position',[1 1 320*3 320*1]);
[ha pos]=tight_subplot(1,3,0.02,0.05,0.05);
PlotTitles={'Miss','FA','RT'};
for nPlot=1:3
    
    hs=subplot(1,3,nPlot); format_fig;
    set(hs,'Position',pos{nPlot})
    switch nPlot
        case 2
            temp_topo=FA_est{1}(:,3);
            temp_clus=FA_clus;
        case 1
            temp_topo=Miss_est{1}(:,3);
            temp_clus=Miss_clus;
        case 3
            temp_topo=Hit_RT_est{1}(:,3);
            temp_clus=Hit_RT_clus;
        case 4
            temp_topo=dprime_est{1}(:,3);
            temp_clus=dprime_clus;
        case 5
            temp_topo=crit_est{1}(:,3);
            temp_clus=crit_clus;
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
%     ft_plot_lay_me(layout, 'chanindx',1:64,'pointsymbol','o','pointcolor',[1 1 1]*0.7,'pointsize',6,'box','no','label','no')
    format_fig;
    caxis([-1 1]*limMax)
    if ~isempty(temp_clus)
        for nclus=1:length(temp_clus)
            if length(match_str(layout.label,temp_clus{nclus}{2}))<limNumClus
                continue;
            end
            ft_plot_lay_me(layout, 'chanindx',match_str(layout.label,temp_clus{nclus}{2}),'pointsymbol','o','pointcolor','k','pointsize',72,'box','no','label','no')
        end
    end
    if nPlot==3
        hb=colorbar('Position',[0.975    0.33   0.02    0.33]);

    end
    colormap(cmap2);
    
    title(PlotTitles{nPlot})
end
print('-dpng', '-r300', '../../Figures/Topo_LME_BehavEffect_byBlock_onTopOfDrug.png')


%% Filter clusters - SW in addition to Drug
load('../../Tables/model_Behav_Alpha_est_v6');
clus_alpha=0.05;
montecarlo_alpha=0.05;

cfg_neighb=[];
cfg_neighb.method = 'template';
cfg_neighb.layout='biosemi64.lay';
cfg_neighb.channel=layout.label;
neighbours = ft_prepare_neighbours(cfg_neighb);

[FA_clus]=get_clusterperm_lme_lsneurom(FA_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
[Miss_clus]=get_clusterperm_lme_lsneurom(Miss_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
[Hit_RT_clus]=get_clusterperm_lme_lsneurom(Hit_RT_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
[dprime_clus]=get_clusterperm_lme_lsneurom(dprime_est,clus_alpha,montecarlo_alpha,totperm,neighbours);
[crit_clus]=get_clusterperm_lme_lsneurom(crit_est,clus_alpha,montecarlo_alpha,totperm,neighbours);



cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);
limNumClus=2;
limMax=6;%max(max(abs(temp_topo_tval)));
figure; set(gcf,'Position',[1 1 320*3 320*1]);
[ha pos]=tight_subplot(1,3,0.02,0.05,0.05);
PlotTitles={'Miss','FA','RT'};
for nPlot=1:3
    
    hs=subplot(1,3,nPlot); format_fig;
    set(hs,'Position',pos{nPlot})
    switch nPlot
        case 2
            temp_topo=FA_est{1}(:,3);
            temp_clus=FA_clus;
        case 1
            temp_topo=Miss_est{1}(:,3);
            temp_clus=Miss_clus;
        case 3
            temp_topo=Hit_RT_est{1}(:,3);
            temp_clus=Hit_RT_clus;
        case 4
            temp_topo=dprime_est{1}(:,3);
            temp_clus=dprime_clus;
        case 5
            temp_topo=crit_est{1}(:,3);
            temp_clus=crit_clus;
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
%     ft_plot_lay_me(layout, 'chanindx',1:64,'pointsymbol','o','pointcolor',[1 1 1]*0.7,'pointsize',6,'box','no','label','no')
    format_fig;
    caxis([-1 1]*limMax)
    if ~isempty(temp_clus)
        for nclus=1:length(temp_clus)
            if length(match_str(layout.label,temp_clus{nclus}{2}))<limNumClus
                continue;
            end
            ft_plot_lay_me(layout, 'chanindx',match_str(layout.label,temp_clus{nclus}{2}),'pointsymbol','o','pointcolor','k','pointsize',72,'box','no','label','no')
        end
    end
    if nPlot==3
        hb=colorbar('Position',[0.975    0.33   0.02    0.33]);

    end
    colormap(cmap2);
    
    title(PlotTitles{nPlot})
end
print('-dpng', '-r300', '../../Figures/Topo_LME_BehavEffect_Alpha_byBlock_onTopOfDrug.png')

