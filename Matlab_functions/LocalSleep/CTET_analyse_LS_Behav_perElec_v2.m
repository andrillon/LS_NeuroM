%%
clear all
close all

run ../localdef.m
addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_LSCPtools));
addpath(genpath([pwd filesep '..']));

% table=readtable('/Users/tand0009/Data/CTET_Dockree/CTET_behav_res.txt');
table_SW=readtable([save_path filesep 'CTET_SWdetection_thr90_allE_P2P_behav_vec.txt']);

%%
cfg = [];
cfg.layout = 'biosemi64.lay';
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

%%
redo=1;
totperm=1000;
Contrasts={{'MPH','PLA'},{'CIT','PLA'},{'ATM','PLA'}};
if redo==1
    for nC=1:3
        % fprintf('%2.0f/%2.0f\n',0,64)
        Miss_est{nC}=cell(1,2);
        FA_est{nC}=cell(1,2);
        Hit_RT_est{nC}=cell(1,2);
        for nE=1:64
            fprintf('%2.0f/%2.0f\n',nE,64)
            sub_table_SW2=table_SW2(ismember(table_SW.Elec,layout.label{nE}) & ismember(table_SW.Drug,Contrasts{nC}),:);
            
            %%%% FA
            [real_out, perm_out]=lme_perm_lsneurom(sub_table_SW2,'SWdens','FA~1+pred+BlockN+(1|SubID)',totperm);
            FA_est{nC}{1}=[FA_est{nC}{1} ; [nE real_out]];
            FA_est{nC}{2}=[FA_est{nC}{2} ; [nE*ones(totperm,1) perm_out]];
            
            %%%% MISS
            [real_out, perm_out]=lme_perm_lsneurom(sub_table_SW2,'SWdens','Miss~1+pred+BlockN+(1|SubID)',totperm);
            Miss_est{nC}{1}=[Miss_est{nC}{1} ; [nE real_out]];
            Miss_est{nC}{2}=[Miss_est{nC}{2} ; [nE*ones(totperm,1) perm_out]];
            
            %%%% RT
            [real_out, perm_out]=lme_perm_lsneurom(sub_table_SW2,'SWdens','Hit_RT~1+pred+BlockN+(1|SubID)',totperm);
            Hit_RT_est{nC}{1}=[Hit_RT_est{nC}{1} ; [nE real_out]];
            Hit_RT_est{nC}{2}=[Hit_RT_est{nC}{2} ; [nE*ones(totperm,1) perm_out]];
        end
    end
    save('../../Tables/model_Behav_perDrug_est','Miss_est','Hit_RT_est','FA_est');
    
else
    load('../../Tables/model_Behav_perDrug_est');
end
%% Filter clusters
clus_alpha=0.05;
montecarlo_alpha=0.05/9;

cfg_neighb=[];
cfg_neighb.method = 'template';
cfg_neighb.layout='biosemi64.lay';
neighbours = ft_prepare_neighbours(cfg_neighb);

for nC=1:3
[FA_clus{nC}]=get_clusterperm_lme_lsneurom(FA_est{nC},clus_alpha,montecarlo_alpha,totperm,neighbours);
[Miss_clus{nC}]=get_clusterperm_lme_lsneurom(Miss_est{nC},clus_alpha,montecarlo_alpha,totperm,neighbours);
[Hit_RT_clus{nC}]=get_clusterperm_lme_lsneurom(Hit_RT_est{nC},clus_alpha,montecarlo_alpha,totperm,neighbours);
end

%%
cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);
limNumClus=2;
limMax=10;%max(max(abs(temp_topo_tval)));
figure; set(gcf,'Position',[213         173        1027         805]);
PlotTitles={'FA','Miss','Hit_RT'};
for nC=1:3
for nPlot=1:3
    
    subplot(3,3,(nC-1)*3+nPlot)
    switch nPlot
        case 1
            temp_topo=FA_est{nC}{1}(:,3);
            temp_clus=FA_clus{nC};
        case 2
            temp_topo=Miss_est{nC}{1}(:,3);
            temp_clus=Miss_clus{nC};
        case 3
            temp_topo=Hit_RT_est{nC}{1}(:,3);
            temp_clus=Hit_RT_clus{nC};
    end
    temp_topo2=zeros(size(temp_topo));
    temp_topo3=zeros(size(temp_topo));
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
    ft_plot_lay_me(layout, 'chanindx',1:64,'pointsymbol','o','pointcolor',[1 1 1]*0.7,'pointsize',6,'box','no','label','no')
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
    if nPlot==3
        h=colorbar;
        set(h,'Position',[0.93 0.4 0.02 0.2])
    end
    colormap(cmap2);
    
    title({PlotTitles{nPlot},sprintf('%s vs %s',Contrasts{nC}{1},Contrasts{nC}{2})})
end
end
% print('-dpng', '-r300', '../../Figures/Topo_LME_BehavEffect.png')
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