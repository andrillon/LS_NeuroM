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
table_SW2=readtable([save_path filesep 'CTET_SWflag_perTrial_byElec_v6.txt']);

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

table_SW2.SubID=categorical(table_SW2.SubID);
table_SW2.SessN=categorical(table_SW2.SessN);
table_SW2.Drug=categorical(table_SW2.Treatment);
table_SW2.Drug=reordercats(table_SW2.Drug,[4 1 2 3]);

%%
redo=1;
totperm=1;
SWdens_est=cell(3,2);
for nBehav=1:3
    SWdens_est{nBehav,1}=nan(length(layout.label)-2,4);
    SWdens_est{nBehav,2}=nan(totperm*(length(layout.label)-2),5);
    if nBehav==1
        Behav_Var='RT';
        sub_table_SW2=table_SW2(ismember(table_SW2.Drug,{'CIT','PLA'}) & ~isnan(table_SW2.RT) & table_SW2.StimType==1,:);
    elseif nBehav==2
        Behav_Var='Miss';
        sub_table_SW2=table_SW2(ismember(table_SW2.Drug,{'CIT','PLA'}) & ~isnan(table_SW2.corrTG) & table_SW2.StimType==1,:);
        sub_table_SW2.Miss=1-sub_table_SW2.corrTG;
    elseif nBehav==3
        Behav_Var='FA';
        sub_table_SW2=table_SW2(ismember(table_SW2.Drug,{'CIT','PLA'}) & ~isnan(table_SW2.corrNT) & table_SW2.StimType==0,:);
        sub_table_SW2.FA=1-sub_table_SW2.corrNT;
    end
    
    mdl_byEle=[];
    % fprintf('%2.0f/%2.0f\n',0,64)
    for nE=1:length(layout.label)-2
        fprintf('B%g: %2.0f/%2.0f\n',nBehav,nE,64)
        if strcmp(Behav_Var,'Miss') || strcmp(Behav_Var,'FA')
            mdl_byEle{nE}=fitglme(sub_table_SW2,sprintf('%s~1+BlockN+%s+(1|SubID)',Behav_Var,layout.label{nE}),'Distribution','binomial');
        else
            mdl_byEle{nE}=fitlme(sub_table_SW2,sprintf('%s~1+BlockN+%s+(1|SubID)',Behav_Var,layout.label{nE}));
        end
%         temp_topo_tval(nE,:)=double(mdl_byEle{nE}.Coefficients(2:end,4));
%         temp_topo_pval(nE,:)=double(mdl_byEle{nE}.Coefficients(2:end,6));
                CITonly_noD_topo_tval(nBehav,nE,:)=double(mdl_byEle{nE}.Coefficients(2:end,4));
        CITonly_noD_topo_pval(nBehav,nE,:)=double(mdl_byEle{nE}.Coefficients(2:end,6));
        if strcmp(Behav_Var,'Miss') || strcmp(Behav_Var,'FA')
            [real_out, perm_out]=glme_perm_lsneurom(sub_table_SW2,layout.label{nE},sprintf('%s~1+BlockN+pred+(1|SubID)',Behav_Var),totperm);
        else
            [real_out, perm_out]=lme_perm_lsneurom(sub_table_SW2,layout.label{nE},sprintf('%s~1+BlockN+pred+(1|SubID)',Behav_Var),totperm);
        end
        SWdens_est{nBehav,1}(nE,:)=[nE real_out];
        SWdens_est{nBehav,2}(((nE-1)*totperm+1):nE*totperm,:)=[nE*ones(totperm,1) perm_out];
    end
end
%%
SWdens_est2=cell(3,2);
CITonly_topo_tval=[];
CITonly_topo_pval=[];
for nBehav=1:3
    SWdens_est2{nBehav,1}=nan(length(layout.label)-2,4);
    SWdens_est2{nBehav,2}=nan(totperm*(length(layout.label)-2),5);
    if nBehav==1
        Behav_Var='RT';
        sub_table_SW2=table_SW2(ismember(table_SW2.Drug,{'CIT','PLA'}) & ~isnan(table_SW2.RT) & table_SW2.StimType==1,:);
    elseif nBehav==2
        Behav_Var='Miss';
        sub_table_SW2=table_SW2(ismember(table_SW2.Drug,{'CIT','PLA'}) & ~isnan(table_SW2.corrTG) & table_SW2.StimType==1,:);
        sub_table_SW2.Miss=1-sub_table_SW2.corrTG;
    elseif nBehav==3
        Behav_Var='FA';
        sub_table_SW2=table_SW2(ismember(table_SW2.Drug,{'CIT','PLA'}) & ~isnan(table_SW2.corrNT) & table_SW2.StimType==0,:);
        sub_table_SW2.FA=1-sub_table_SW2.corrNT;
    end
    sub_table_SW2.Drug=removecats(sub_table_SW2.Drug);
    

    mdl_byEle=[];
    % fprintf('%2.0f/%2.0f\n',0,64)
    for nE=1:length(layout.label)-2
        fprintf('B%g: %2.0f/%2.0f\n',nBehav,nE,64)
        if strcmp(Behav_Var,'Miss') || strcmp(Behav_Var,'FA')
            mdl_byEle{nE}=fitglme(sub_table_SW2,sprintf('%s~1+BlockN+Drug+%s+(1|SubID)',Behav_Var,layout.label{nE}),'Distribution','binomial');
        else
            mdl_byEle{nE}=fitlme(sub_table_SW2,sprintf('%s~1+BlockN+Drug+%s+(1|SubID)',Behav_Var,layout.label{nE}));
        end
        CITonly_topo_tval(nBehav,nE,:)=double(mdl_byEle{nE}.Coefficients(2:end,4));
        CITonly_topo_pval(nBehav,nE,:)=double(mdl_byEle{nE}.Coefficients(2:end,6));
        
        if strcmp(Behav_Var,'Miss') || strcmp(Behav_Var,'FA')
            [real_out, perm_out]=glme_perm_lsneurom(sub_table_SW2,layout.label{nE},sprintf('%s~1+BlockN+Drug+pred+(1|SubID)',Behav_Var),totperm);
        else
            [real_out, perm_out]=lme_perm_lsneurom(sub_table_SW2,layout.label{nE},sprintf('%s~1+BlockN+Drug+pred+(1|SubID)',Behav_Var),totperm);
        end
        SWdens_est2{nBehav,1}(nE,:)=[nE real_out];
        SWdens_est2{nBehav,2}(((nE-1)*totperm+1):nE*totperm,:)=[nE*ones(totperm,1) perm_out];
    end
end

%% CHECK
CITonly_topo_tval2=[];
CITonly_topo_pval2=[];
for nBehav=1:3
    SWdens_est2{nBehav,1}=nan(length(layout.label)-2,4);
    SWdens_est2{nBehav,2}=nan(totperm*(length(layout.label)-2),5);
    if nBehav==1
        Behav_Var='RT';
        sub_table_SW2=table_SW2(ismember(table_SW2.Drug,{'CIT','PLA'}) & ~isnan(table_SW2.RT) & table_SW2.StimType==1,:);
    elseif nBehav==2
        Behav_Var='Miss';
        sub_table_SW2=table_SW2(ismember(table_SW2.Drug,{'CIT','PLA'}) & ~isnan(table_SW2.corrTG) & table_SW2.StimType==1,:);
        sub_table_SW2.Miss=1-sub_table_SW2.corrTG;
    elseif nBehav==3
        Behav_Var='FA';
        sub_table_SW2=table_SW2(ismember(table_SW2.Drug,{'CIT','PLA'}) & ~isnan(table_SW2.corrNT) & table_SW2.StimType==0,:);
        sub_table_SW2.FA=1-sub_table_SW2.corrNT;
    end
    sub_table_SW2.Drug=removecats(sub_table_SW2.Drug);
    

    mdl_byEle=[];
    % fprintf('%2.0f/%2.0f\n',0,64)
    for nE=1:length(layout.label)-2
        fprintf('B%g: %2.0f/%2.0f\n',nBehav,nE,64)
            mdl_byEle{nE}=fitglme(sub_table_SW2,sprintf('%s~1+BlockN+Drug+(1|SubID)',layout.label{nE}),'Distribution','binomial');
        CITonly_topo_tval2(nE,:)=double(mdl_byEle{nE}.Coefficients(2:end,4));
        CITonly_topo_pval2(nE,:)=double(mdl_byEle{nE}.Coefficients(2:end,6));

    end
end


%% Filter clusters
clus_alpha=0.05;
montecarlo_alpha=0.05;

cfg_neighb=[];
cfg_neighb.method = 'template';
cfg_neighb.layout='biosemi64.lay';
cfg.channel=unique(table_SW.Elec);
neighbours = ft_prepare_neighbours(cfg_neighb);
neighbours(~ismember({neighbours.label},unique(table_SW.Elec)))=[];
[SWdens_clus_RT]=get_clusterperm_lme_lsneurom(SWdens_est(1,:),clus_alpha,montecarlo_alpha,totperm,neighbours,1);
[SWdens_clus_Miss]=get_clusterperm_lme_lsneurom(SWdens_est(2,:),clus_alpha,montecarlo_alpha,totperm,neighbours,1);
[SWdens_clus_FA]=get_clusterperm_lme_lsneurom(SWdens_est(3,:),clus_alpha,montecarlo_alpha,totperm,neighbours,1);
%%
cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);
limNumClus=0;
limMax=10;%max(max(abs(temp_topo_tval)));
figure; set(gcf,'Position',[213         173        1027         805]);
ClustersByDrugs=cell(2,3);
Behav_Vars={'RT','Miss','FA'};
for nBehav=1:3
    subplot(1,3,nBehav)
    if nBehav==1
        SWdens_clus=SWdens_clus_RT;
    elseif nBehav==2
        SWdens_clus=SWdens_clus_Miss;
    elseif nBehav==3
        SWdens_clus=SWdens_clus_FA;
    end
    temp_topo=SWdens_est{nBehav,1}(:,3);
    temp_topo2=zeros(size(temp_topo));
    temp_topo3=zeros(size(temp_topo));
    temp_clus=SWdens_clus;
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
            %             if strcmp(temp_clus{nclus}{1},'pos')
            %                 ClustersByDrugs{2,nDrug}=[ClustersByDrugs{2,nDrug} ; temp_clus{nclus}{2}];
            %             elseif strcmp(temp_clus{nclus}{1},'neg')
            %                 ClustersByDrugs{1,nDrug}=[ClustersByDrugs{1,nDrug} ; temp_clus{nclus}{2}];
            %             end
        end
    end
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    
    if nBehav==4
        h=colorbar;
        set(h,'Position',[0.93 0.4 0.02 0.2])
    end
    colormap(cmap2);
    
    title(Behav_Vars{nBehav})
end
