%%
clear all
close all

run ../localdef.m
addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_LSCPtools));

% table=readtable('/Users/tand0009/Data/CTET_Dockree/CTET_behav_res.txt');
table_SW=readtable([save_path filesep 'CTET_SWdetection_thr90_allE_P2P_vec.txt']);

%%
cfg = [];
cfg.layout = 'biosemi64.lay';
layout=ft_prepare_layout(cfg);

cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
Drugs={'PLA','ATM','CIT','MPH'};

%%
table_SW2=table_SW;
table_SW2.SubID=categorical(table_SW2.SubID);
table_SW2.SessN=categorical(table_SW2.SessN);
table_SW2.Elec=categorical(table_SW2.Elec);
table_SW2.Drug=categorical(table_SW2.Drug);
table_SW2.Drug=reordercats(table_SW2.Drug,[4 1 2 3]);

temp_topo_tval=[];
temp_topo_pval=[];
for nE=1:64
    sub_table_SW2=table_SW2(find_trials(table_SW.Elec,layout.label{nE}),:);
    mdl_byEle{nE}=fitlme(sub_table_SW2,'SWdens~1+Drug+(1|SubID)');
    temp_topo_tval(nE,:)=double(mdl_byEle{nE}.Coefficients(2:4,4));
    temp_topo_pval(nE,:)=double(mdl_byEle{nE}.Coefficients(2:4,6));
end


%%
cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);

limMax=max(max(abs(temp_topo_tval)));
figure; set(gcf,'position',[16         321        1224         657]);
for nDrug=1:3
    subplot(1,3,nDrug)
    simpleTopoPlot_ft(temp_topo_tval(:,nDrug), layout,'on',[],0,1);
    format_fig;
    caxis([-1 1]*limMax)
    if nDrug==3
        h=colorbar;
        set(h,'Position',[0.93 0.4 0.02 0.2])
    end
    colormap(cmap2);
    title(Drugs{nDrug+1})
    
    sigElec=find(temp_topo_pval(:,nDrug)<fdr(temp_topo_pval(:,nDrug),0.05));
    if ~isempty(sigElec)
        ft_plot_lay_me(layout, 'chanindx',sigElec,'pointsymbol','o','pointcolor','r','pointsize',64,'box','no','label','no')
    end
end