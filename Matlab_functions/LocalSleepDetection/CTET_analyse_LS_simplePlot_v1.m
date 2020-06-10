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

limMax=13;

figure;
temp_topo=[];
for nE=1:64
    temp_topo(nE)=mean(table_SW.SWdens(find_trials(table_SW.Elec,layout.label{nE})));
end
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title('All'); 
caxis([0 1]*limMax)
h=colorbar;
colormap(cmap);
ylabel(h, 'waves/min')
%     set(h,'Position',[0.85 0.7 0.04 0.2])
format_fig;

%% By Block
limMax=13;

figure;
for nBl=1:10
    subplot(2,5,nBl)
    sub_table_SW=table_SW(table_SW.BlockN==nBl,:);
    temp_topo=[];
    for nE=1:64
        temp_topo(nE)=mean(sub_table_SW.SWdens(find_trials(sub_table_SW.Elec,layout.label{nE})));
    end
    simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
    title(sprintf('Block %g',nBl));
    caxis([0 1]*limMax)
%     h=colorbar;
    colormap(cmap);
%     ylabel(h, 'waves/min')
    %     set(h,'Position',[0.85 0.7 0.04 0.2])
    format_fig;
end

figure;
temp_plot=[];
for nBl=1:10
    sub_table_SW=table_SW(table_SW.BlockN==nBl,:);
    temp_plot(1,nBl)=mean(sub_table_SW.SWdens);
    temp_plot(2,nBl)=sem(sub_table_SW.SWdens);
end
plot(1:10,temp_plot(1,:),'Color','k','Marker','o','LineWidth',3)
hold on;
errorbar(1:10,temp_plot(1,:),temp_plot(2,:),'Color','k','LineWidth',2)
    format_fig;

%% By Drug
figure; hold on; set(gcf,'Position',[201         428        1135         500]);
for nD=1:4
    subplot(1,4,nD);
    sub_table_SW=table_SW(ismember(table_SW.Drug,Drugs{nD}),:);
    temp_topo=[];
    for nE=1:64
        temp_topo(nE)=mean(sub_table_SW.SWdens(~cellfun(@isempty,regexp(sub_table_SW.Elec,layout.label{nE}))));
    end
    simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
    format_fig;
    caxis([0 1]*limMax)
    if nD==4
        h=colorbar;
        set(h,'Position',[0.93 0.4 0.02 0.2])
    end
    colormap(cmap);
    title(Drugs{nD})
end

% %% By Drug and block
% figure; hold on;
% for nD=1:4
%     for nBl=1:10
%         subplot(4,10,10*(nD-1)+nBl);
%         sub_table_SW=table_SW(ismember(table_SW.Drug,Drugs{nD}) & table_SW.BlockN==nBl,:);
%         temp_topo=[];
%         for nE=1:64
%             temp_topo(nE)=mean(sub_table_SW.SWdens(~cellfun(@isempty,regexp(sub_table_SW.Elec,layout.label{nE}))));
%         end
%         simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
%         format_fig;
%         caxis([0 1]*limMax)
%         h=colorbar;
%         colormap(cmap);
%     end
% end


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
figure;
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
end