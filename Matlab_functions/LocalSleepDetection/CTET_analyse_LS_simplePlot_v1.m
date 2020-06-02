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

cmap=cbrewer('seq','YlOrRd',64);
Drugs={'PLA','ATM','CIT','MPH'};

limMax=13;

figure;
temp_topo=[];
for nE=1:64
    temp_topo(nE)=mean(table_SW.SWdens(~cellfun(@isempty,regexp(table_SW.Elec,layout.label{nE}))));
end
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
title('All'); h=colorbar;
caxis([0 1]*limMax)
h=colorbar;
colormap(cmap);
ylabel(h, 'waves/min')
%     set(h,'Position',[0.85 0.7 0.04 0.2])
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
table_SW.SubID=categorical(table_SW.SubID);
table_SW.SessN=categorical(table_SW.SessN);
table_SW.Elec=categorical(table_SW.Elec);
table_SW.Drug=categorical(table_SW.Drug);
table_SW.Drug=reordercats(table_SW.Drug,[4 1 2 3]);

mdl0=fitlme(table_SW,'SWdens~1+(1|SubID)');
mdl1=fitlme(table_SW,'SWdens~1+BlockN+(1|SubID)');
mdl2=fitlme(table_SW,'SWdens~1+BlockN+Drug+(1|SubID)');
mdl3=fitlme(table_SW,'SWdens~1+BlockN*Drug+(1|SubID)');
