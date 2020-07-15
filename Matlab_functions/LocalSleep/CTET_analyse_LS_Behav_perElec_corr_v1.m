%%
clear all
close all

run ../localdef.m
addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_LSCPtools));
addpath(genpath([pwd filesep '..']));

% table=readtable('/Users/tand0009/Data/CTET_Dockree/CTET_behav_res.txt');
table_SW=readtable([save_path filesep 'CTET_SWdetection_thr90_byE_P2P_behav_vec_full_v2.txt']);

%%
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel=unique(table_SW.Elec);
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

uniqueS=unique(table_SW2.SubID);
for nS=1:length(uniqueS)
    for nBl=1:10
        CIT_SWdens=table_SW2.SWdens(table_SW2.BlockN==nBl & table_SW2.SubID==uniqueS(nS) & table_SW2.Drug=='CIT');
    end
end

%%
Contrasts={{'MPH','PLA'},{'CIT','PLA'},{'ATM','PLA'}};
contrast_table=table_SW2(ismember(table_SW2.Drug,{'CIT','MPH','ATM'}),:);
contrast_table.Drug=removecats(contrast_table.Drug);
for nC=1:3
    % fprintf('%2.0f/%2.0f\n',0,64)
    for nE=1:length(layout.label)-2
        fprintf('%2.0f/%2.0f\n',nE,length(layout.label)-2)
        sub_table1=contrast_table(ismember(contrast_table.Elec,layout.label{nE}) & ismember(contrast_table.Drug,Contrasts{nC}{1}),:);
        sub_table2=table_SW2(ismember(table_SW2.Elec,layout.label{nE}) & ismember(table_SW2.Drug,Contrasts{nC}{2}),:);

        sub_table1.P2P=sub_table1.P2P-sub_table2.P2P;
        sub_table1.Miss=sub_table1.Miss-sub_table2.Miss;
        sub_table1.FA=sub_table1.FA-sub_table2.FA;
        sub_table1.Hit_RT=sub_table1.Hit_RT-sub_table2.Hit_RT;

        contrast_table(ismember(contrast_table.Elec,layout.label{nE}) & ismember(contrast_table.Drug,Contrasts{nC}{1}),:)=sub_table1;

                sub_table1=table_SW2(ismember(table_SW2.Elec,layout.label{nE}) & ismember(table_SW2.Drug,Contrasts{nC}{1}),:);
        sub_table2=table_SW2(ismember(table_SW2.Elec,layout.label{nE}) & ismember(table_SW2.Drug,Contrasts{nC}{2}),:);
        sub_table1.diffSW=sub_table1.P2P-sub_table2.P2P;
        sub_table1.diffMiss=sub_table1.Miss-sub_table2.Miss;
        sub_table1.diffFA=sub_table1.FA-sub_table2.FA;
        sub_table1.diffRT=sub_table1.Hit_RT-sub_table2.Hit_RT;
        
        mdl=fitlme(sub_table1,'diffFA~1+diffSW+(1|SubID)');
        FA_est{nC}(nE,:)=double(mdl.Coefficients(2,[2 4 6]));
        mdl=fitlme(sub_table1,'diffMiss~1+diffSW+(1|SubID)');
        Miss_est{nC}(nE,:)=double(mdl.Coefficients(2,[2 4 6]));
        mdl=fitlme(sub_table1,'diffRT~1+diffSW+(1|SubID)');
        Hit_RT_est{nC}(nE,:)=double(mdl.Coefficients(2,[2 4 6]));
    end
end

for nE=1:length(layout.label)-2
    fprintf('%2.0f/%2.0f\n',nE,length(layout.label)-2)
    sub_table1=contrast_table(ismember(contrast_table.Elec,layout.label{nE}),:);
    mdl=fitlme(sub_table1,'FA~1+P2P+(1|SubID)');
    FA_est{4}(nE,:)=double(mdl.Coefficients(2,[2 4 6]));
    mdl=fitlme(sub_table1,'Miss~1+P2P+(1|SubID)');
    Miss_est{4}(nE,:)=double(mdl.Coefficients(2,[2 4 6]));
    mdl=fitlme(sub_table1,'Hit_RT~1+P2P+(1|SubID)');
    Hit_RT_est{4}(nE,:)=double(mdl.Coefficients(2,[2 4 6]));
end

%%
%%
cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);
limNumClus=0;
limMax=4;%max(max(abs(temp_topo_tval)));
figure; set(gcf,'Position',[213         173        1027         805/3]);
PlotTitles={'FA','Miss','Hit_RT'};
    for nPlot=1:3
        
        subplot(1,3,nPlot)
        switch nPlot
            case 1
                temp_topo=FA_est{4}(:,3);
            case 2
                temp_topo=Miss_est{4}(:,3);
            case 3
                temp_topo=Hit_RT_est{4}(:,3);
        end
        
        simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
        ft_plot_lay_me(layout, 'chanindx',1:64,'pointsymbol','o','pointcolor',[1 1 1]*0.7,'pointsize',6,'box','no','label','no')
        format_fig;
        caxis([-1 1]*limMax)
        
        if nPlot==3
            h=colorbar;
            set(h,'Position',[0.93 0.4 0.02 0.2])
        end
        colormap(cmap2);
        title({PlotTitles{nPlot},sprintf('%s vs %s',Contrasts{nC}{1},Contrasts{nC}{2})})
    end

%%
cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);
limNumClus=2;
limMax=4;%max(max(abs(temp_topo_tval)));
figure; set(gcf,'Position',[213         173        1027         805]);
PlotTitles={'FA','Miss','Hit_RT'};
for nC=1:3
    for nPlot=1:3
        
        subplot(3,3,(nC-1)*3+nPlot)
        switch nPlot
            case 1
                temp_topo=FA_est{nC}(:,3);
            case 2
                temp_topo=Miss_est{nC}(:,3);
            case 3
                temp_topo=Hit_RT_est{nC}(:,3);
        end
        
        simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
        ft_plot_lay_me(layout, 'chanindx',1:64,'pointsymbol','o','pointcolor',[1 1 1]*0.7,'pointsize',6,'box','no','label','no')
        format_fig;
        caxis([-1 1]*limMax)
        
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