%%
clear all
% close all

run ../localdef.m
addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_LSCPtools));
addpath(genpath([pwd filesep '..']));
path_mediation=[pwd filesep '..' filesep '..' filesep 'Tables'];

table_SW=readtable([save_path filesep 'CTET_SWdetection_thr90_allE_P2P_behav_vec_full_v3.txt']);

cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel=unique(table_SW.Elec);
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)

%% Load all drugs
filename='CTET_Mediation_allE_ATM_RT_Est_perTrial_v3.txt';
mediation_ATM_RT_Est = CTET_import_mediationoutput([path_mediation filesep filename]); mediation_ATM_RT_Est(~ismember(mediation_ATM_RT_Est.Channels,layout.label),:)=[];
filename='CTET_Mediation_allE_ATM_RT_pV_perTrial_v2.txt';
mediation_ATM_RT_pV = CTET_import_mediationoutput([path_mediation filesep filename]); mediation_ATM_RT_pV(~ismember(mediation_ATM_RT_pV.Channels,layout.label),:)=[];

filename='CTET_Mediation_allE_ATM_Miss_Est_perTrial_v3.txt';
mediation_ATM_Miss_Est = CTET_import_mediationoutput([path_mediation filesep filename]); mediation_ATM_Miss_Est(~ismember(mediation_ATM_Miss_Est.Channels,layout.label),:)=[];
filename='CTET_Mediation_allE_ATM_Miss_pV_perTrial_v3.txt';
mediation_ATM_Miss_pV = CTET_import_mediationoutput([path_mediation filesep filename]); mediation_ATM_Miss_pV(~ismember(mediation_ATM_Miss_pV.Channels,layout.label),:)=[];

myChannels=cellstr(mediation_ATM_Miss_Est.Channels);

%%%
filename='CTET_Mediation_allE_MPH_RT_Est_perTrial_v3.txt';
mediation_MPH_RT_Est = CTET_import_mediationoutput([path_mediation filesep filename]); mediation_MPH_RT_Est(~ismember(mediation_MPH_RT_Est.Channels,layout.label),:)=[];
filename='CTET_Mediation_allE_MPH_RT_pV_perTrial_v3.txt';
mediation_MPH_RT_pV = CTET_import_mediationoutput([path_mediation filesep filename]); mediation_MPH_RT_pV(~ismember(mediation_MPH_RT_pV.Channels,layout.label),:)=[];

filename='CTET_Mediation_allE_MPH_Miss_Est_perTrial_v3.txt';
mediation_MPH_Miss_Est = CTET_import_mediationoutput([path_mediation filesep filename]); mediation_MPH_Miss_Est(~ismember(mediation_MPH_Miss_Est.Channels,layout.label),:)=[];
filename='CTET_Mediation_allE_MPH_Miss_pV_perTrial_v3.txt';
mediation_MPH_Miss_pV = CTET_import_mediationoutput([path_mediation filesep filename]); mediation_MPH_Miss_pV(~ismember(mediation_MPH_Miss_pV.Channels,layout.label),:)=[];

%%%
filename='CTET_Mediation_allE_CIT_RT_Est_perTrial_v3.txt';
mediation_CIT_RT_Est = CTET_import_mediationoutput([path_mediation filesep filename]); mediation_CIT_RT_Est(~ismember(mediation_CIT_RT_Est.Channels,layout.label),:)=[];
filename='CTET_Mediation_allE_CIT_RT_pV_perTrial_v3.txt';
mediation_CIT_RT_pV = CTET_import_mediationoutput([path_mediation filesep filename]); mediation_CIT_RT_pV(~ismember(mediation_CIT_RT_pV.Channels,layout.label),:)=[];

filename='CTET_Mediation_allE_CIT_Miss_Est_perTrial_v3.txt';
mediation_CIT_Miss_Est = CTET_import_mediationoutput([path_mediation filesep filename]); mediation_CIT_Miss_Est(~ismember(mediation_CIT_Miss_Est.Channels,layout.label),:)=[];
filename='CTET_Mediation_allE_CIT_Miss_pV_perTrial_v3.txt';
mediation_CIT_Miss_pV = CTET_import_mediationoutput([path_mediation filesep filename]); mediation_CIT_Miss_pV(~ismember(mediation_CIT_Miss_pV.Channels,layout.label),:)=[];


all_pVs=[mediation_ATM_RT_pV.ACME_treated ; mediation_ATM_Miss_pV.ACME_treated ; 
    mediation_MPH_RT_pV.ACME_treated ; mediation_MPH_Miss_pV.ACME_treated ;
    mediation_CIT_RT_pV.ACME_treated ; mediation_CIT_Miss_pV.ACME_treated ; 
    ];

FDR_Thr=fdr(all_pVs,0.05);
fprintf('... FDR threshold: %g\n',FDR_Thr)
%% ATM
figure;
subplot(1,2,1);
temp_topo=[];
temp_pV=[];
for nE=1:length(layout.label)-2
    temp_topo(nE)=(mediation_ATM_RT_Est.ACME_treated(find(ismember(myChannels,layout.label{nE}))));
    temp_pV(nE)=(mediation_ATM_RT_pV.ACME_treated(find(ismember(myChannels,layout.label{nE}))));
end
% temp_topo(temp_pV>=fdr(temp_pV,0.05))=0;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
if ~isempty(find(temp_pV<FDR_Thr)) %fdr(temp_pV,0.05)))
ft_plot_lay_me(layout, 'chanindx',find(temp_pV<FDR_Thr),'pointsymbol','o','pointcolor','r','pointsize',36,'box','no','label','no')
end
title('ATM RT');
h=colorbar;
colormap('parula');
format_fig;
caxis([-1 1]*max(abs(temp_topo)))

subplot(1,2,2);
temp_topo=[];
temp_pV=[];
for nE=1:length(layout.label)-2
    temp_topo(nE)=(mediation_ATM_Miss_Est.ACME_treated(find(ismember(myChannels,layout.label{nE}))));
    temp_pV(nE)=(mediation_ATM_Miss_pV.ACME_treated(find(ismember(myChannels,layout.label{nE}))));
end
% temp_topo(temp_pV>=FDR_Thr)=0;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
if ~isempty(find(temp_pV<FDR_Thr))
ft_plot_lay_me(layout, 'chanindx',find(temp_pV<FDR_Thr),'pointsymbol','o','pointcolor','r','pointsize',36,'box','no','label','no')
end
title('ATM Miss');
h=colorbar;
colormap('parula');
format_fig;
caxis([-1 1]*max(abs(temp_topo)))

%% MPH
figure;
subplot(1,2,1);
temp_topo=[];
temp_pV=[];
for nE=1:length(layout.label)-2
    temp_topo(nE)=(mediation_MPH_RT_Est.ACME_treated(find(ismember(myChannels,layout.label{nE}))));
    temp_pV(nE)=(mediation_MPH_RT_pV.ACME_treated(find(ismember(myChannels,layout.label{nE}))));
end
% temp_topo(temp_pV>=FDR_Thr)=0;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
if ~isempty(find(temp_pV<FDR_Thr))
ft_plot_lay_me(layout, 'chanindx',find(temp_pV<FDR_Thr),'pointsymbol','o','pointcolor','r','pointsize',36,'box','no','label','no')
end
title('MPH RT');
h=colorbar;
colormap('parula');
format_fig;
caxis([-1 1]*max(abs(temp_topo)))

subplot(1,2,2);
temp_topo=[];
temp_pV=[];
for nE=1:length(layout.label)-2
    temp_topo(nE)=(mediation_MPH_Miss_Est.ACME_treated(find(ismember(myChannels,layout.label{nE}))));
    temp_pV(nE)=(mediation_MPH_Miss_pV.ACME_treated(find(ismember(myChannels,layout.label{nE}))));
end
% temp_topo(temp_pV>=FDR_Thr)=0;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
if ~isempty(find(temp_pV<FDR_Thr))
ft_plot_lay_me(layout, 'chanindx',find(temp_pV<FDR_Thr),'pointsymbol','o','pointcolor','r','pointsize',36,'box','no','label','no')
end
title('MPH Miss');
h=colorbar;
colormap('parula');
format_fig;
caxis([-1 1]*max(abs(temp_topo)))

%% CIT
figure;
subplot(1,2,1);
temp_topo=[];
temp_pV=[];
for nE=1:length(layout.label)-2
    temp_topo(nE)=(mediation_CIT_RT_Est.ACME_treated(find(ismember(myChannels,layout.label{nE}))));
    temp_pV(nE)=(mediation_CIT_RT_pV.ACME_treated(find(ismember(myChannels,layout.label{nE}))));
end
% temp_topo(temp_pV>=FDR_Thr)=0;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
if ~isempty(find(temp_pV<FDR_Thr))
ft_plot_lay_me(layout, 'chanindx',find(temp_pV<FDR_Thr),'pointsymbol','o','pointcolor','r','pointsize',36,'box','no','label','no')
end
title('CIT RT');
h=colorbar;
colormap('parula');
format_fig;
caxis([-1 1]*max(abs(temp_topo)))

subplot(1,2,2);
temp_topo=[];
temp_pV=[];
for nE=1:length(layout.label)-2
    temp_topo(nE)=(mediation_CIT_Miss_Est.ACME_treated(find(ismember(myChannels,layout.label{nE}))));
    temp_pV(nE)=(mediation_CIT_Miss_pV.ACME_treated(find(ismember(myChannels,layout.label{nE}))));
end
% temp_topo(temp_pV>=FDR_Thr)=0;
simpleTopoPlot_ft(temp_topo', layout,'on',[],0,1);
if ~isempty(find(temp_pV<FDR_Thr))
ft_plot_lay_me(layout, 'chanindx',find(temp_pV<FDR_Thr),'pointsymbol','o','pointcolor','r','pointsize',36,'box','no','label','no')
end
title('CIT Miss');
h=colorbar;
colormap('parula');
format_fig;
caxis([-1 1]*max(abs(temp_topo)))
