%%
clear all
close all

run ../localdef.m
addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_LSCPtools));
addpath(genpath([pwd filesep '..']));
addpath(genpath(path_mediation));
addpath(genpath(path_CanlabCore));

% table=readtable('/Users/tand0009/Data/CTET_Dockree/CTET_behav_res.txt');
table_SW=readtable([save_path filesep 'CTET_SWdetection_thr90_byE_P2P_behav_vec_full_v2.txt']);

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
nc=0;
Subs=unique(table_SW2.SubID);
Labels=unique(table_SW2.Elec);
sub_table_SW=table_SW2((table_SW2.Drug=='CIT' | table_SW2.Drug=='MPH'),:);
sub_table_SW.Drug=removecats(sub_table_SW.Drug);
all_stats=cell(1,4);
for nE=1:length(Labels)
    clear X Y M
    nSc=0;
    for nS=1:length(Subs)
        %     if length(unique(sub_table_SW.Drug(sub_table_SW.SubID==Subs(nS))))==2
        %         nc=nc+1;
        %         X{nc}={sub_table_SW.Miss(sub_table_SW.Elec=='Fz' & sub_table_SW.SubID==Subs(nS))};
        %         Y{nc}={sub_table_SW.Drug(sub_table_SW.Elec=='Fz' & sub_table_SW.SubID==Subs(nS))=='MPH'};
        %         M{nc}={sub_table_SW.SWdens(sub_table_SW.Elec=='Fz' & sub_table_SW.SubID==Subs(nS))};
        if max(sub_table_SW.SWdens(sub_table_SW.Elec==Labels(nE)  & (sub_table_SW.SubID)==Subs(nS)))==0
            continue;
        end
        nSc=nSc+1;
        %         for nBl=1:10
        X{1,nSc}=sub_table_SW.Miss(sub_table_SW.Elec=='Fz'  & (sub_table_SW.SubID)==Subs(nS));
        Y{1,nSc}=double(sub_table_SW.Drug(sub_table_SW.Elec=='Fz'  & (sub_table_SW.SubID)==Subs(nS))=='CIT');
        M{1,nSc}=sub_table_SW.SWdens(sub_table_SW.Elec==Labels(nE)  & (sub_table_SW.SubID)==Subs(nS));
        %         end
        %     end
    end
    [paths, stats2] = mediation(X, Y, M, 'verbose', 'bootsamples', 1000);
    all_stats{1}=[all_stats{1} ; [stats2.mean nE]];
    all_stats{2}=[all_stats{2} ; [stats2.t nE]];
    all_stats{3}=[all_stats{3} ; [stats2.p nE]];
    all_stats{4}=[all_stats{4} ; [paths nE*ones(size(paths,1),1) (1:size(paths,1))']];
end

%%
figure;
% subplot(1,3,1)
% simpleTopoPlot_ft(all_stats{2}(:,4), layout,'on',[],0,1);
subplot(1,2,1)
temp=all_stats{2}(:,3);
tempPv=all_stats{3}(:,3);
temp(tempPv>0.05)=1;
simpleTopoPlot_ft(temp, layout,'on',[],0,1);
caxis([-1 1]*3)

subplot(1,2,2)
temp=all_stats{2}(:,5);
tempPv=all_stats{3}(:,5);
temp(tempPv>0.05)=1;
simpleTopoPlot_ft(temp, layout,'on',[],0,1);
caxis([-1 1]*3)

%%
figure;
for nleg=1:5
    subplot(1,5,nleg)
    temp=grpstats(all_stats{4}(:,nleg),all_stats{4}(:,6));
%     tempPv=all_stats{3}(:,3);
%     temp(tempPv>0.05)=1;
    simpleTopoPlot_ft(temp, layout,'on',[],0,1);
    caxis([-1 1]*max(abs(temp)))
    colorbar;
end