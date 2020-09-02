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
table_SW=readtable([save_path filesep 'CTET_SWdetection_thr90_byE_P2P_behav_vec_full_v3.txt']);

%%
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel=unique(table_SW.Elec);
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

%%
nc=0;
Subs=unique(table_SW2.SubID);
Labels=unique(table_SW2.Elec);
rmpath('/Users/tand0009/Work/local/fieldtrip/external/stats/');
    all_stats=cell(3,3);
    all_Pvs=[];
for nD=1:3
    all_stats{nD,1}=cell(1,4);
    all_stats{nD,2}=cell(1,4);
    all_stats{nD,3}=cell(1,4);
    
    sub_table_SW=table_SW2((table_SW2.Drug=='PLA' | table_SW2.Drug==Drugs{nD+1}),:);
    sub_table_SW.Drug=removecats(sub_table_SW.Drug);
    for nE=1:length(layout.label)-2
        clear Y1 Y2 Y3 X M
        nSc=0;
        for nS=1:length(Subs)
            %     if length(unique(sub_table_SW.Drug(sub_table_SW.SubID==Subs(nS))))==2
            %         nc=nc+1;
            %         X{nc}={sub_table_SW.Miss(sub_table_SW.Elec=='Fz' & sub_table_SW.SubID==Subs(nS))};
            %         Y{nc}={sub_table_SW.Drug(sub_table_SW.Elec=='Fz' & sub_table_SW.SubID==Subs(nS))=='MPH'};
            %         M{nc}={sub_table_SW.SWdens(sub_table_SW.Elec=='Fz' & sub_table_SW.SubID==Subs(nS))};
            if max(sub_table_SW.SWdens(sub_table_SW.Elec==layout.label{nE}  & (sub_table_SW.SubID)==Subs(nS)))==0
                continue;
            end
            nSc=nSc+1;
            %         for nBl=1:10
            Y1{1,nSc}=sub_table_SW.Hit_RT(sub_table_SW.Elec==layout.label{nE}  & (sub_table_SW.SubID)==Subs(nS));
            Y2{1,nSc}=sub_table_SW.Miss(sub_table_SW.Elec==layout.label{nE}  & (sub_table_SW.SubID)==Subs(nS));
            Y3{1,nSc}=sub_table_SW.FA(sub_table_SW.Elec==layout.label{nE}  & (sub_table_SW.SubID)==Subs(nS));
            X{1,nSc}=double(sub_table_SW.Drug(sub_table_SW.Elec==layout.label{nE}  & (sub_table_SW.SubID)==Subs(nS))==Drugs{nD+1});
            M{1,nSc}=sub_table_SW.SWdens(sub_table_SW.Elec==layout.label{nE}  & (sub_table_SW.SubID)==Subs(nS));
            %         end
            %     end
        end
        [paths, stats2] = mediation(X, Y1, M, 'verbose', 'bootstrapfirst','bootsamples', 10000);
        all_stats{nD,1}{1}=[all_stats{nD,1}{1} ; [stats2.mean nE]];
        all_stats{nD,1}{2}=[all_stats{nD,1}{2} ; [stats2.t nE]];
        all_stats{nD,1}{3}=[all_stats{nD,1}{3} ; [stats2.p nE]];
        all_stats{nD,1}{4}=[all_stats{nD,1}{4} ; [paths nE*ones(size(paths,1),1) (1:size(paths,1))']];
        all_Pvs=[all_Pvs stats2.p];
        
        [paths, stats2] = mediation(X, Y2, M, 'verbose', 'bootstrapfirst','bootsamples', 10000);
        all_stats{nD,2}{1}=[all_stats{nD,2}{1} ; [stats2.mean nE]];
        all_stats{nD,2}{2}=[all_stats{nD,2}{2} ; [stats2.t nE]];
        all_stats{nD,2}{3}=[all_stats{nD,2}{3} ; [stats2.p nE]];
        all_stats{nD,2}{4}=[all_stats{nD,2}{4} ; [paths nE*ones(size(paths,1),1) (1:size(paths,1))']];
                all_Pvs=[all_Pvs stats2.p];

        [paths, stats2] = mediation(X, Y3, M, 'verbose', 'bootstrapfirst','bootsamples', 10000);
        all_stats{nD,3}{1}=[all_stats{nD,3}{1} ; [stats2.mean nE]];
        all_stats{nD,3}{2}=[all_stats{nD,3}{2} ; [stats2.t nE]];
        all_stats{nD,3}{3}=[all_stats{nD,3}{3} ; [stats2.p nE]];
        all_stats{nD,3}{4}=[all_stats{nD,3}{4} ; [paths nE*ones(size(paths,1),1) (1:size(paths,1))']];
         all_Pvs=[all_Pvs stats2.p];
   end
end
%%
for nD=1:3
    figure; set(gcf,'Name',Drugs{nD+1})
    % subplot(1,3,1)
    % simpleTopoPlot_ft(all_stats{2}(:,4), layout,'on',[],0,1);
    subplot(2,3,1)
    temp=all_stats{nD,1}{2}(:,1);
    tempPv=all_stats{nD,1}{3}(:,1);
%     temp(tempPv>fdr(all_Pvs,0.05))=0;
    simpleTopoPlot_ft(temp, layout,'on',[],0,1);
    caxis([-1 1]*5)
    title('RT - c''')
    
    subplot(2,3,4)
    temp=all_stats{nD,1}{2}(:,2);
    tempPv=all_stats{nD,1}{3}(:,2);
%     temp(tempPv>fdr(all_Pvs,0.05))=0;
    simpleTopoPlot_ft(temp, layout,'on',[],0,1);
    caxis([-1 1]*5)
    title('RT - ab')
    
    subplot(2,3,2)
    temp=all_stats{nD,2}{2}(:,1);
    tempPv=all_stats{nD,2}{3}(:,1);
%     temp(tempPv>fdr(all_Pvs,0.05))=0;
    simpleTopoPlot_ft(temp, layout,'on',[],0,1);
    caxis([-1 1]*5)
    title('MISS - c''')
    
    subplot(2,3,5)
    temp=all_stats{nD,2}{2}(:,2);
    tempPv=all_stats{nD,2}{3}(:,2);
%     temp(tempPv>fdr(all_Pvs,0.05))=0;
    simpleTopoPlot_ft(temp, layout,'on',[],0,1);
    caxis([-1 1]*5)
    title('MISS - ab')
    
    subplot(2,3,3)
    temp=all_stats{nD,3}{2}(:,1);
    tempPv=all_stats{nD,3}{3}(:,1);
%     temp(tempPv>fdr(all_Pvs,0.05))=0;
    simpleTopoPlot_ft(temp, layout,'on',[],0,1);
    caxis([-1 1]*5)
    title('FA - c''')
    
    subplot(2,3,6)
    temp=all_stats{nD,3}{2}(:,2);
    tempPv=all_stats{nD,3}{3}(:,2);
%     temp(tempPv>fdr(all_Pvs,0.05))=0;
    simpleTopoPlot_ft(temp, layout,'on',[],0,1);
    caxis([-1 1]*5)
    title('FA - ab')
end
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