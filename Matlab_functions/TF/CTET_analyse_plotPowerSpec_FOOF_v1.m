%%
clear all
close all
run ../localdef.m
addpath(genpath([pwd filesep '..']))

path_fieldtrip='/Users/tand0009/Work/local/fieldtrip/';
path_localsleep='/Users/tand0009/WorkGit/projects/inprogress/wanderIM/localsleep';
addpath(path_fieldtrip);
addpath(path_localsleep);
ft_defaults;
fooof_path='/Users/tand0009/WorkGit/projects/ext/fooof_mat/';
addpath(genpath(fooof_path))


path_LSCPtools='/Users/tand0009/WorkGit/LSCPtools/';
addpath(genpath(path_LSCPtools));

data_path='/Volumes/tLab_BackUp1/Monash/CTET_Dockree/EEG_CTET/';
save_path='/Users/tand0009/Data/CTET_Dockree/';
% data_path='/Volumes/shared/R-MNHS-SPP/Bellgrove-data/Jess Barnes EEG Backup Data/EEG_CTET/';

table=readtable('/Users/tand0009/Data/CTET_Dockree/CTET_behav_res.txt');

load(['/Users/tand0009/Data/CTET_Dockree/headers_BDFfiles']);
%%
files=dir(['/Users/tand0009/Data/CTET_Dockree/CTET_*_FOOOF_FFT_perBlock_byElec_avMast.mat']);

nFc=0;
av_bg=[];
av_peaks=[];
all_Drugs=[];
for nF=1:length(files)
    File_Name=files(nF).name;
    fprintf('... processing %s\n',File_Name);
    septag=findstr(File_Name,'_');
    SubN=str2num(File_Name(septag(1)+1:septag(2)-1));
    SessN=str2num(File_Name(septag(4)-1));
    DrugC=(File_Name(septag(4)+1:septag(4)+3));
    
     load(['/Users/tand0009/Data/CTET_Dockree/' File_Name]);
        nFc=nFc+1;
        av_pow(nFc,:,:,:)=all_pow;
        all_Drugs{nF}=DrugC;
        
        av_bg=[av_bg ; all_bg];
        av_peaks=[av_peaks ; all_peaks];
end

%%
thiCh=match_str(hdr.label,'Pz');
figure;
    hold on;
plot(myfreqs,mean(squeeze(mean(av_pow(:,:,thiCh,:),2)))); %,0,Colors(nDrug,:),0,'-',.5,1,0,[],2);
 xlim([2 30])

%%
thiCh=match_str(hdr.label,'Pz');
figure;
    hold on;
Drugs={'PLA','MPH','CIT','ATM'};
for nDrug=1:4
[~,hp(nDrug)]=simpleTplot(myfreqs,squeeze(mean(av_pow(ismember(all_Drugs,Drugs{nDrug}),:,thiCh,:),2)),0,Colors(nDrug,:),0,'-',.5,1,0,0,2);
end
 xlim([2 30])
legend(hp,Drugs);

%%
table_bg=array2table(av_bg,'VariableNames',{'SubID','SessN','nFile','BlockN','ElecN','offset','slope'});
table_bg.SubID=categorical(table_bg.SubID);
table_bg.Drug=cell(size(table_bg,1),1);
nFiles=unique(table_bg.nFile);
for k=1:length(nFiles)
    table_bg.Drug(table_bg.nFile==nFiles(k))=repmat(all_Drugs(nFiles(k)),sum(table_bg.nFile==nFiles(k)),1);
end
table_bg.Drug=categorical(table_bg.Drug);

table_bg.Elec=cell(size(table_bg,1),1);
for k=1:64
    table_bg.Elec(table_bg.ElecN==k)=repmat(hdr.label(k),sum(table_bg.ElecN==k),1);
end
table_bg.Elec=categorical(table_bg.Elec);
table_bg.Drug=reordercats(table_bg.Drug,[4 1 2 3]);

table_pk=array2table(av_peaks,'VariableNames',{'SubID','SessN','nFile','BlockN','ElecN','peak_freq','peak_amp','peak_width'});
table_pk.SubID=categorical(table_pk.SubID);
uS=unique(table_pk.SubID);
table_pk.Drug=cell(size(table_pk,1),1);
nFiles=unique(table_pk.nFile);
for k=1:length(nFiles)
    table_pk.Drug(table_pk.nFile==nFiles(k))=repmat(all_Drugs(nFiles(k)),sum(table_pk.nFile==nFiles(k)),1);
end
table_pk.Drug=categorical(table_pk.Drug);

table_pk.Elec=cell(size(table_pk,1),1);
for k=1:64
    table_pk.Elec(table_pk.ElecN==k)=repmat(hdr.label(k),sum(table_pk.ElecN==k),1);
end
table_pk.Elec=categorical(table_pk.Elec);
table_pk.Drug=reordercats(table_pk.Drug,[4 1 2 3]);


%% PARAM PLOTS
table_SW=readtable([save_path filesep 'CTET_SWdetection_thr90_allE_P2P_behav_vec_full_v3.txt']);
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.channel=unique(table_SW.Elec);
cfg.center      = 'yes';
layout=ft_prepare_layout(cfg);

Drugs={'PLA','ATM','CIT','MPH'};
cmap=cbrewer('seq','YlOrRd',64); % select a sequential colorscale from yellow to red (64)
cmap3=cbrewer('div','Spectral',64); cmap3=flipud(cmap3);

clus_alpha=0.05;
montecarlo_alpha=0.05;

cfg_neighb=[];
cfg_neighb.method = 'template';
cfg_neighb.layout='biosemi64.lay';
cfg_neighb.channel=unique(table_SW.Elec);
neighbours = ft_prepare_neighbours(cfg_neighb);

redo=1;
totperm=1000;
%% BACKGROUND: PLOTS
figure; set(gcf,'Position',[213         173        1027         805]);
for nDrug=1:4
    subplot(2,4,nDrug)
    temp_topo=[];
    for nE=1:length(layout.label)-2
        temp_topo(nE,1)=mean(table_bg.offset(table_bg.Drug==Drugs{nDrug} & table_bg.Elec==layout.label{nE}));
    end
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    
    
    colormap(cmap3);
    
    title(Drugs{nDrug});
    caxis([2.5 4.5]);
    format_fig;
end

for nDrug=1:4
    subplot(2,4,4+nDrug)
    temp_topo=[];
    for nE=1:length(layout.label)-2
        temp_topo(nE,1)=mean(table_bg.slope(table_bg.Drug==Drugs{nDrug} & table_bg.Elec==layout.label{nE}));
    end
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    
    
    colormap(cmap3);
    
    title(Drugs{nDrug});
    caxis([0.5 2]);
    format_fig;
end

%% PEAKS: PLOTS
table_pkap=table_pk(table_pk.peak_freq>8 & table_pk.peak_freq<12,:);

figure; set(gcf,'Position',[213         173        1027         805]);
for nDrug=1:4
    subplot(2,4,nDrug)
    temp_topo=[];
    for nE=1:length(layout.label)-2
        temp_topo(nE,1)=mean(table_pkap.peak_freq(table_pkap.Drug==Drugs{nDrug} & table_pkap.Elec==layout.label{nE}));
    end
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    
    
    colormap(cmap3);
    
    title(Drugs{nDrug});
    caxis([9.8 10.2]);
    format_fig;
end

for nDrug=1:4
    subplot(2,4,4+nDrug)
    temp_topo=[];
    for nE=1:length(layout.label)-2
        temp_topo(nE,1)=mean(table_pkap.peak_amp(table_pkap.Drug==Drugs{nDrug} & table_pkap.Elec==layout.label{nE}));
    end
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    
    
    colormap(cmap3);
    
    title(Drugs{nDrug});
    caxis([.6 1]);
    format_fig;
end

%% PEAKS: PLOTS
table_pktag=table_pk(table_pk.peak_freq>=24.8 & table_pk.peak_freq<=25.2,:);

figure; set(gcf,'Position',[213         173        1027         805]);
for nDrug=1:4
    subplot(1,4,nDrug)
    temp_topo=[];
    for nE=1:length(layout.label)-2
        temp_topo(nE,1)=mean(table_pktag.peak_amp(table_pktag.Drug==Drugs{nDrug} & table_pktag.Elec==layout.label{nE}));
    end
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    
    
    colormap(cmap3);
    
    title(Drugs{nDrug});
    caxis([0.4 0.75]);
    format_fig;
end

%% BACKGROUND: STATS
if redo==1
    temp_topo_tval=[];
    temp_topo_pval=[];
    % fprintf('%2.0f/%2.0f\n',0,64)
    Slope_est=cell(1,2);
    for nE=1:length(layout.label)-2
        fprintf('%2.0f/%2.0f\n',nE,64)
        sub_table_bg=table_bg((table_bg.Elec==layout.label{nE}),:);
        mdl_byEle{nE}=fitlme(sub_table_bg,'slope~1+BlockN+Drug+(1|SubID)');
        temp_topo_tval(nE,:)=double(mdl_byEle{nE}.Coefficients(2:end,4));
        temp_topo_pval(nE,:)=double(mdl_byEle{nE}.Coefficients(2:end,6));
        
        [real_out, perm_out]=lme_perm_lsneurom(sub_table_bg,'Drug','slope~1+pred+BlockN+(1|SubID)',totperm);
        Slope_est{1}=[Slope_est{1} ; [nE*ones(3,1) real_out]];
        for nDrug=1:3
            Slope_est{2}=[Slope_est{2} ; [nE*ones(totperm,1) perm_out{nDrug} nDrug*ones(totperm,1)]];
        end
    end
    save('../../Tables/model_FOOF_slope_est','Slope_est');
else
    load('../../Tables/model_FOOF_slope_est');
end

%%% Filter clusters
[Slope_clus]=get_clusterperm_lme_lsneurom(Slope_est,clus_alpha,montecarlo_alpha,totperm,neighbours);

if redo==1
    temp_topo_tval=[];
    temp_topo_pval=[];
    % fprintf('%2.0f/%2.0f\n',0,64)
    Offset_est=cell(1,2);
    for nE=1:length(layout.label)-2
        fprintf('%2.0f/%2.0f\n',nE,64)
        sub_table_bg=table_bg((table_bg.Elec==layout.label{nE}),:);
        mdl_byEle{nE}=fitlme(sub_table_bg,'offset~1+BlockN+Drug+(1|SubID)');
        temp_topo_tval(nE,:)=double(mdl_byEle{nE}.Coefficients(2:end,4));
        temp_topo_pval(nE,:)=double(mdl_byEle{nE}.Coefficients(2:end,6));
        
        [real_out, perm_out]=lme_perm_lsneurom(sub_table_bg,'Drug','offset~1+pred+BlockN+(1|SubID)',totperm);
        Offset_est{1}=[Offset_est{1} ; [nE*ones(3,1) real_out]];
        for nDrug=1:3
            Offset_est{2}=[Offset_est{2} ; [nE*ones(totperm,1) perm_out{nDrug} nDrug*ones(totperm,1)]];
        end
    end
    save('../../Tables/model_FOOF_offset_est','Offset_est');
else
    load('../../Tables/model_FOOF_offset_est');
end

%%% Filter clusters
[Offset_clus]=get_clusterperm_lme_lsneurom(Offset_est,clus_alpha,montecarlo_alpha,totperm,neighbours);

%% PEAK: STATS
if redo==1
    temp_topo_tval=[];
    temp_topo_pval=[];
    % fprintf('%2.0f/%2.0f\n',0,64)
    AlphaFreq_est=cell(1,2);
    for nE=1:length(layout.label)-2
        fprintf('%2.0f/%2.0f\n',nE,64)
        sub_table_pkap=table_pkap((table_pkap.Elec==layout.label{nE}),:);
        mdl_byEle{nE}=fitlme(sub_table_pkap,'peak_freq~1+BlockN+Drug+(1|SubID)');
        temp_topo_tval(nE,:)=double(mdl_byEle{nE}.Coefficients(2:end,4));
        temp_topo_pval(nE,:)=double(mdl_byEle{nE}.Coefficients(2:end,6));
        
        [real_out, perm_out]=lme_perm_lsneurom(sub_table_pkap,'Drug','peak_freq~1+pred+BlockN+(1|SubID)',totperm);
        AlphaFreq_est{1}=[AlphaFreq_est{1} ; [nE*ones(3,1) real_out]];
        for nDrug=1:3
            AlphaFreq_est{2}=[AlphaFreq_est{2} ; [nE*ones(totperm,1) perm_out{nDrug} nDrug*ones(totperm,1)]];
        end
    end
    save('../../Tables/model_FOOF_AlphaFreq_est','AlphaFreq_est');
else
    load('../../Tables/model_FOOF_AlphaFreq_est');
end

%%% Filter clusters
[AlphaFreq_clus]=get_clusterperm_lme_lsneurom(AlphaFreq_est,clus_alpha,montecarlo_alpha,totperm,neighbours);

if redo==1
    temp_topo_tval=[];
    temp_topo_pval=[];
    % fprintf('%2.0f/%2.0f\n',0,64)
    AlphaAmp_est=cell(1,2);
    for nE=1:length(layout.label)-2
        fprintf('%2.0f/%2.0f\n',nE,64)
        sub_table_pkap=table_pkap((table_pkap.Elec==layout.label{nE}),:);
        mdl_byEle{nE}=fitlme(sub_table_pkap,'peak_amp~1+BlockN+Drug+(1|SubID)');
        temp_topo_tval(nE,:)=double(mdl_byEle{nE}.Coefficients(2:end,4));
        temp_topo_pval(nE,:)=double(mdl_byEle{nE}.Coefficients(2:end,6));
        
        [real_out, perm_out]=lme_perm_lsneurom(sub_table_pkap,'Drug','peak_amp~1+pred+BlockN+(1|SubID)',totperm);
        AlphaAmp_est{1}=[AlphaAmp_est{1} ; [nE*ones(3,1) real_out]];
        for nDrug=1:3
            AlphaAmp_est{2}=[AlphaAmp_est{2} ; [nE*ones(totperm,1) perm_out{nDrug} nDrug*ones(totperm,1)]];
        end
    end
    save('../../Tables/model_FOOF_AlphaAmp_est','AlphaAmp_est');
else
    load('../../Tables/model_FOOF_AlphaAmp_est');
end

%%% Filter clusters
[AlphaAmp_clus]=get_clusterperm_lme_lsneurom(AlphaAmp_est,clus_alpha,montecarlo_alpha,totperm,neighbours);

%% BACKGROUND: CLUSTER
cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);
limNumClus=0;
limMax=10;%max(max(abs(temp_topo_tval)));
figure; set(gcf,'Position',[213         173        1027         805]);
for nDrug=1:3
    subplot(1,3,nDrug)
    
    temp_topo=Slope_est{1}(Slope_est{1}(:,5)==nDrug,3);
    temp_topo2=zeros(size(temp_topo));
    temp_topo3=zeros(size(temp_topo));
    temp_clus=Slope_clus{nDrug};
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
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    ft_plot_lay_me(layout, 'chanindx',1:length(layout.label),'pointsymbol','o','pointcolor',[1 1 1]*0.7,'pointsize',6,'box','no','label','no')
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
    if nDrug==3
        h=colorbar;
        set(h,'Position',[0.93 0.4 0.02 0.2])
    end
    colormap(cmap2);
    
    title(Drugs{nDrug+1})
end

figure; set(gcf,'Position',[213         173        1027         805]);
for nDrug=1:3
    subplot(1,3,nDrug)
    
    temp_topo=Offset_est{1}(Offset_est{1}(:,5)==nDrug,3);
    temp_topo2=zeros(size(temp_topo));
    temp_topo3=zeros(size(temp_topo));
    temp_clus=Offset_clus{nDrug};
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
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    ft_plot_lay_me(layout, 'chanindx',1:length(layout.label),'pointsymbol','o','pointcolor',[1 1 1]*0.7,'pointsize',6,'box','no','label','no')
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
    if nDrug==3
        h=colorbar;
        set(h,'Position',[0.93 0.4 0.02 0.2])
    end
    colormap(cmap2);
    
    title(Drugs{nDrug+1})
end

%% PEAK: CLUSTER
cmap2=cbrewer('div','RdBu',64); % select a sequential colorscale from yellow to red (64)
cmap2=flipud(cmap2);
limNumClus=0;
limMax=10;%max(max(abs(temp_topo_tval)));
figure; set(gcf,'Position',[213         173        1027         805]);
for nDrug=1:3
    subplot(1,3,nDrug)
    
    temp_topo=AlphaFreq_est{1}(AlphaFreq_est{1}(:,5)==nDrug,3);
    temp_topo2=zeros(size(temp_topo));
    temp_topo3=zeros(size(temp_topo));
    temp_clus=AlphaFreq_clus{nDrug};
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
    ft_plot_lay_me(layout, 'chanindx',1:length(layout.label),'pointsymbol','o','pointcolor',[1 1 1]*0.7,'pointsize',6,'box','no','label','no')
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
    if nDrug==3
        h=colorbar;
        set(h,'Position',[0.93 0.4 0.02 0.2])
    end
    colormap(cmap2);
    
    title(Drugs{nDrug+1})
end

figure; set(gcf,'Position',[213         173        1027         805]);
for nDrug=1:3
    subplot(1,3,nDrug)
    
    temp_topo=AlphaAmp_est{1}(AlphaAmp_est{1}(:,5)==nDrug,3);
    temp_topo2=zeros(size(temp_topo));
    temp_topo3=zeros(size(temp_topo));
    temp_clus=AlphaAmp_clus{nDrug};
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
    simpleTopoPlot_ft(temp_topo, layout,'on',[],0,1);
    ft_plot_lay_me(layout, 'chanindx',1:length(layout.label),'pointsymbol','o','pointcolor',[1 1 1]*0.7,'pointsize',6,'box','no','label','no')
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
    if nDrug==3
        h=colorbar;
        set(h,'Position',[0.93 0.4 0.02 0.2])
    end
    colormap(cmap2);
    
    title(Drugs{nDrug+1})
end


