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
mdl0_byEle=[];
mdl_byEle=[];
for nBehav=1:3
    
    if nBehav==1
        Behav_Var='RT';
        sub_table_SW2=table_SW2(~isnan(table_SW2.RT) & table_SW2.StimType==1,:);
    elseif nBehav==2
        Behav_Var='Miss';
        sub_table_SW2=table_SW2(~isnan(table_SW2.corrTG) & table_SW2.StimType==1,:);
        sub_table_SW2.Miss=1-sub_table_SW2.corrTG;
    elseif nBehav==3
        Behav_Var='FA';
        sub_table_SW2=table_SW2(~isnan(table_SW2.corrNT) & table_SW2.StimType==0,:);
        sub_table_SW2.FA=1-sub_table_SW2.corrNT;
    end
    
    
    
    % fprintf('%2.0f/%2.0f\n',0,64)
    for nE=1:length(layout.label)-2
        fprintf('B%g: %2.0f/%2.0f\n',nBehav,nE,64)
        if strcmp(Behav_Var,'Miss') || strcmp(Behav_Var,'FA')
            mdl0_byEle{nBehav,nE}=fitglme(sub_table_SW2,sprintf('%s~1+Drug+(1|SubID)+(1|Drug)',Behav_Var),'Distribution','binomial');
            mdl_byEle{nBehav,nE}=fitglme(sub_table_SW2,sprintf('%s~1+Drug+%s+(1|SubID)+(1|Drug)',Behav_Var,layout.label{nE}),'Distribution','binomial');
        else
            mdl0_byEle{nBehav,nE}=fitlme(sub_table_SW2,sprintf('%s~1+Drug+(1|SubID)+(1|Drug)',Behav_Var));
            mdl_byEle{nBehav,nE}=fitlme(sub_table_SW2,sprintf('%s~1+Drug+%s+(1|SubID)+(1|Drug)',Behav_Var,layout.label{nE}));
        end
        topo_tval(nBehav,nE,:)=double(mdl_byEle{nBehav,nE}.Coefficients(2:end,4));
        topo_pval(nBehav,nE,:)=double(mdl_byEle{nBehav,nE}.Coefficients(2:end,6));
        try
            table_comp=compare(mdl0_byEle{nBehav,nE},mdl_byEle{nBehav,nE});
            topo_comp_Xval(nBehav,nE,:)=[table_comp.LRStat(2) table_comp.pValue(2)];
        catch
            topo_comp_Xval(nBehav,nE,:)=nan(1,2);
        end
    end
end


