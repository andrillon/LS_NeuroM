%%
clear all
close all

path_fieldtrip='/Users/tand0009/Work/local/fieldtrip/';
addpath(path_fieldtrip);
ft_defaults;

path_LSCPtools='/Users/tand0009/WorkGit/LSCPtools/';
addpath(genpath(path_LSCPtools));

data_path='/Volumes/tLab_BackUp1/Monash/CTET_Dockree/EEG_CTET/';
% data_path='/Volumes/shared/R-MNHS-SPP/Bellgrove-data/Jess Barnes EEG Backup Data/EEG_CTET/';
files=dir([data_path filesep '*.bdf']);

wrong_events=[{'28_ctet_session3_PLA.bdf'},{'28_ctet_session4_CIT.bdf'}];
%%
res_mat=[];
drug_cond=[];
problematic_files=[];
resblock_mat=[];
drugblock_cond=[];
nc=0;
for nF=1:length(files)
    File_Name=files(nF).name;
    if ismember(File_Name,wrong_events)
        fprintf('... skipping %s\n',File_Name);
        continue;
    else
        fprintf('... processing %s\n',File_Name);
    end
    septag=findstr(File_Name,'_');
    SubN=str2num(File_Name(1:septag(1)-1));
    SessN=str2num(File_Name(septag(3)-1));
    DrugC=(File_Name(septag(3)+1:septag(3)+3));
    
    % data=ft_read_data('/Volumes/tLab_BackUp1/Monash/CTET_Dockree/EEG_CTET/01_ctet_session1_ATM.bdf');
    hdr=ft_read_header([data_path filesep File_Name]);
    events=ft_read_event([data_path filesep File_Name]);

    my_events=events(find_trials({events.type},'STATUS'));
    unique_values=unique([my_events.value]);
    fprintf('... %g events found:',length([my_events.value]));
    for nU=1:length(unique_values)
        fprintf(' %g (n=%g) ',unique_values(nU),sum([my_events.value]==unique_values(nU)))
    end
    fprintf('\n');


    %%% find blocks
    all_idx=[my_events.sample];
    all_val=[my_events.value];
    blocks_boundaries=all_idx(find(diff(all_idx/hdr.Fs)>5 & all_val(2:end)==1)+1);
    if strcmp(File_Name,'33_ctet_session2_ATM.bdf') || strcmp(File_Name,'36_ctet_session4_CIT.bdf')
        blocks_boundaries=all_idx(find(diff(all_idx/hdr.Fs)>3 & all_val(2:end)==1)+1);
    end
    if strcmp(File_Name,'15_ctet_session2_MPH.bdf') && length(blocks_boundaries)==11
        blocks_boundaries(end)=[];
    end
    if all_val(1)==1 && length(blocks_boundaries)==9
        blocks_boundaries=[all_idx(1) blocks_boundaries];
    elseif length(blocks_boundaries)==9
        blocks_boundaries=[all_idx(1) blocks_boundaries];
    end
    fprintf('... %g block transitions found\n',length(blocks_boundaries));
    if length(blocks_boundaries)~=10
        nc=nc+1;
       problematic_files{nc}= File_Name;
    end
    all_val(ismember(all_idx,blocks_boundaries))=100;
    
    %%% find relevant events
    resp_idx=[my_events(all_val==1).sample];
    
    targets_idx=[my_events(all_val==max(unique_values)).sample];
    nontargets_idx=[my_events(all_val==20).sample];
    stim_idx=[my_events(all_val==max(unique_values) | all_val==20).sample];
    stim_val=[my_events(all_val==max(unique_values) | all_val==20).value];
    
    %%% Gather info per targets
    this_behav=nan(length(stim_idx),7);
    for nSt=1:length(all_idx)
        this_FA=NaN;
        this_RT=NaN;
        if all_val(nSt)==max(unique_values) % we have a target
            count=0; resprec=0;
            while nSt+count<length(all_val) && all_val(nSt+count)~=1
                count=count+1;
            end
            if all_val(nSt+count)==1
                this_RT=(all_idx(nSt+count)-all_idx(nSt))/hdr.Fs;
            else
                this_RT=NaN;
            end
        elseif all_val(nSt)==20 % we have a non-target
            if nSt<length(all_val) & nSt>5
                if all_val(nSt+1)==1 && max(all_val(nSt+(-5:-1)))~=max(unique_values)
                    this_FA=1;
                else
                    this_FA=0;
                end
            end
        else
            continue;
        end
        %%%% sort things out
        this_block=max(find(all_idx(nSt)>blocks_boundaries));
        if isempty(this_block)
            this_block=NaN;
        end
        if ~isnan(this_RT)
            if this_RT>2.6
                this_Miss=1;
            else
                this_Miss=0;
            end
        elseif isnan(this_RT) && all_val(nSt)==max(unique_values) && nSt~=length(all_idx)
            this_Miss=1;
        else
            this_Miss=NaN;
        end
        this_behav(find(stim_idx==all_idx(nSt)),:)=[find(stim_idx==all_idx(nSt)) stim_idx(find(stim_idx==all_idx(nSt))) this_block all_val(nSt)==max(unique_values) this_RT this_FA this_Miss];
    end
    this_behav(this_behav(:,5)<mode(diff(nontargets_idx)/hdr.Fs),7)=1;
    this_behav(find(this_behav(:,5)<mode(diff(nontargets_idx)/hdr.Fs))-1,6)=1;
    this_cleanRT=this_behav(:,5);
    this_cleanRT(this_cleanRT>2.6 | this_cleanRT<mode(diff(nontargets_idx)/hdr.Fs))=NaN;
    res_mat=[res_mat; [SubN*ones(size(this_behav,1),1) SessN*ones(size(this_behav,1),1) this_behav this_cleanRT]];
    drug_cond=[drug_cond ; repmat({DrugC},size(this_behav,1),1)];

    for nbl=1:10
        resblock_mat=[resblock_mat; [SubN SessN nbl nanmean(this_behav(this_behav(:,3)==nbl & this_behav(:,4)==0,6)==0) nanmean(this_behav(this_behav(:,3)==nbl & this_behav(:,4)==0,6)==1) ...
            nanmean(this_behav(this_behav(:,3)==nbl & this_behav(:,4)==1,7)==0) nanmean(this_behav(this_behav(:,3)==nbl & this_behav(:,4)==1,7)==1) nanmean(this_cleanRT(this_behav(:,3)==nbl & this_behav(:,4)==1))]];
        drugblock_cond=[drugblock_cond ; {DrugC}];
    end
    
end

%%
table=array2table(res_mat,'VariableNames',{'SubID','SessN','TrialN','Sample','BlockN','StimType','All_RT','FA','Miss','Hit_RT'});
table.Treatment=drug_cond(~isnan(res_mat(:,3)));
% table.Treatment(match_str(table.Treatment,'MPH_run-at-wrong_Hz-rate'))={'MPH'};
% table.Treatment(match_str(table.Treatment,'ATM_run-at-wrong-Hz-rate'))={'ATM'};
table.SubID=categorical(table.SubID);
table.SessN=categorical(table.SessN);
table.Treatment=categorical(table.Treatment);
table.Treatment=reordercats(table.Treatment,[4,1,2,3]);

table.TrialN=ordinal(table.TrialN);
table.BlockN=ordinal(table.BlockN);

mdl1=fitlme(table,'Hit_RT~1+Treatment+(1|SubID)');
mdl2=fitlme(table,'FA~1+Treatment+(1|SubID)');
mdl3=fitlme(table,'Miss~1+Treatment+(1|SubID)');

writetable(table,'/Users/tand0009/Data/CTET_Dockree/CTET_behav_res.txt');

%%
table2=array2table(resblock_mat,'VariableNames',{'SubID','SessN','BlockN','CR','FA','Hit','Miss','Hit_RT'});
table2.Treatment=drugblock_cond;
table2.SubID=categorical(table2.SubID);
table2.SessN=categorical(table2.SessN);
table2.Treatment=categorical(table2.Treatment);
table2.Treatment=reordercats(table2.Treatment,[4,1,2,3]);

table2.BlockN=ordinal(table2.BlockN);

mdl1=fitlme(table2,'Hit_RT~1+Treatment+(1|SubID)');
mdl2=fitlme(table2,'FA~1+Treatment+(1|SubID)');
mdl3=fitlme(table2,'Miss~1+Treatment+(1|SubID)');

writetable(table2,'/Users/tand0009/Data/CTET_Dockree/CTET_behav_resblock.txt');
%%
uniqueIDs=unique(table.SubID);
uniqueDrugs=unique(table.Treatment);
labelsDrugs=[];
for nSubj=1:length(uniqueIDs)
    for nDrug=1:length(uniqueDrugs)
        Hit_RT(nSubj,nDrug)=nanmean(table.Hit_RT(table.SubID==uniqueIDs(nSubj) & table.Treatment==uniqueDrugs(nDrug)));
        FA(nSubj,nDrug)=nanmean(table.FA(table.SubID==uniqueIDs(nSubj) & table.Treatment==uniqueDrugs(nDrug)));
        Miss(nSubj,nDrug)=nanmean(table.Miss(table.SubID==uniqueIDs(nSubj) & table.Treatment==uniqueDrugs(nDrug)));
        labelsDrugs{nDrug}=char(uniqueDrugs(nDrug));
    end
end

Colors=[[1 1 1]*0.5; [253,174,97]/256 ; [44,123,182]/256; [215,25,28]/256];
figure; set(gcf,'Position',[82         600        1158         378]);
for nPlot=1:3
    subplot(1,3,nPlot); format_fig; hold on;
    if nPlot==1
        temp=Hit_RT;
    elseif nPlot==2
        temp=FA;
    else
        temp=Miss;
    end
    format_fig;
    for nDrug=1:length(uniqueDrugs)
        simpleBarPlot(nDrug,temp(:,nDrug),Colors(nDrug,:),0.8,'k',[],3);
    end
    set(gca,'XTick',1:length(uniqueDrugs),'XTickLabel',labelsDrugs);
    if nPlot==1
    ylim([1.2 1.5])
        title('Hit RT')
    elseif nPlot==2
        title('FA')
    else
        title('Miss')
    end
end

%%
load('/Users/tand0009/Data/CTET_Dockree/VAS_alertness_exclusions.mat');
figure; %set(gcf,'Position',[82         600        1158         378]);
for nTime=1:3
    subplot(2,1,1); format_fig
    for nDrug=1:4
        eval(sprintf('temp=VASalertnessexclusions.%s_VS%g;',labelsDrugs{nDrug},nTime));
        simpleBarPlot(nDrug+(nTime-2)*0.2,temp,Colors(nDrug,:),0.15,'k',[],3);
    end
        set(gca,'XTick',1:length(uniqueDrugs),'XTickLabel',labelsDrugs);

        subplot(2,1,2); format_fig
        eval(sprintf('temp0=VASalertnessexclusions.%s_VS%g;',labelsDrugs{1},nTime));
    for nDrug=2:4
        eval(sprintf('temp=VASalertnessexclusions.%s_VS%g;',labelsDrugs{nDrug},nTime));
        simpleBarPlot(nDrug+(nTime-2)*0.2,temp./temp0,Colors(nDrug,:),0.15,'k',[],3);
    end
        set(gca,'XTick',2:length(uniqueDrugs),'XTickLabel',labelsDrugs(2:4));
line(xlim,[1 1],'LineStyle','--','Color','k');
end



