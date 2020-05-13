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

%%
res_mat=[];
drug_cond=[];
for nF=1:length(files)
    File_Name=files(nF).name;
    fprintf('... processing %s\n',File_Name);
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


    %%%
    % find STATUS events
    
    all_idx=[my_events.sample];
    all_val=[my_events.value];
    targets_idx=[my_events([my_events.value]==max(unique_values)).sample];
    nontargets_idx=[my_events([my_events.value]==20).sample];
    resp_idx=[my_events([my_events.value]==1).sample];
    stim_idx=[my_events([my_events.value]==max(unique_values) | [my_events.value]==20).sample];
    stim_val=[my_events([my_events.value]==max(unique_values) | [my_events.value]==20).value];
    
    %%% Gather info per targets
    behav_RT=nan(length(targets_idx),1);
    behav_Lapse=nan(length(targets_idx),1);
    behav_Impulse=nan(length(targets_idx),1);
    for nTa=1:length(targets_idx)
        this_target_idx=targets_idx(nTa);
        count=0; prev_nontarget_idx=NaN;
        while stim_val(find(ismember(stim_idx,targets_idx(nTa)))-count)~=20 && find(ismember(stim_idx,targets_idx(nTa)))-count>1
            count=count+1;
        end
        if stim_val(find(ismember(stim_idx,targets_idx(nTa)))-count)==20
            prev_nontarget_idx=stim_idx(find(ismember(stim_idx,targets_idx(nTa)))-count);
        end
        
        count=0; next_nontarget_idx=NaN;
        while stim_val(find(ismember(stim_idx,targets_idx(nTa)))+count)~=20 && (find(ismember(stim_idx,targets_idx(nTa)))+count)<length(stim_val)
            count=count+1;
        end
        if stim_val(find(ismember(stim_idx,targets_idx(nTa)))+count)==20
            next_nontarget_idx=stim_idx(find(ismember(stim_idx,targets_idx(nTa)))+count);
        end
        
        count=0; next_resp_idx=NaN;
        while all_val(find(ismember(all_idx,targets_idx(nTa)))+count)~=1 && (find(ismember(all_idx,targets_idx(nTa)))+count)<length(all_val)
            count=count+1;
        end
        if all_val(find(ismember(all_idx,targets_idx(nTa)))+count)==1
            next_resp_idx=all_idx(find(ismember(all_idx,targets_idx(nTa)))+count);
        end
        
        %%% RT and lapses
        this_RT=(next_resp_idx-this_target_idx)/hdr.Fs;
        behav_RT(nTa)=this_RT;
        if this_RT<2.6
            if this_RT<mode(diff(nontargets_idx)/hdr.Fs)
                behav_Impulse(nTa)=1;
            else
                behav_Impulse(nTa)=0;
            end
            behav_Lapse(nTa)=0;
        elseif this_RT>=2.6 && ~isnan(this_RT)
            behav_Lapse(nTa)=1;
            behav_RT(nTa)=NaN;
            behav_Impulse(nTa)=0;
        end
%              if this_RT<-20
%                 pause;
%             end
    end
    res_mat=[res_mat; [SubN*ones(size(behav_Lapse,1),1) SessN*ones(size(behav_Lapse,1),1) behav_Lapse behav_RT behav_Impulse (1:size(behav_Lapse,1))']];
    drug_cond=[drug_cond ; repmat({DrugC},size(behav_Lapse,1),1)];
end

%%
res_mat2=res_mat(~isnan(res_mat(:,3)),:);
table=array2table(res_mat2,'VariableNames',{'SubID','SessN','Lapse','RT','Impulse','TrialN'});
table.Treatment=drug_cond(~isnan(res_mat(:,3)));
% table.Treatment(match_str(table.Treatment,'MPH_run-at-wrong_Hz-rate'))={'MPH'};
% table.Treatment(match_str(table.Treatment,'ATM_run-at-wrong-Hz-rate'))={'ATM'};
table.Treatment=categorical(table.Treatment);
table.Treatment=reordercats(table.Treatment,[4,1,2,3]);
