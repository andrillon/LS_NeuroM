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

prob_files=[%{'05_ctet_session4_MPH.bdf'}
    %{'15_ctet_session2_MPH.bdf'}
    %%{'28_ctet_session3_PLA.bdf'}
    %%{'28_ctet_session4_CIT.bdf'}
    %{'33_ctet_session2_ATM.bdf'}
    {'36_ctet_session4_CIT.bdf'}];
%%
res_mat=[];
drug_cond=[];
problematic_files=[];
nc=0;
for nF=1:length(files)
    File_Name=files(nF).name;
    fprintf('... processing %s\n',File_Name);
    if ~ismember(File_Name,prob_files)
        continue;
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

end

%%


