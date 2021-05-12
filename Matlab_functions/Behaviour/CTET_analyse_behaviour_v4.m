%%
clear all
close all
run ../localdef.m

addpath(path_fieldtrip);
ft_defaults;
addpath(genpath(path_LSCPtools));

files=dir([data_path filesep '*.bdf']);

wrong_events=[{'28_ctet_session3_PLA.bdf'},{'28_ctet_session4_CIT.bdf'}];

%%
redo=0;
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
    
    if ~ismember(SubN,ListSubjectsID)
        fprintf('... %s not in subject list\n',File_Name);
        continue;
    end
    if redo==1 || exist([save_path filesep 'CTET_behav_' File_Name(1:end-4) '.txt'])==0
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
        %%% Gather info per targets
        this_behav=nan(length(stim_idx),9);
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
                    count=0; resprec=0;
                    while nSt+count<length(all_val) && all_val(nSt+count)~=1
                        count=count+1;
                    end
                    if all_val(nSt+count)==1
                        this_RT=(all_idx(nSt+count)-all_idx(nSt))/hdr.Fs;
                    else
                        this_RT=NaN;
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
            this_behav(find(stim_idx==all_idx(nSt)),:)=[find(stim_idx==all_idx(nSt)) stim_idx(find(stim_idx==all_idx(nSt))) this_block all_val(nSt)==max(unique_values) this_RT NaN NaN NaN NaN];
        end
        this_behav(:,8)=[diff(this_behav(:,2))/hdr.Fs ; NaN];
        this_behav(this_behav(:,5)>this_behav(:,8),5)=NaN;
        
        respidx=find(~isnan(this_behav(:,5)));
        this_behav(this_behav(:,4)==0,6)=1;
        this_behav(this_behav(:,4)==1,7)=0;
        for m=1:length(respidx)
            %         post_resp=this_behav(this_behav(:,2)>this_behav(respidx(m),2),3);
            %         post_resp=post_resp(1);
            %         if this_behav(respidx(m),3)~=post_resp
            %             this_behav(respidx(m),5:9)=NaN;
            %             continue;
            %         end
            if this_behav(respidx(m),5)>5*mode(this_behav(this_behav(:,8)>min(this_behav(:,8))+0.1,8))
                this_behav(respidx(m),5:9)=NaN;
                %             continue;
            end
            if this_behav(respidx(m),4)==1 % resp on target
                if this_behav(respidx(m),5)>min(this_behav(:,8)) %&  this_behav(respidx(m),5)<mode(this_behav(this_behav(:,8)>min(this_behav(:,8))+0.1,8))+2*min(this_behav(:,8))
                    this_behav(respidx(m),7)=1;
                    this_behav(respidx(m),9)=this_behav(respidx(m),5);
                else
                    this_behav(respidx(m),7)=0;
                    this_behav(respidx(m)-1,6)=0;
                    this_behav(respidx(m)-1,9)=this_behav(respidx(m),5)+this_behav(respidx(m)-1,8);
                end
            elseif this_behav(respidx(m),4)==0 % resp on non target
                if this_behav(respidx(m)-1,4)==1 %&& this_behav(respidx(m)-1,3)==this_behav(respidx(m),3)
                    this_behav(respidx(m)-1,7)=1;
                    this_behav(respidx(m)-1,9)=this_behav(respidx(m),5)+this_behav(respidx(m)-1,8);
                else
                    if this_behav(respidx(m)-2,4)==1 %&& this_behav(respidx(m)-2,3)==this_behav(respidx(m),3)
                        this_behav(respidx(m)-2,7)=1;
                        this_behav(respidx(m)-2,9)=this_behav(respidx(m),5)+this_behav(respidx(m)-1,8)+this_behav(respidx(m)-2,8);
                    else
                        this_behav(respidx(m)-1,6)=0;
                        this_behav(respidx(m)-1,9)=this_behav(respidx(m),5)+this_behav(respidx(m)-1,8);
                    end
                end
            end
        end
        
        this_table=array2table(this_behav,...
            'VariableNames',{'TrialN','Sample','BlockN','StimType','rawRT','corrNT','corrTG','ITI','RT'});
        writetable(this_table,[save_path filesep 'CTET_behav_' File_Name(1:end-4) '.txt']);
    else
        this_table=readtable([save_path filesep 'CTET_behav_' File_Name(1:end-4) '.txt']);
        this_behav=table2array(this_table);
    end
    res_mat=[res_mat; [SubN*ones(size(this_behav,1),1) SessN*ones(size(this_behav,1),1) this_behav]];
    drug_cond=[drug_cond ; repmat({DrugC},size(this_behav,1),1)];
    
    for nbl=1:10
        [dprime,crit]=calc_dprime(this_behav(this_behav(:,3)==nbl & this_behav(:,4)==0,6)==1,this_behav(this_behav(:,3)==nbl & this_behav(:,4)==1,7)==0);
        resblock_mat=[resblock_mat; [SubN SessN nbl nanmean(this_behav(this_behav(:,3)==nbl & this_behav(:,4)==0,6)==1) nanmean(this_behav(this_behav(:,3)==nbl & this_behav(:,4)==0,6)==0) ...
            nanmean(this_behav(this_behav(:,3)==nbl & this_behav(:,4)==1,7)==1) nanmean(this_behav(this_behav(:,3)==nbl & this_behav(:,4)==1,7)==0) nanmean(this_behav(this_behav(:,3)==nbl & this_behav(:,4)==1 & this_behav(:,7)==1,9)) nanstd(this_behav(this_behav(:,3)==nbl & this_behav(:,4)==1 & this_behav(:,7)==1,9)) ...
            dprime crit]];
        drugblock_cond=[drugblock_cond ; {DrugC}];
    end
    
end

%%
table=array2table(res_mat,'VariableNames',{'SubID','SessN','TrialN','Sample','BlockN','StimType','rawRT','corrNT','corrTG','ITI','RT'});
table.Treatment=drug_cond(~isnan(res_mat(:,3)));
% table.Treatment(match_str(table.Treatment,'MPH_run-at-wrong_Hz-rate'))={'MPH'};
% table.Treatment(match_str(table.Treatment,'ATM_run-at-wrong-Hz-rate'))={'ATM'};
table.SubID=categorical(table.SubID);
table.SessN=categorical(table.SessN);
table.Treatment=categorical(table.Treatment);
table.Treatment=reordercats(table.Treatment,[4,1,2,3]);

% table.TrialN=ordinal(table.TrialN);
% table.BlockN=ordinal(table.BlockN);

mdl1=fitlme(table,'RT~1+Treatment+(1|SubID)');
% mdl2=fitlme(table,'FA~1+Treatment+(1|SubID)');
% mdl3=fitlme(table,'Miss~1+Treatment+(1|SubID)');

writetable(table,[save_path filesep 'CTET_behav_res.txt']);

myS=unique(table.SubID);
myD=unique(table.Treatment);
for nS=1:length(myS)
    for nD=1:length(myD)
        nSession(nS,nD)=sum(table.SubID==myS(nS) & table.Treatment==myD(nD))~=0;
    end
end
FullSubID=myS(sum(nSession,2)==4);
table2=table(ismember(table.SubID,FullSubID),:);
writetable(table2,[save_path filesep 'CTET_behav_res_full.txt']);

%%
table2=array2table(resblock_mat,'VariableNames',{'SubID','SessN','BlockN','CR','FA','Hit','Miss','Hit_RT','STD_RT','dprime','crit'});
table2.Treatment=drugblock_cond;
table2.SubID=categorical(table2.SubID);
table2.SessN=categorical(table2.SessN);
table2.Treatment=categorical(table2.Treatment);
table2.Treatment=reordercats(table2.Treatment,[4,1,2,3]);

table2.BlockN=ordinal(table2.BlockN);

mdl1=fitlme(table2,'Hit_RT~1+Treatment+(1|SubID)');
mdl2=fitlme(table2,'FA~1+Treatment+(1|SubID)');
mdl3=fitlme(table2,'Miss~1+Treatment+(1|SubID)');

writetable(table2,[save_path filesep 'CTET_behav_resblock.txt']);

myS=unique(table2.SubID);
myD=unique(table2.Treatment);
for nS=1:length(myS)
    for nD=1:length(myD)
        nSession(nS,nD)=sum(table2.SubID==myS(nS) & table2.Treatment==myD(nD))~=0;
    end
end
FullSubID=myS(sum(nSession,2)==4);
table3=table2(ismember(table2.SubID,FullSubID),:);
writetable(table3,[save_path filesep 'CTET_behav_resblock_full.txt']);


