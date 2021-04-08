function trl = CTET_blockfun(cfg);

% this function requires the following fields to be specified
% cfg.dataset
% cfg.trialdef.eventtype
% cfg.trialdef.eventvalue
% cfg.trialdef.prestim
% cfg.trialdef.poststim

hdr=ft_read_header(cfg.dataset);
trl = [];

% clean events
table=cfg.table;
subtable=table(table.SubID==cfg.SubID & table.SessN==cfg.SessN,:);

for i=1:10 % expecting 10 blocks
    temp_trial=subtable.Sample(ismember(subtable.BlockN,num2str(i)));
    first_trial=temp_trial(1);
    last_trial=temp_trial(end);
    
    % add this to the trl definition
    begsample     = first_trial - cfg.trialdef.prestim*hdr.Fs;
    endsample     = last_trial + cfg.trialdef.poststim*hdr.Fs - 1;
    offset        = -cfg.trialdef.prestim*hdr.Fs;
    trl(end+1, :) = [round([begsample endsample offset])];
end
