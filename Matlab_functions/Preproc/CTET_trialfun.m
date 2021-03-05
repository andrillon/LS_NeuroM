function trl = CTET_trialfun(cfg);

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

for i=1:size(subtable,1)
    % add this to the trl definition
    begsample     = subtable.Sample(i) - cfg.trialdef.prestim*hdr.Fs;
    endsample     = subtable.Sample(i) + cfg.trialdef.poststim*hdr.Fs - 1;
    offset        = -cfg.trialdef.prestim*hdr.Fs;
    trl(end+1, :) = [round([begsample endsample offset])];
end
