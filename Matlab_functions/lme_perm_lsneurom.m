function [real_out, perm_out]=lme_perm_lsneurom(table,predictor,formula,totperm)
% table=GO_table;
%
% formula='RT~1+pred+Task+(1|SubID)';
% permvars={'Task','SubID'};
% totperm=100;

% run real model
if strcmp(predictor,'Drug')
    eval(sprintf('table.pred=table.%s;',predictor));
    model= fitlme(table,formula);
    real_out=[];
    for nDrug=1:3
        real_out=[real_out ; double(model.Coefficients(2+nDrug,2)) double(model.Coefficients(2+nDrug,4)) double(model.Coefficients(2+nDrug,6)) nDrug];
    end
    
    uniqueIDs=unique(table.SubID);
    uniqueBlocks=unique(table.BlockN);
    perm_out={nan(totperm,4),nan(totperm,4),nan(totperm,4)};
    fprintf('%4.0f/%4.0f\n',0,totperm)
    for np=1:totperm
        pred_perm=table.pred;
        for nT=1:length(uniqueBlocks)
            for nS=1:length(uniqueIDs)
                idx=find(table.SubID==uniqueIDs(nS) & table.BlockN==uniqueBlocks(nT));
                temp=pred_perm(idx);
                pred_perm(idx)=temp(randperm(length(temp)));
            end
        end
        
        table2=table;
        table2.pred=pred_perm;
        model= fitlme(table2,formula);
        for nDrug=1:3
            perm_out{nDrug}(np,:)=[double(model.Coefficients(2+nDrug,2)) double(model.Coefficients(2+nDrug,4)) double(model.Coefficients(2+nDrug,6)) np];
        end
        fprintf('\b\b\b\b\b\b\b\b\b\b%4.0f/%4.0f\n',np,totperm)
    end
    fprintf('\n');
    
else
    eval(sprintf('table.pred=table.%s;',predictor));
    model= fitlme(table,formula);
    real_out=[double(model.Coefficients(3,2)) double(model.Coefficients(3,4)) double(model.Coefficients(3,6))];
    
    uniqueIDs=unique(table.SubID);
    uniqueBlocks=unique(table.BlockN);
    perm_out=nan(totperm,4);
    fprintf('%4.0f/%4.0f\n',0,totperm)
    for np=1:totperm
        pred_perm=table.pred;
        for nT=1:length(uniqueBlocks)
            for nS=1:length(uniqueIDs)
                idx=find(table.SubID==uniqueIDs(nS) & table.BlockN==uniqueBlocks(nT));
                temp=pred_perm(idx);
                pred_perm(idx)=temp(randperm(length(temp)));
            end
        end
        
        table2=table;
        table2.pred=pred_perm;
        model= fitlme(table2,formula);
        perm_out(np,:)=[double(model.Coefficients(3,2)) double(model.Coefficients(3,4)) double(model.Coefficients(3,6)) np];
        fprintf('\b\b\b\b\b\b\b\b\b\b%4.0f/%4.0f\n',np,totperm)
    end
    fprintf('\n');
    
end