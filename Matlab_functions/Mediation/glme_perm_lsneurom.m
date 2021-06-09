function [real_out, perm_out,out_pred_perm]=glme_perm_lsneurom(table,predictor,formula,totperm,all_pred_perm)
% table=GO_table;
%
% formula='RT~1+pred+Task+(1|SubID)';
% permvars={'Task','SubID'};
% totperm=100;
if nargin<5
    all_pred_perm=[];
end
out_pred_perm=[];
% run real model
if strcmp(predictor,'Drug')
    eval(sprintf('table.pred=table.%s;',predictor));
    model= fitglme(table,formula,'Distribution','binomial');
    real_out=[];
    for nDrug=1:3
        real_out=[real_out ; double(model.Coefficients(2+nDrug,2)) double(model.Coefficients(2+nDrug,4)) double(model.Coefficients(2+nDrug,6)) nDrug];
    end
    
    uniqueIDs=unique(table.SubID);
    uniqueBlocks=unique(table.BlockN);
    perm_out={nan(totperm,4),nan(totperm,4),nan(totperm,4)};
    fprintf('%4.0f/%4.0f\n',0,totperm)
    for np=1:totperm
        if size(all_pred_perm,1)~=totperm
            pred_perm=table.pred;
            for nS=1:length(uniqueIDs)
                idx=find(table.SubID==uniqueIDs(nS));
                temp=pred_perm(idx);
                pred_perm(idx)=temp(randperm(length(temp)));
            end
            out_pred_perm(np,:)=pred_perm;
        else
            pred_perm=all_pred_perm(np,:)';
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
    model= fitglme(table,formula,'Distribution','binomial');
    real_out=[double(model.Coefficients(match_str(model.CoefficientNames,'pred'),2)) double(model.Coefficients(match_str(model.CoefficientNames,'pred'),4)) double(model.Coefficients(match_str(model.CoefficientNames,'pred'),6))];
    
    uniqueIDs=unique(table.SubID);
    uniqueBlocks=unique(table.BlockN);
    perm_out=nan(totperm,4);
    fprintf('%4.0f/%4.0f\n',0,totperm)
    for np=1:totperm
        if size(all_pred_perm,1)~=totperm
            pred_perm=table.pred;
            for nS=1:length(uniqueIDs)
                idx=find(table.SubID==uniqueIDs(nS));
                temp=pred_perm(idx);
                pred_perm(idx)=temp(randperm(length(temp)));
            end
            out_pred_perm(np,:)=pred_perm;
        else
            pred_perm=all_pred_perm(np,:)';
        end
        
        table2=table;
        table2.pred=pred_perm;
        model= fitglme(table2,formula,'Distribution','binomial');
        perm_out(np,:)=[double(model.Coefficients(match_str(model.CoefficientNames,'pred'),2)) double(model.Coefficients(match_str(model.CoefficientNames,'pred'),4)) double(model.Coefficients(match_str(model.CoefficientNames,'pred'),6)) np];
        fprintf('\b\b\b\b\b\b\b\b\b\b%4.0f/%4.0f\n',np,totperm)
    end
    fprintf('\n');
    
end