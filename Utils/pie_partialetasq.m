function[pe] = pie_partialetasq(lme)
tbl = lme.Variables;
%Compute effect size pe (partial eta squared) for all fixed effects predictors.
predicted = tbl.(lme.ResponseName);
[beta, betanames] = fixedEffects(lme);
bn = table2cell(betanames);
%With categorical predictors, the name is appended, remove and check:
vals = cell(1,size(betanames,1));
for i = 2:size(betanames,1)
    if ~isempty(findstr(bn{i},'_'))
        bnt = bn{i}(1:(findstr(bn{i},'_')-1));
        valt = bn{i}((findstr(bn{i},'_')+1):end);
        if sum(strcmp(lme.PredictorNames, bnt))>0
            bn{i} = bnt;
            vals{i} = valt;
        end
    end
end
all_predictors = [ones(size(tbl,1),1)];
pe = nan(1,(size(betanames,1)-1));
for iPred = 2:size(betanames,1)
    if iscategorical(tbl.(bn{iPred}))
    all_predictors = [all_predictors double(tbl.(bn{iPred})==vals{iPred})];
    elseif iscell(tbl.(bn{iPred}))
    all_predictors = [all_predictors double(strcmp(tbl.(bn{iPred}),vals{iPred}))];    
    else
    all_predictors = [all_predictors tbl.(bn{iPred})];
    end
end
SSR = nanvar( predicted - all_predictors*beta);
for iPred = 2:size(betanames,1)
    predictor_test = [all_predictors(:,1) all_predictors(:,iPred)];
    SSEp = nanvar(predictor_test*[beta(1);beta(iPred)]);
    pe(iPred) = SSEp / (SSEp + SSR);
end
end
