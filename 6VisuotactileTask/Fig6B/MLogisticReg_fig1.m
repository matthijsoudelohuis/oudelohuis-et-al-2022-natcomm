%% Multinomial Logistic Regression - Two regressors
% To be used on trialData preprocessed with PieTaskBeh_fittingfig1
% JPie 22-10-2020
%% Input

%Select predictors
    SelPred = 'Visual'; %!! 'Visual' or 'Tactile' for now
%Select trials to fit:
    if strcmp(SelPred, 'Visual')
    	SelConds = [1 2 3 4 5 6 7];
    elseif strcmp(SelPred, 'Tactile')
    	SelConds = [1 8 9 10 11];
    end
    

%% Data preparation    
iscond = ismember(trialData.conditionnumber,SelConds);
tfields = fieldnames(trialData);
for k=1:length(tfields) %clean trials that didnt happen
        if (isnumeric( trialData.(tfields{k}) ) || iscell( trialData.(tfields{k}) ) || islogical(trialData.(tfields{k}) ))
            trialData.(tfields{k}) = trialData.(tfields{k})(iscond);
        end
end
%
if strcmp(SelPred, 'Visual')
    predictor = 100*trialData.stimContrast;
elseif strcmp(SelPred, 'Tactile')
    predictor = trialData.stimDeflection;
end
%
stimR = predictor; stimR(stimR<0)=0;
stimL = predictor; stimL(stimL>0)=0; stimL = -stimL;
%
Response = 3*ones(size(trialData.firstRespLFR)); %Reference response is "no lick" as number 3.
Response(trialData.firstRespLFR==76) = 1; %lick Left
Response(trialData.firstRespLFR==82) = 2; %lick Right
sp = categorical(Response);
%Build Y (used in likelihood calc)
[stimCombs, ~, iS] = unique([stimL stimR],'rows');
totalTrials_percomb =  accumarray(iS,1); %! not in "plotting order"
Y = nan(size(stimCombs,1),3);
for k = 1:size(stimCombs,1)
    Y(k,:) = [ sum(Response(iS==k)==1) sum(Response(iS==k)==2) sum(Response(iS==k)==3)];
end

%% Finding best mle fit for n
N = [0.1:0.01:1];
LL = nan(1,length(N)); %loglikelihood

for it = 1:length(N)
    n = N(it);
    X = stimCombs.^n;
    [B] = mnrfit(X,Y); 
    X2 = [ones(size(X,1),1) X];
    Y2 = Y(:,1:2);
    LL(it) = sum(sum(Y2.*(X2*B))) - sum( totalTrials_percomb.*log(1+(  sum(exp(X2*B),2)  )) );
end
%Get best n
[~,bestni] = max(LL);

%% Do regression for best n and compute probabilities
n = N(bestni);
X = stimCombs.^n;

xL = [0: 0.1 : 100^(n)+0.1];
xR = [0: 0.1 : 100^(n)+0.1];
[B,dev,stats] = mnrfit(X,Y); 
pL=[];
pR=[];

for i = 1:length(xL)
    for j = 1:length(xR)
        RHS_L(i,j) = B(1,1) + B(2,1)*xL(i) + B(3,1)*xR(j); 
        RHS_R(i,j) = B(1,2) + B(2,2)*xL(i) + B(3,2)*xR(j); 
        pL(i,j) = exp(RHS_L(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
        pR(i,j) = exp(RHS_R(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
    end
end

%% Plot

figure()
hold on
plot_stimCombs = sum([stimCombs(:,1)*-1 stimCombs(:,2)],2); %!! 
R = Y ./ repmat(totalTrials_percomb,1,3);
plot(plot_stimCombs,R(:,1),'s','Color',[0 0.8 0],'MarkerSize',5)
plot(plot_stimCombs,R(:,2),'s','Color',[0.8 0 0],'MarkerSize',5)
% plot(100*Contrasts,HitRight(:,(Deflections==0)),'s','Color',[0.8 0 0],'MarkerSize',5)
% plot(100*Contrasts,HitLeft(:,(Deflections==0)),'s','Color',[0 0.8 0],'MarkerSize',5)
ylim([0 1])
plot(xR.^(1/n),pR(1,:),'r') %
plot(-xL.^(1/n),pR(:,1),'r')
plot(-xL.^(1/n),pL(:,1),'g')
plot(xR.^(1/n),pL(1,:),'g')
% plot(xR,pR(1,:),'r') %
% plot(-xL,pR(:,1),'r')
% plot(-xL,pL(:,1),'g')
% plot(xR,pL(1,:),'g')
% plot(xR,pN(:,1),'Color',[0.6 0.6 0.6])
% plot(-xL,pN(1,:),'Color',[0.6 0.6 0.6]) 
xlabel('Contrast')
% xlabel('Contrast')
ylabel('Choice rate')
legend({'Right','Left'})

%% Save
%Also for fast plotting: xR,xL,pR,pL (model) trC, HR, HL. (data)
%save(strcat('MLogisticReg_',SelPred,'_',sessionData.Mouse(1:4)),'n','stimCombs','Y','B','trialData','sessionData','MergedStr','xR','xL','pR','pL');