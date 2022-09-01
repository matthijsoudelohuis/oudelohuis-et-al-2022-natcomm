function MOL_plotOptoBehavior_Dprime(params,dVis_UST,dAud_UST,dVis_MST,dAud_MST,Mice_UST,Mice_MST)
%% Dprime figure:
offset = 4;

% Make figure:
figure; set(gcf,'color','w','units','normalized','Position', [0.4 0.3 .28 .376]); hold all;

xpos = [1 2];
fprintf('Visual thr change UST, Ctrl vs Early:\n')
MOL_subPlot_Dprime_Condition(xpos,params.colors_experiments{2},[1 1 1],squeeze(dVis_UST(:,1,[1 2])),Mice_UST);

xpos = [3 4];
fprintf('Visual max change UST, Ctrl vs Early:\n')
MOL_subPlot_Dprime_Condition(xpos,params.colors_experiments{2},params.colors_experiments{2},squeeze(dVis_UST(:,2,[1 2])),Mice_UST);

xpos = [1 2] - 0.2;
fprintf('Visual thr change MST, Ctrl vs Early:\n')
MOL_subPlot_Dprime_Condition(xpos,params.colors_experiments{3},[1 1 1],squeeze(dVis_MST(:,1,[1 2])),Mice_MST);

xpos = [3 4] - 0.2;
fprintf('Visual max change MST, Ctrl vs Early::\n')
MOL_subPlot_Dprime_Condition(xpos,params.colors_experiments{3},params.colors_experiments{3},squeeze(dVis_MST(:,2,[1 2])),Mice_MST);

xpos = [1 2] + 1 * offset;
fprintf('Visual thr change UST, Ctrl vs Late:\n')
MOL_subPlot_Dprime_Condition(xpos,params.colors_experiments{2},[1 1 1],squeeze(dVis_UST(:,1,[1 3])),Mice_UST);

xpos = [3 4] + 1 * offset;
fprintf('Visual max change UST, Ctrl vs Late:\n')
MOL_subPlot_Dprime_Condition(xpos,params.colors_experiments{2},params.colors_experiments{2},squeeze(dVis_UST(:,2,[1 3])),Mice_UST);

xpos = [1 2]  + 1 * offset - 0.2;
fprintf('Visual thr change MST, Ctrl vs Late:\n')
MOL_subPlot_Dprime_Condition(xpos,params.colors_experiments{3},[1 1 1],squeeze(dVis_MST(:,1,[1 3])),Mice_MST);

xpos = [3 4] + 1 * offset - 0.2;
fprintf('Visual max change MST, Ctrl vs Late:\n')
MOL_subPlot_Dprime_Condition(xpos,params.colors_experiments{3},params.colors_experiments{3},squeeze(dVis_MST(:,2,[1 3])),Mice_MST);

%Make up:
ylabel('Dprime')
XTickLabels = repmat({'Ctrl' 'Early' 'Ctrl' 'Late'},1,2);
set(gca,'XTick',1:length(XTickLabels),'XTickLabels',XTickLabels,'XTickLabelRotation',60);
set(gca,'YTick',[0 1 2])
ylim([-0.3 2.2])
xlim([0.4 length(XTickLabels)+0.2])
grid on 

%Auditory figure;
figure; set(gcf,'color','w','units','normalized','Position', [0.4 0.3 .28 .376]); hold all;
xpos = [1 2];
fprintf('Audio thr change MST, Ctrl vs Early:\n')
MOL_subPlot_Dprime_Condition(xpos,params.colors_experiments{3},[1 1 1],squeeze(dAud_MST(:,1,[1 2])),Mice_MST);
xpos = [3 4];
fprintf('Audio thr change MST, Ctrl vs Late:\n')
MOL_subPlot_Dprime_Condition(xpos,params.colors_experiments{3},[1 1 1],squeeze(dAud_MST(:,1,[1 3])),Mice_MST);

xpos = [1 2];
fprintf('Audio max change MST, Ctrl vs Early:\n')
MOL_subPlot_Dprime_Condition(xpos+4,params.colors_experiments{3},params.colors_experiments{3},squeeze(dAud_MST(:,2,[1 2])),Mice_MST);
xpos = [3 4];
fprintf('Audio max change MST, Ctrl vs Late:\n')
MOL_subPlot_Dprime_Condition(xpos+4,params.colors_experiments{3},params.colors_experiments{3},squeeze(dAud_MST(:,2,[1 3])),Mice_MST);

%Make up:
ylabel('Dprime')
XTickLabels = repmat({'Ctrl' 'Early' 'Ctrl' 'Late'},1,2);
set(gca,'XTick',1:length(XTickLabels),'XTickLabels',XTickLabels,'XTickLabelRotation',60);
set(gca,'YTick',[0 1 2])
ylim([-0.3 2.6])
xlim([0.4 length(XTickLabels)+0.2])
grid on 

end


function MOL_subPlot_Dprime_Condition(xpos,coloredge,colorfill,dMat,Mice)

y_mean      = squeeze(nanmedian(dMat(:,:))); y_mean = y_mean(~isnan(y_mean));
y_std_lower = nanmedian(dMat(:,:)) - prctile(dMat(:,:),25); % / sqrt(sum(~isnan(dMat(:,1)))); y_std = y_std(~isnan(y_std));
y_std_upper = prctile(dMat(:,:),75) - nanmedian(dMat(:,:)); %/ sqrt(sum(~isnan(dMat(:,1)))); y_std = y_std(~isnan(y_std));
errorbar(xpos, y_mean,y_std_lower,y_std_upper,'-','Color','k','MarkerSize',0.001,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',3);
errorbar(xpos, y_mean,[0 0],[0 0],'o','MarkerSize',12,'MarkerEdgeColor',coloredge,'MarkerFaceColor',colorfill,'LineWidth',3);

% errorbar(xpos, y_mean,[0 0],[0 0],'-','Color','k','MarkerSize',0.001,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',5);
% errorbar(xpos, y_mean,y_std_lower,y_std_upper,'o','Color',color,'MarkerSize',12,'MarkerEdgeColor',color,'MarkerFaceColor',color,'LineWidth',0.0001);

%Statistical testing:    
Y               = [dMat(:,1); dMat(:,2)];
X_inac          = [ones(size(dMat,1),1); ones(size(dMat,1),1)*2];
G_mou           = [Mice; Mice;];
tbl             = table(Y,X_inac,G_mou,'VariableNames',{'Dprime','Inac','Mouse'}); %Create table for mixed model

lme             = fitlme(tbl,'Dprime~Inac+(1|Mouse)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
stats{2,5}      = stats{2,5}*4; %Bonferroni Holm multiple comparison correction
fprintf('F(%d,%2.0f)=%1.2f, p=%1.2f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
p = stats{2,5};
if p<0.05
    sigstar(xpos([1 2]),p) %use sigstar function to identify significance 
end

fprintf('n=%d sessions, p=%1.6f\n',size(dMat,1),p)


% %Statistical testing:    signed rank test with bonferroni correction:
% [p] = signrank(dMat(~isnan(dMat(:,2)),1),dMat(~isnan(dMat(:,2)),2)); 
% p = min([p*3 1]); %bonferroni correction
% if p<0.05
%     sigstar(xpos([1 2]),p) %use sigstar function to identify significance 
% end
% fprintf('n=%d sessions, p=%1.6f\n',size(dMat,1),p)

end
