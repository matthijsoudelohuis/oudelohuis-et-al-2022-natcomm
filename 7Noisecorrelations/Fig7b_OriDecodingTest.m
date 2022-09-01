%Load exported orientation decoding scores from python: 
datadir         = 'E:\SurfDrive\Documents\Manuscript - 2nd Bump\Submission - NatComm\rev2\OudeLohuisetal_2022_NatComms';
filename        = fullfile(datadir,'7Noisecorrelations','SourceData_7b_OriDecoding_PrePost.csv');
tbl             = readtable(filename);

%Fit LMM:
lme             = fitlme(tbl,'score~group+(1|animal_id)'); %construct linear mixed effects model with fixed effect of temporal window and random intercept for different mice
stats           = dataset2table(anova(lme,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3e; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
p = stats{2,5};


%% Figure; 
scores = table2array(tbl(:,3));
prescores = scores(strcmp(table2array(tbl(:,2)),{'pre'}));
postscores = scores(strcmp(table2array(tbl(:,2)),{'post'}));

nSessions = length(prescores);
cohorts     = table2array(tbl(:,6));

plotcohorts = cohorts(strcmp(table2array(tbl(:,2)),{'pre'}));

params = MOL_getColors_CHDET();

% Figure of dominance dependence on last trial:
figure; hold all; set(gcf,'units','normalized','Position',[0.2 0.3 0.07 0.2],'color','w');
for iSes = 1:nSessions
    
    plot([1 2],[prescores(iSes) postscores(iSes)],'k-','LineWidth',1); %scatter the individual animals
    switch plotcohorts{iSes}
        case 'ChangeDetectionConflictDecor'
            clr = params.colors_experiments{1};
        case 'VisOnlyTwolevels'
            clr = params.colors_experiments{2};
        case 'ChangeDetectionConflict'
            clr = params.colors_experiments{3};
    end
    plot([1 2],[prescores(iSes) postscores(iSes)],'.','MarkerSize',20,'Color',clr); %scatter the individual animals
    
end
set(gca,'XTick',[1 2],'XTickLabel',{'Pre' 'Post'},'XTickLabelRotation',45,'YAxisLocation','right') %xtick
sigstar([1 2],p)
xlim([0.7 2.3])
ylim([-0.05 0.4])
% export_fig(fullfile(params.savedir,filename),gcf)


