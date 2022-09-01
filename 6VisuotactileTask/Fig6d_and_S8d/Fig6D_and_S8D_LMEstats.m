% This code is an upgrade of the statistics done for Fig6d and S8d
% It applies Linear Mixed Effects statistical models to the silencing data
% BehLat in Data6_3

%% Prep
sessionnames = categorical({BehLat.SessionName});
sessionnum = double(sessionnames);
for i=1:length(BehLat)
    SN = BehLat(i).SessionName;
    BehLat(i).Mouse = SN(strfind(SN,'J'):strfind(SN,'J')+3);
end

%% Low saliency visual contra
%Subselect MST Early
tbl = struct2table( BehLat(strcmp({BehLat.Side},'Contra') & strcmp({BehLat.Modality},'Visual') & strcmp({BehLat.Task},'tMulti') & strcmp({BehLat.Silencing},'Early') & strcmp({BehLat.Saliency},'Low')) );
tbl2 = tbl; tbl2.Dprime = tbl2.DprimeCtrl; tbl2.Silencing = repmat({'OFF'},size(tbl,1),1); tbl2.LaserOnset = repmat([1],size(tbl,1),1);
tbl = [tbl;tbl2]; 
clear tbl2
tbl.Silencing   = categorical(tbl.Silencing);
tbl.Mouse       = categorical(tbl.Mouse);
lme1 = fitlme(tbl,'Dprime ~ Silencing + (1|Mouse)');

%Subselect MST Late
tbl = struct2table( BehLat(strcmp({BehLat.Side},'Contra') & strcmp({BehLat.Modality},'Visual') & strcmp({BehLat.Task},'tMulti') & strcmp({BehLat.Silencing},'Late') & strcmp({BehLat.Saliency},'Low')) );
tbl2 = tbl; tbl2.Dprime = tbl2.DprimeCtrl; tbl2.Silencing = repmat({'OFF'},size(tbl,1),1); tbl2.LaserOnset = repmat([1],size(tbl,1),1);
tbl = [tbl;tbl2]; 
clear tbl2
tbl.Silencing   = categorical(tbl.Silencing);
tbl.Mouse       = categorical(tbl.Mouse);
lme2 = fitlme(tbl,'Dprime ~ Silencing + (1|Mouse)');

%Subselect UST Early
tbl = struct2table( BehLat(strcmp({BehLat.Side},'Contra') & strcmp({BehLat.Modality},'Visual') & strcmp({BehLat.Task},'tUni') & strcmp({BehLat.Silencing},'Early') & strcmp({BehLat.Saliency},'Low')) );
tbl2 = tbl; tbl2.Dprime = tbl2.DprimeCtrl; tbl2.Silencing = repmat({'OFF'},size(tbl,1),1); tbl2.LaserOnset = repmat([1],size(tbl,1),1);
tbl = [tbl;tbl2]; 
clear tbl2
tbl.Silencing   = categorical(tbl.Silencing);
tbl.Mouse       = categorical(tbl.Mouse);
lme3 = fitlme(tbl,'Dprime ~ Silencing + (1|Mouse)');

%Subselect UST Late
tbl = struct2table( BehLat(strcmp({BehLat.Side},'Contra') & strcmp({BehLat.Modality},'Visual') & strcmp({BehLat.Task},'tUni') & strcmp({BehLat.Silencing},'Late') & strcmp({BehLat.Saliency},'Low')) );
tbl2 = tbl; tbl2.Dprime = tbl2.DprimeCtrl; tbl2.Silencing = repmat({'OFF'},size(tbl,1),1); tbl2.LaserOnset = repmat([1],size(tbl,1),1);
tbl = [tbl;tbl2]; 
clear tbl2
tbl.Silencing   = categorical(tbl.Silencing);
tbl.Mouse       = categorical(tbl.Mouse);
lme4 = fitlme(tbl,'Dprime ~ Silencing + (1|Mouse)');


%% High saliency visual contra
%Subselect MST Early
tbl = struct2table( BehLat(strcmp({BehLat.Side},'Contra') & strcmp({BehLat.Modality},'Visual') & strcmp({BehLat.Task},'tMulti') & strcmp({BehLat.Silencing},'Early') & strcmp({BehLat.Saliency},'High')) );
tbl2 = tbl; tbl2.Dprime = tbl2.DprimeCtrl; tbl2.Silencing = repmat({'OFF'},size(tbl,1),1); tbl2.LaserOnset = repmat([1],size(tbl,1),1);
tbl = [tbl;tbl2]; 
clear tbl2
tbl.Silencing   = categorical(tbl.Silencing);
tbl.Mouse       = categorical(tbl.Mouse);
lme5 = fitlme(tbl,'Dprime ~ Silencing + (1|Mouse)');

%Subselect MST Late
tbl = struct2table( BehLat(strcmp({BehLat.Side},'Contra') & strcmp({BehLat.Modality},'Visual') & strcmp({BehLat.Task},'tMulti') & strcmp({BehLat.Silencing},'Late') & strcmp({BehLat.Saliency},'High')) );
tbl2 = tbl; tbl2.Dprime = tbl2.DprimeCtrl; tbl2.Silencing = repmat({'OFF'},size(tbl,1),1); tbl2.LaserOnset = repmat([1],size(tbl,1),1);
tbl = [tbl;tbl2]; 
clear tbl2
tbl.Silencing   = categorical(tbl.Silencing);
tbl.Mouse       = categorical(tbl.Mouse);
lme6 = fitlme(tbl,'Dprime ~ Silencing + (1|Mouse)');

%Subselect UST Early
tbl = struct2table( BehLat(strcmp({BehLat.Side},'Contra') & strcmp({BehLat.Modality},'Visual') & strcmp({BehLat.Task},'tUni') & strcmp({BehLat.Silencing},'Early') & strcmp({BehLat.Saliency},'High')) );
tbl2 = tbl; tbl2.Dprime = tbl2.DprimeCtrl; tbl2.Silencing = repmat({'OFF'},size(tbl,1),1); tbl2.LaserOnset = repmat([1],size(tbl,1),1);
tbl = [tbl;tbl2]; 
clear tbl2
tbl.Silencing   = categorical(tbl.Silencing);
tbl.Mouse       = categorical(tbl.Mouse);
lme7 = fitlme(tbl,'Dprime ~ Silencing + (1|Mouse)');

%Subselect UST Late
tbl = struct2table( BehLat(strcmp({BehLat.Side},'Contra') & strcmp({BehLat.Modality},'Visual') & strcmp({BehLat.Task},'tUni') & strcmp({BehLat.Silencing},'Late') & strcmp({BehLat.Saliency},'High')) );
tbl2 = tbl; tbl2.Dprime = tbl2.DprimeCtrl; tbl2.Silencing = repmat({'OFF'},size(tbl,1),1); tbl2.LaserOnset = repmat([1],size(tbl,1),1);
tbl = [tbl;tbl2]; 
clear tbl2
tbl.Silencing   = categorical(tbl.Silencing);
tbl.Mouse       = categorical(tbl.Mouse);
lme8 = fitlme(tbl,'Dprime ~ Silencing + (1|Mouse)');

%% Low-saliency ipsi 
%Subselect MST Early
tbl = struct2table( BehLat(strcmp({BehLat.Side},'Ipsi') & strcmp({BehLat.Modality},'Visual') & strcmp({BehLat.Task},'tMulti') & strcmp({BehLat.Silencing},'Early') & strcmp({BehLat.Saliency},'Low')) );
tbl2 = tbl; tbl2.Dprime = tbl2.DprimeCtrl; tbl2.Silencing = repmat({'OFF'},size(tbl,1),1); tbl2.LaserOnset = repmat([1],size(tbl,1),1);
tbl = [tbl;tbl2]; 
clear tbl2
tbl.Silencing   = categorical(tbl.Silencing);
tbl.Mouse       = categorical(tbl.Mouse);
lme9 = fitlme(tbl,'Dprime ~ Silencing + (1|Mouse)');

%Subselect MST Late
tbl = struct2table( BehLat(strcmp({BehLat.Side},'Ipsi') & strcmp({BehLat.Modality},'Visual') & strcmp({BehLat.Task},'tMulti') & strcmp({BehLat.Silencing},'Late') & strcmp({BehLat.Saliency},'Low')) );
tbl2 = tbl; tbl2.Dprime = tbl2.DprimeCtrl; tbl2.Silencing = repmat({'OFF'},size(tbl,1),1); tbl2.LaserOnset = repmat([1],size(tbl,1),1);
tbl = [tbl;tbl2]; 
clear tbl2
tbl.Silencing   = categorical(tbl.Silencing);
tbl.Mouse       = categorical(tbl.Mouse);
lme10 = fitlme(tbl,'Dprime ~ Silencing + (1|Mouse)');

%Subselect UST Early
tbl = struct2table( BehLat(strcmp({BehLat.Side},'Ipsi') & strcmp({BehLat.Modality},'Visual') & strcmp({BehLat.Task},'tUni') & strcmp({BehLat.Silencing},'Early') & strcmp({BehLat.Saliency},'Low')) );
tbl2 = tbl; tbl2.Dprime = tbl2.DprimeCtrl; tbl2.Silencing = repmat({'OFF'},size(tbl,1),1); tbl2.LaserOnset = repmat([1],size(tbl,1),1);
tbl = [tbl;tbl2]; 
clear tbl2
tbl.Silencing   = categorical(tbl.Silencing);
tbl.Mouse       = categorical(tbl.Mouse);
lme11 = fitlme(tbl,'Dprime ~ Silencing + (1|Mouse)');

%Subselect UST Late
tbl = struct2table( BehLat(strcmp({BehLat.Side},'Ipsi') & strcmp({BehLat.Modality},'Visual') & strcmp({BehLat.Task},'tUni') & strcmp({BehLat.Silencing},'Late') & strcmp({BehLat.Saliency},'Low')) );
tbl2 = tbl; tbl2.Dprime = tbl2.DprimeCtrl; tbl2.Silencing = repmat({'OFF'},size(tbl,1),1); tbl2.LaserOnset = repmat([1],size(tbl,1),1);
tbl = [tbl;tbl2]; 
clear tbl2
tbl.Silencing   = categorical(tbl.Silencing);
tbl.Mouse       = categorical(tbl.Mouse);
lme12 = fitlme(tbl,'Dprime ~ Silencing + (1|Mouse)');

%% High-saliency ipsi was not silenced

%% Obtain stats

%Obtain stats - Vis Conta Low ---------
stats = dataset2table(anova(lme1,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Multivariate Linear Mixed Model (Vis Contra Low - MST - Early Inac:)\n')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
p_me1 = stats{2,5};

stats = dataset2table(anova(lme2,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Multivariate Linear Mixed Model (Vis Contra Low - MST - Late Inac:)\n')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
p_ml1 = stats{2,5};

stats = dataset2table(anova(lme3,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Multivariate Linear Mixed Model (Vis Contra Low - UST - Early Inac:)\n')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
p_ue1 = stats{2,5};

stats = dataset2table(anova(lme4,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Multivariate Linear Mixed Model (Vis Contra Low - UST - Late Inac:)\n')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
p_ul1 = stats{2,5};

% Vis Contra High ---------
stats = dataset2table(anova(lme5,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Multivariate Linear Mixed Model (Vis Contra High - MST - Early Inac:)\n')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
p_me2 = stats{2,5};

stats = dataset2table(anova(lme6,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Multivariate Linear Mixed Model (Vis Contra High - MST - Late Inac:)\n')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
p_ml2 = stats{2,5};

stats = dataset2table(anova(lme7,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Multivariate Linear Mixed Model (Vis Contra High - UST - Early Inac:)\n')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
p_ue2 = stats{2,5};

stats = dataset2table(anova(lme8,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Multivariate Linear Mixed Model (Vis Contra High - UST - Late Inac:)\n')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
p_ul2 = stats{2,5};

% Vis Ipsi Low ---------
stats = dataset2table(anova(lme9,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Multivariate Linear Mixed Model (Vis Contra High - MST - Early Inac:)\n')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
p_me3 = stats{2,5};

stats = dataset2table(anova(lme10,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Multivariate Linear Mixed Model (Vis Contra High - MST - Late Inac:)\n')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
p_ml3 = stats{2,5};

stats = dataset2table(anova(lme11,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Multivariate Linear Mixed Model (Vis Contra High - UST - Early Inac:)\n')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
p_ue3 = stats{2,5};

stats = dataset2table(anova(lme12,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Multivariate Linear Mixed Model (Vis Contra High - UST - Late Inac:)\n')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
p_ul3 = stats{2,5};


%% Stats summary without correction
pvalues = [p_me1 p_ml1 p_ue1 p_ul1 p_me2 p_ml2 p_ue2 p_ul2 p_me3 p_ml3 p_ue3 p_ul3]
pvalues<0.05

%% Stats Bonferroni-Holm correction
% All together (12 comps)
% corrected_pvalues = bonf_holm(pvalues);
% corrected_pvalues(corrected_pvalues>1) = 1
% corrected_pvalues<0.05

% Grouped by trial type (4 comps per)
% Group 1 : Vis Contra Low
pvalues1 = [p_me1 p_ml1 p_ue1 p_ul1];
corrected_pv1 = bonf_holm(pvalues1); corrected_pv1(corrected_pv1>1) = 1;
% Group 2 : Vis Contra High
pvalues2 = [p_me2 p_ml2 p_ue2 p_ul2];
corrected_pv2 = bonf_holm(pvalues2); corrected_pv2(corrected_pv2>1) = 1;
% Group 3 : Vis Ipsi Low
pvalues3 = [p_me3 p_ml3 p_ue3 p_ul3];
corrected_pv3 = bonf_holm(pvalues3); corrected_pv3(corrected_pv3>1) = 1;

%% Figure S8D: Low saliency tactile contra (High-saliency tactile was not silenced)
%Subselect MST Early
tbl = struct2table( BehLat(strcmp({BehLat.Side},'Contra') & strcmp({BehLat.Modality},'Tactile') & strcmp({BehLat.Task},'tMulti') & strcmp({BehLat.Silencing},'Early') & strcmp({BehLat.Saliency},'Low')) );
tbl2 = tbl; tbl2.Dprime = tbl2.DprimeCtrl; tbl2.Silencing = repmat({'OFF'},size(tbl,1),1); tbl2.LaserOnset = repmat([1],size(tbl,1),1);
tbl = [tbl;tbl2]; 
clear tbl2
tbl.Silencing   = categorical(tbl.Silencing);
tbl.Mouse       = categorical(tbl.Mouse);
lme13 = fitlme(tbl,'Dprime ~ Silencing + (1|Mouse)');
stats = dataset2table(anova(lme13,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Multivariate Linear Mixed Model (Tac Contra Low - MST - Early Inac:)\n')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
p_mtace = stats{2,5};

%Subselect MST Late
tbl = struct2table( BehLat(strcmp({BehLat.Side},'Contra') & strcmp({BehLat.Modality},'Tactile') & strcmp({BehLat.Task},'tMulti') & strcmp({BehLat.Silencing},'Late') & strcmp({BehLat.Saliency},'Low')) );
tbl2 = tbl; tbl2.Dprime = tbl2.DprimeCtrl; tbl2.Silencing = repmat({'OFF'},size(tbl,1),1); tbl2.LaserOnset = repmat([1],size(tbl,1),1);
tbl = [tbl;tbl2]; 
clear tbl2
tbl.Silencing   = categorical(tbl.Silencing);
tbl.Mouse       = categorical(tbl.Mouse);
lme14 = fitlme(tbl,'Dprime ~ Silencing + (1|Mouse)');
stats = dataset2table(anova(lme14,'DFMethod','Satterthwaite')); %Perform ANOVA on model and output as matrix
fprintf('Multivariate Linear Mixed Model (Tac Contra Low - MST - Late Inac:)\n')
fprintf('F(%d,%2.0f)=%1.2f, p=%1.3f; \n',stats{2,3},stats{2,4},stats{2,2},stats{2,5})
p_mtacl = stats{2,5};

%Bonferroni correction
pvalues4 = [p_mtace p_mtacl];
corrected_pv4 = bonf_holm(pvalues4); corrected_pv4(corrected_pv4>1) = 1;

