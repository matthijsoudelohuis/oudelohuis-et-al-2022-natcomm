% With right and left side plots
%% UST 
folder = strcat(maindirectory,'\OudeLohuisetal_2022_NatComms\6VisuotactileTask\Data6_1\UST');
%Visual
files = dir(strcat(folder,'\**\All\MLogisticReg_Visual*.mat'));
%fR
fR = figure(); %Figure right
hold on
%scaling factor
scf = 4;
xlimmin = 0.05;

for iMuis = 1:length(files)
    load(fullfile(files(iMuis).folder,files(iMuis).name)) 
    %Need to scale to 0.01 till maximum
    xR = [0:0.01:0.99 1:0.1:25.1].^n;
    xL = xR;
    %Recalculate pR and pL for new x-axis
    pL=[];pR=[];RHS_L=[];RHS_R=[];
    for i = 1%:length(xL)
        for j = 1:length(xR)
            RHS_L(i,j) = B(1,1) + B(2,1)*xL(i) + B(3,1)*xR(j); 
            RHS_R(i,j) = B(1,2) + B(2,2)*xL(i) + B(3,2)*xR(j); 
            pL(i,j) = exp(RHS_L(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
            pR(i,j) = exp(RHS_R(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
        end
    end    
    xxR = scf*xR.^(1/n); 
    plot(xxR,100*pR(1,:),'r','LineWidth',1) 
    plot(xxR,100*pL(1,:),'g','LineWidth',1)
    %Plot catch
    plot(xlimmin, 100*pR(1,1), 'or')
    plot(xlimmin, 100*pL(1,1), 'og')
    
end
ylim([0 100]); 
xlabel('Contrast')
ylabel('Response Rate')
set(gca,'XScale','log')
% xlim([0.01 25]);
xlim([xlimmin 100]);
yticks([0 25 50 75 100])
%%
%fL
fL = figure(); %Figure left
hold on
for iMuis = 1:length(files)
    load(fullfile(files(iMuis).folder,files(iMuis).name)) 
    %Need to scale to 0.01 till maximum
    xR = [0:0.01:0.99 1:0.1:25.1].^n;
    xL = xR;
    %Recalculate pR and pL for new x-axis
    pL=[];pR=[];RHS_L=[];RHS_R=[];
    for i = 1:length(xL)
        for j = 1:length(xR)
            RHS_L(i,j) = B(1,1) + B(2,1)*xL(i) + B(3,1)*xR(j); 
            RHS_R(i,j) = B(1,2) + B(2,2)*xL(i) + B(3,2)*xR(j); 
            pL(i,j) = exp(RHS_L(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
            pR(i,j) = exp(RHS_R(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
        end
    end    
    xxL = scf*xL.^(1/n); 
    plot(xxL,100*pR(:,1),'r','LineWidth',1) 
    plot(xxL,100*pL(:,1),'g','LineWidth',1)
    
    %Plot catch
    plot(xlimmin, 100*pR(1,1), 'or')
    plot(xlimmin, 100*pL(1,1), 'og')
    
end
ylim([0 100]); 
xlabel('Contrast')
ylabel('Response Rate')
set(gca,'XScale','log')
set(gca,'XDir','reverse')
xlim([xlimmin 100]);
yticks([0 25 50 75 100])
%%
%Tactile
folder = strcat(maindirectory,'\OudeLohuisetal_2022_NatComms\6VisuotactileTask\Data6_1\UST');
%fR
figure(); %Figure right
hold on
for iMuis = 1:length(files)
    load(fullfile(files(iMuis).folder,files(iMuis).name)) 
    xR = [0:0.01:0.99 1:0.1:9.9 10:1:100].^n; xL=xR;
    %Recalculate pR and pL for new x-axis
    pL=[];pR=[];RHS_L=[];RHS_R=[];
    for i = 1:length(xL)
        for j = 1:length(xR)
            RHS_L(i,j) = B(1,1) + B(2,1)*xL(i) + B(3,1)*xR(j); 
            RHS_R(i,j) = B(1,2) + B(2,2)*xL(i) + B(3,2)*xR(j); 
            pL(i,j) = exp(RHS_L(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
            pR(i,j) = exp(RHS_R(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
        end
    end 
    xxR = xR.^(1/n);
    plot(xxR,100*pR(1,:),'Color',[.6 0 0],'LineWidth',1) 
    plot(xxR,100*pL(1,:),'Color',[0 .6 0],'LineWidth',1)
end
ylim([0 100]); 
xlabel('Deflection')
ylabel('Response Rate')
set(gca,'XScale','log')
xlim([0.05 100]);
yticks([0 25 50 75 100])
%xticks([0 25 50 75 100])
%fL
figure(); %Figure left
hold on
for iMuis = 1:length(files)
    load(fullfile(files(iMuis).folder,files(iMuis).name)) 
    xR = [0:0.01:0.99 1:0.1:9.9 10:1:100].^n; xL=xR;
    %Recalculate pR and pL for new x-axis
    pL=[];pR=[];RHS_L=[];RHS_R=[];
    for i = 1:length(xL)
        for j = 1:length(xR)
            RHS_L(i,j) = B(1,1) + B(2,1)*xL(i) + B(3,1)*xR(j); 
            RHS_R(i,j) = B(1,2) + B(2,2)*xL(i) + B(3,2)*xR(j); 
            pL(i,j) = exp(RHS_L(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
            pR(i,j) = exp(RHS_R(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
        end
    end 
    xxL = xL.^(1/n);
    plot(xxL,100*pR(:,1),'Color',[.6 0 0],'LineWidth',1) 
    plot(xxL,100*pL(:,1),'Color',[0 .6 0],'LineWidth',1)
end
ylim([0 100]); 
xlabel('Deflection')
ylabel('Response Rate')
set(gca,'XScale','log')
xlim([0.05 100]);
yticks([0 25 50 75 100])
%xticks([0 25 50 75 100])
set(gca,'XDir','reverse')

%% MST
folder = strcat(maindirectory,'\OudeLohuisetal_2022_NatComms\6VisuotactileTask\Data6_1\MST');
files = dir(strcat(folder,'\**\all\MLogisticReg_Visual*.mat'));
%fR
fR = figure(); %Figure right
hold on
for iMuis = 1:length(files)
    load(fullfile(files(iMuis).folder,files(iMuis).name)) 
    %Need to scale to 0.01 till maximum
    xR = [0:0.01:0.99 1:0.1:25.1].^n;
    xL = xR;
    %Recalculate pR and pL for new x-axis
    pL=[];pR=[];RHS_L=[];RHS_R=[];
    for i = 1:length(xL)
        for j = 1:length(xR)
            RHS_L(i,j) = B(1,1) + B(2,1)*xL(i) + B(3,1)*xR(j); 
            RHS_R(i,j) = B(1,2) + B(2,2)*xL(i) + B(3,2)*xR(j); 
            pL(i,j) = exp(RHS_L(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
            pR(i,j) = exp(RHS_R(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
        end
    end    
    xxR = 4*xR.^(1/n); 
    plot(xxR,100*pR(1,:),'r','LineWidth',1) 
    plot(xxR,100*pL(1,:),'g','LineWidth',1)
end
ylim([0 100]); 
xlabel('Contrast')
ylabel('Response Rate')
set(gca,'XScale','log')
% xlim([0.01 25]);
xlim([0.05 100]);
yticks([0 25 50 75 100])
%fL
fL = figure(); %Figure left
hold on
for iMuis = 1:length(files)
    load(fullfile(files(iMuis).folder,files(iMuis).name)) 
    %Need to scale to 0.01 till maximum
    xR = [0:0.01:0.99 1:0.1:25.1].^n;
    xL = xR;
    %Recalculate pR and pL for new x-axis
    pL=[];pR=[];RHS_L=[];RHS_R=[];
    for i = 1:length(xL)
        for j = 1:length(xR)
            RHS_L(i,j) = B(1,1) + B(2,1)*xL(i) + B(3,1)*xR(j); 
            RHS_R(i,j) = B(1,2) + B(2,2)*xL(i) + B(3,2)*xR(j); 
            pL(i,j) = exp(RHS_L(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
            pR(i,j) = exp(RHS_R(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
        end
    end    
    xxL = 4*xL.^(1/n); 
    plot(xxL,100*pR(:,1),'r','LineWidth',1) 
    plot(xxL,100*pL(:,1),'g','LineWidth',1)
end
ylim([0 100]); 
xlabel('Contrast')
ylabel('Response Rate')
set(gca,'XScale','log')
set(gca,'XDir','reverse')
xlim([0.05 100]);
yticks([0 25 50 75 100])

%Tactile
files = dir(strcat(folder,'\**\all\MLogisticReg_Tactile*.mat'));
%fR
figure(); %Figure right
hold on
for iMuis = 1:length(files)
    load(fullfile(files(iMuis).folder,files(iMuis).name)) 
    xR = [0:0.01:0.99 1:0.1:9.9 10:1:100].^n; xL=xR;
    %Recalculate pR and pL for new x-axis
    pL=[];pR=[];RHS_L=[];RHS_R=[];
    for i = 1:length(xL)
        for j = 1:length(xR)
            RHS_L(i,j) = B(1,1) + B(2,1)*xL(i) + B(3,1)*xR(j); 
            RHS_R(i,j) = B(1,2) + B(2,2)*xL(i) + B(3,2)*xR(j); 
            pL(i,j) = exp(RHS_L(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
            pR(i,j) = exp(RHS_R(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
        end
    end 
    xxR = xR.^(1/n);
    plot(xxR,100*pR(1,:),'Color',[.6 0 0],'LineWidth',1) 
    plot(xxR,100*pL(1,:),'Color',[0 .6 0],'LineWidth',1)
end
ylim([0 100]); 
xlabel('Deflection')
ylabel('Response Rate')
set(gca,'XScale','log')
xlim([0.05 100]);
yticks([0 25 50 75 100])
%xticks([0 25 50 75 100])
%fL
figure(); %Figure left
hold on
for iMuis = 1:length(files)
    load(fullfile(files(iMuis).folder,files(iMuis).name)) 
    xR = [0:0.01:0.99 1:0.1:9.9 10:1:100].^n; xL=xR;
    %Recalculate pR and pL for new x-axis
    pL=[];pR=[];RHS_L=[];RHS_R=[];
    for i = 1:length(xL)
        for j = 1:length(xR)
            RHS_L(i,j) = B(1,1) + B(2,1)*xL(i) + B(3,1)*xR(j); 
            RHS_R(i,j) = B(1,2) + B(2,2)*xL(i) + B(3,2)*xR(j); 
            pL(i,j) = exp(RHS_L(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
            pR(i,j) = exp(RHS_R(i,j)) ./ (1+ exp(RHS_L(i,j)) + exp(RHS_R(i,j)));
        end
    end 
    xxL = xL.^(1/n);
    plot(xxL,100*pR(:,1),'Color',[.6 0 0],'LineWidth',1) 
    plot(xxL,100*pL(:,1),'Color',[0 .6 0],'LineWidth',1)
end
ylim([0 100]); 
xlabel('Deflection')
ylabel('Response Rate')
set(gca,'XScale','log')
xlim([0.05 100]);
yticks([0 25 50 75 100])
%xticks([0 25 50 75 100])
set(gca,'XDir','reverse')