function MOL_plotRaster(sessionData,trialData,spikeData,cell_IDs,params)

% params              = params_histresponse();
params.t_pre        = -0.5e6;
params.t_post       = 1e6;
params.zscore       = 0;
% params.conv_win             = 'gaussian';

%%
nNeurons                = length(cell_IDs);

params.nSplits          = 4;

params.AlignOn          = 'stimChange';

for iNeuron = 1:nNeurons %Loop over all neurons:
    
    cell_idx = strcmp(spikeData.cell_ID,cell_IDs(iNeuron));
    
    %Get the relevant data for each session individually:
    temptrialData        = MOL_getTempPerSes(spikeData.session_ID(cell_idx),trialData);
    
    switch params.trialcategories
        case 'ProbeVis'
            splits          = {};
            splits{1}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.vecResponse==3 & temptrialData.visualOriChangeNorm==3;
            splits{2}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.vecResponse==2 & temptrialData.visualOriChangeNorm==3;
            splits{3}       = ismember(temptrialData.trialType,{'P'}) & temptrialData.vecResponse==3;
            splits{4}       = ismember(temptrialData.trialType,{'P'}) & temptrialData.vecResponse==2;
            params.SortBy   = 'Lat';
        case 'ThrMax'
            splits          = {};
            splits{1}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==2 & temptrialData.vecResponse==3;
            splits{2}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==2 & temptrialData.vecResponse==2;
            splits{3}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.vecResponse==3;
            splits{4}       = ismember(temptrialData.trialType,{'X'}) & temptrialData.visualOriChangeNorm==3 & temptrialData.vecResponse==2;
            params.SortBy   = 'Ori';
    end
    
    temptrialData.visualOriPostChangeNorm = round(temptrialData.visualOriPostChangeNorm/2);
    
    if all(cellfun(@sum,splits)>5)
        
        %Get times:
        events_ts               = temptrialData.(params.AlignOn);
        spikes_ts               = spikeData.ts{cell_idx};
        
        %sort spikes
        vecTrialPerSpike        = nan(size(spikes_ts));
        vecTimePerSpike         = nan(size(spikes_ts));
        
        % Get spiking times:
        for ev = 1:length(events_ts) % Get spikes within window around event ev:
            idx                     = spikes_ts>events_ts(ev)+params.t_pre & spikes_ts<=events_ts(ev)+params.t_post;
            vecTrialPerSpike(idx)   = ev;
            vecTimePerSpike(idx)    = spikes_ts(idx)-events_ts(ev);
        end
        
        %Compute histogram:
        events_ts               = temptrialData.(params.AlignOn);
        hist_mat                = calc_psth(events_ts,spikes_ts,params);    %Construct histogram matrix:
        
        splits{4}               = splits{4} & ~(ismember(temptrialData.trialNum,[452 607]) & strcmp(temptrialData.session_ID,'20342019121961'));
        
        for iSplit = 1:params.nSplits %Store the mean response for each of these splits
            ratemat(iSplit,:) = nanmean(hist_mat(splits{iSplit},:),1); %#ok<*AGROW>
        end
        
        figure; set(gcf,'units','normalized','Position',[0.1 0.3 0.3 0.6],'color','w'); hold all;
        subplot(5,1,1); hold all;
        for iSplit = 1:params.nSplits %Store the mean response for each of these splits
            handles(iSplit) = plot(params.xtime,ratemat(iSplit,:),params.lines_ztrials{iSplit},'Color',params.colors_ztrials{iSplit},'LineWidth',3);
        end

        xlabel('')
        xlim([params.t_pre params.t_post]);
        ylim([0 ceil(max(max(ratemat)))])
        set(gca,'XTick',[],'YTick',get(gca,'YLim'))
        plot([0 0],[0 max(max(ratemat))],'k--','LineWidth',2)
        if iNeuron ==1 
            ylabel('Firing rate (sp s-1)','FontSize', 15)
            legend(handles,params.labels_ztrials,'FontSize',10); legend boxoff
        else ylabel('')
        end
        
%         for iVar = 1:length(params.AUC_varselec)
%             temp = find(squeeze(outputmat_sign(params.AUC_varselec(iVar),cell_idx,1:end-1)));
%             for i=1:numel(temp)
%                 plot([params.xtime(temp(i)) params.xtime(temp(i)+1)],[ymax+iVar ymax+iVar],'Color',params.colors_AUC{params.AUC_varselec(iVar)});
%             end
%         end
        
        
        for i = 1:params.nSplits
            subplot(5,1,i+1); hold all;

            %Sort trials based on post-change orientation and response latency (orientation is sorted by adding a lot of latency)
            trialselec          = find(splits{i}');
            
            switch params.SortBy
                case 'Ori'
                    temp                = temptrialData.responseLatency(trialselec);
                    temp(isnan(temp))   = 0;
                    [~,idx]             = sort(temp + 1e13*temptrialData.visualOriPostChangeNorm(trialselec),'descend');
                case 'Lat'
                    temp                = temptrialData.responseLatency(trialselec);
                    temp(isnan(temp))   = 0;
                    [~,idx]             = sort(temp + rand(size(temp)),'descend');
            end
            trialselec          = trialselec(idx);
            
            %plot per trial
%             h = patch([params.t_pre params.t_post params.t_post params.t_pre],[0 0 length(trialselec) length(trialselec)],'k');
%             h.FaceAlpha = 0.1;
%             h.FaceAlpha = 1;
%             h.EdgeColor = [1 1 1];
            
            for intTrial = 1:length(trialselec)%find(trialselec') %1:numel(vecTrialStarts)
                
                vecTimes = vecTimePerSpike(vecTrialPerSpike==trialselec(intTrial));
                lat = temptrialData.responseLatency(trialselec(intTrial));
                if ~isnan(lat)
                    line([lat lat],[intTrial-0.5;intTrial+0.5],'Color','k','LineWidth',2);
                end
                line([vecTimes(:)';vecTimes(:)'],[intTrial*ones(1,numel(vecTimes))-1;intTrial*ones(1,numel(vecTimes))+1],'Color',params.colors_ztrials{i},'LineWidth',1);
%                 line([vecTimes(:)';vecTimes(:)'],[intTrial*ones(1,numel(vecTimes))-0.5;intTrial*ones(1,numel(vecTimes))+0.5],'Color',params.colors_ztrials{i},'LineWidth',2);
            end
            plot([0 0],[0 length(trialselec)],'k--','LineWidth',2)
           
            %set fig props
            ylim([0 length(trialselec)]);
            xlim([params.t_pre params.t_post]);
            if i==4
                xlabel('Time from stimulus change (s)');
                set(gca,'XTick',[-0.5 0 0.5 1]*1e6,'XTickLabel',[-0.5 0 0.5 1],'YTick',[])
            else
                xlabel('')
                set(gca,'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
            end
            ylabel(params.labels_ztrials{i},'FontSize',15);
            
        end
        explabel = sessionData.Experiment(strcmp(sessionData.session_ID,spikeData.session_ID(cell_idx)));
        title([cell_IDs(iNeuron) explabel])
        if isfield(params,'export_fig') && params.export_fig
            export_fig(fullfile(params.savedir,sprintf('ExNeuron%d_%s',iNeuron,spikeData.session_ID{cell_idx})),'-eps','-nocrop')
        end
    end
    
end



end