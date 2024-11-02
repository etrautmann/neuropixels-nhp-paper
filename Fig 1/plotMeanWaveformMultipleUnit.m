sdfunction [] = plotMeanWaveformMultipleUnit(ks, varargin)
%PLOTMEANWAVEFORMMULTIPLEUNIT Summary of this function goes here
%   Detailed explanation goes here
p = inputParser();
p.addParameter('figh',[], @ishandle);
p.addParameter('plotClusters',[],@isnumeric);
p.addParameter('nChan',30,@isscalar);
p.addParameter('yGain',.3,@isscalar);
p.addParameter('spikeWindow',[-15 25],@isnumeric);
p.addParameter('showChannelLabels',false,@islogical);
p.addParameter('chanRange',[],@isnumeric);
p.parse(varargin{:});


nNeurons = length(p.Results.plotClusters);
cmap = (colorcet('R4','N',nNeurons));

if isempty(p.Results.figh)
    figh = figure();
else
    figh = p.Results.figh;
end

set(figh,'color','w');
set(0, 'CurrentFigure', figh)


chanIds = [];
allaxes = [];
for ii = 1:nNeurons
    axh = subtightplot(1,nNeurons, ii, [.1 .02]);

    % Note: can do this all in one call, but dividing it into two allows
    % for calculating mean waveforms from a larger num_waveforms without
    % cluttering the plot with all of those waveforms. 

    % plot individual waveforms.

    ss1 = ks.getWaveformsFromRawData('cluster_id', p.Results.plotClusters(ii),...
        'num_waveforms', 10, ...
        'best_n_channels', p.Results.nChan,  ...
        'subtractOtherClusters', false, ...
        'window', p.Results.spikeWindow);

    maskChan = find(ismember(ss1.channel_ids_by_snippet(:,1),p.Results.chanRange));

    [hdata, settings, scale_info] = ss1.plotAtProbeLocations('showIndividual', true, ...
                                    'showMean',false, ...
                                    'maskChannels', maskChan, ...
                                    'showChannelLabels', false, ...
                                    'plotArgs',{'color',[.75*[1,1,1], 1] ,'linewidth',.4}, ...
                                    'xmag', 4, ...
                                    'axh', axh, ...
                                    'gain',p.Results.yGain, ...
                                    'colormap',cmap(ii,:));
    
    chanIds(:,ii) = ss1.channel_ids_by_snippet(:,1);

    % Plot mean waveforms
    ss2 = ks.getWaveformsFromRawData('cluster_id', p.Results.plotClusters(ii), ...
        'num_waveforms', 100, ...
        'best_n_channels', p.Results.nChan, ...
        'subtractOtherClusters', false, ...
        'window', p.Results.spikeWindow);
    
    hold on
    out2 = ss2.plotAtProbeLocations('showIndividual', false, ...
                                    'showMedian',true, ...
                                    'maskChannels', maskChan, ...
                                    'showChannelLabels', p.Results.showChannelLabels, ...
                                    'plotArgs',{'color',.4*[1 1 1],'linewidth',3}, ...
                                    'xmag', 4, ...
                                    'axh', axh, ...
                                    'gain',p.Results.yGain, ...
                                    'colormap',cmap(ii,:));
    lines = out2.waveformMedian;  
    if ii == 1
        ylimits = ylim;
    end

    set(axh,'Color',.85*[1 1 1])
    set(axh, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[],'xcolor','none','ycolor','none')
    title(['unit ' num2str(p.Results.plotClusters(ii))])
    axis on

    allaxes(ii) = axh;
end
linkaxes(allaxes,'y')

end

