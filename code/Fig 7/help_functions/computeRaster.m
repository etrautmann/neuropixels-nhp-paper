function [spikeRaster,b_TriggerTimes_valid] = computeRaster(TriggerTimes,spikeTimecourses,RasterTimeBeforeMS,RasterTimeAfterMS)
% spikeTimecourses needs to be n_cells x n_time_points
if ~iscell(spikeTimecourses)
    if size(spikeTimecourses,2)==1
        spikeTimecourses = spikeTimecourses';
    end
    spikeTimecourses = num2cell(spikeTimecourses,2);
end
b_TriggerTimes_valid = and(TriggerTimes-RasterTimeBeforeMS>0,TriggerTimes+RasterTimeAfterMS<=length(spikeTimecourses{1}));
TriggerTimes = TriggerTimes(b_TriggerTimes_valid); %Exclude Trigger times that are out of range.

if all(cellfun(@islogical,spikeTimecourses))
    spikeRaster = false(length(TriggerTimes),length(-RasterTimeBeforeMS:RasterTimeAfterMS),length(spikeTimecourses));
else
    spikeRaster = nan(length(TriggerTimes),length(-RasterTimeBeforeMS:RasterTimeAfterMS),length(spikeTimecourses));
end
TriggerTimesIntervals = zeros(length(TriggerTimes),length((-RasterTimeBeforeMS):RasterTimeAfterMS));
for TriggerTimesIndex=1:length(TriggerTimes)
    TriggerTimesIntervals(TriggerTimesIndex,:) = (TriggerTimes(TriggerTimesIndex)-RasterTimeBeforeMS):(TriggerTimes(TriggerTimesIndex)+RasterTimeAfterMS);
end
for cellIndex = 1:length(spikeTimecourses)
    spikeRaster(:,:,cellIndex) = spikeTimecourses{cellIndex}(round(TriggerTimesIntervals));
end
end
