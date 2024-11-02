classdef CCGTool < handle
    properties(SetAccess=protected)
        windowMs
        binWidthMs
        fs % samples per second
        binsMs % bin centers
    end
    
    properties(Dependent)
        windowSamples
        binWidthSamples
        
        nbins
        binMinMs
        binMaxMs
        ind_bin_zero
    end
    
    methods
        function tool = CCGTool(varargin)
            p = inputParser();
            p.addParameter('fs', 30000, @isscalar);
            p.addParameter('windowMs', 100, @isscalar);
            p.addParameter('binWidthMs', 1, @isscalar);
            p.parse(varargin{:});
            
            tool.fs = p.Results.fs;
            tool.windowMs = p.Results.windowMs;
            tool.binWidthMs = p.Results.binWidthMs;
            
            tool.windowMs = ceil(tool.windowMs / tool.binWidthMs) * tool.binWidthMs;
            tool.binsMs = (-tool.windowMs : tool.binWidthMs : tool.windowMs)';
        end
        
        function nbins = get.nbins(tool)
            nbins = numel(tool.binsMs);
        end
            
        function bin_min = get.binMinMs(tool)
            half_bin = tool.binWidthMs/2;
            bin_min = tool.binsMs - half_bin;
        end
        
        function bin_max = get.binMaxMs(tool)
            half_bin = tool.binWidthMs/2;
            bin_max = tool.binsMs + half_bin;
        end
        
        function ind_bin_zero = get.ind_bin_zero(tool)
            [~, ind_bin_zero] = min(abs(tool.binsMs));
        end
        
        function windowSamples = get.windowSamples(tool)
            windowSamples =  tool.windowMs * tool.fs / 1000;
        end
        
        function binWidthSamples = get.binWidthSamples(tool)
            binWidthSamples =  tool.binWidthMs * tool.fs / 1000;
        end
        
        function [ccg_tensor, pair_counts] = compute_ccg_normalized(tool, times_cell, clusters_cell, varargin)
            % compute single trial raw CCG, normalized by triangle function T - |tau|
            % that compensates for edges where some bins are undercounted
            % c.f. https://www.jneurosci.org/content/25/14/3661#sec-2
            %
            % times_cell and clusters_cell are both nTrials x 1 cells with sorted spike times / cluster inds within them as vectors
            %
            % ccg_tensor(tbin, i1, i2) counts (ignoring normalization) spikes of i2 that occured at tbin latency relative to i1 as the reference spike
            % normalized by geometric mean of firing rates, averaged over trials
            
            p = inputParser();
            p.addParameter('trial_durations_samples', [], @(x) isempty(x) || isvector(x));
            p.addParameter('nClusters', NaN, @isscalar);
            p.addParameter('smoothWidthMs', NaN, @isscalar);
            p.addParameter('showProgress', false, @islogical);
            p.addParameter('progressMessage', 'Computing multi-CCG over trials', @ischar);
            
            p.addParameter('accumulateEachTrial', false, @islogical); % true is potentially slower but requires far less memory
            p.addParameter('auto_split_samples', 1000000, @isscalar); % split very long single trial or full dataset into multiple shorter trials to run more efficiently
            p.parse(varargin{:});
            accumulateEachTrial = p.Results.accumulateEachTrial;
            
            trial_durations_samples = double(p.Results.trial_durations_samples);
            if isempty(trial_durations_samples)
                trial_durations_samples = double(cellfun(@max, times_cell));
            end
            nClusters = p.Results.nClusters;
            if isnan(nClusters)
                nClusters = max(cellfun(@max, clusters_cell));
            end
            smoothWidthMs = p.Results.smoothWidthMs;
            if isnan(smoothWidthMs)
                smoothWidthMs = 0;
            end
            
            if ~iscell(times_cell)
                times_cell = {times_cell};
                clusters_cell = {clusters_cell};
            end
            nTrials = numel(times_cell);
            
            auto_split_samples = p.Results.auto_split_samples;
            if nTrials == 1 && ~isinf(auto_split_samples)
                % rather than deal with a single enormous trial, we automatically split it into multiple trials
                [min_time, max_time] = bounds(times_cell{1});
                bins = min_time : uint64(auto_split_samples) : max_time + uint64(auto_split_samples);
                n_per_bin = histcounts(times_cell{1}, bins);
                
                times_cell = mat2cell(times_cell{1}, n_per_bin, 1);
                clusters_cell = mat2cell(clusters_cell{1}, n_per_bin, 1);
                
                trial_durations_samples = repmat(auto_split_samples, numel(times_cell), 1);
                trial_durations_samples(end) = max_time - bins(end-1) + uint64(1);
                nTrials = numel(times_cell);
                
                accumulateEachTrial = true; % trials will be very large so this is a good idea
            end
            
            showProgress = p.Results.showProgress && nTrials > 1;
            
            % in samples
            nbins = tool.nbins;
            counts = zeros(nClusters, 1); % accumulate the total per-cluster spike counts for rate normalization, must be column
            bin_opportunities = zeros(nbins, 1); % accumulate the number of opportunities to see a given lag based on trial durations
            
            % we loop through the spikes once with two counters. ihigh is the leading spike, ilow is the 
            if showProgress
                prog = ProgressBar(nTrials, p.Results.progressMessage);
%                 prog.enableParallel();
            end
            
            nbinsfromzero = (tool.nbins - 1) / 2;
            binWidth_samples = tool.binWidthSamples;
            dt_samples = (tool.binsMs(end) + tool.binWidthMs/2) * tool.fs / 1000; % max offset in samples to include in ccg
            ind_bin_zero = tool.ind_bin_zero; %#ok<*PROPLC>
            
            bin_indices = (-nbinsfromzero : nbinsfromzero)'; % must be column for division below to work
            assert(numel(bin_indices) == nbins);
            
            if accumulateEachTrial
                ccg_counts = zeros([tool.nbins, nClusters, nClusters], 'uint32');
            else
                subs = cell(nTrials, 1);
            end
%             parfor (iTrial = 1:nTrials, 8)
            for iTrial = 1:nTrials
                times = times_cell{iTrial};
                clusters = clusters_cell{iTrial};
                trial_dur_nbins = trial_durations_samples(iTrial) / tool.binWidthSamples;
                bin_opportunities = bin_opportunities + (trial_dur_nbins - abs(bin_indices));
                counts = counts + accumarray(clusters, 1, [nClusters, 1]);
                
                % use range search for all times against all other times
                idx = rangesearch(single(times), single(times), dt_samples);
                [i2, i1] = TensorUtils.catWhich(2, idx{:});
                i2 = i2';
                
                ibin = ind_bin_zero + round((double(times(i2)) - double(times(i1)))/binWidth_samples);
                c1 = clusters(i1);
                c2 = clusters(i2);
                % range search occasionally returns times with offsets == dt_samples which would land in the immediately
                % flanking bins, but anything further out than that is problematic
                if any(ibin < 0 | ibin > nbins+1)
                    error('Issue with bins');
                end
                mask_keep = i1 ~= i2 & ibin >= 1 & ibin <= nbins;
               
                % build contribution to ccg_tensor directly using accumarray
                subs_this_trial = [ibin(mask_keep), c1(mask_keep), c2(mask_keep)];
                
                if accumulateEachTrial
                    ccg_counts = ccg_counts + uint32(accumarray(subs_this_trial, 1, [tool.nbins, nClusters, nClusters]));
                else
                    subs{iTrial}  = subs_this_trial;
                end
                if showProgress, prog.update(iTrial); end
            end
            if showProgress, prog.finish(); end
            
            if ~accumulateEachTrial
                ccg_counts = accumarray(cat(1, subs{:}), 1, [tool.nbins, nClusters, nClusters]);
                clear subs
            end
                
            rates = counts ./ sum(trial_durations_samples); % nClusters x 1
            geom_mean_rates = shiftdim(sqrt(rates .* rates'), -1); % 1 x nClusters x nClusters
            geom_mean_rates(geom_mean_rates == 0) = 1;
            
            pair_counts = shiftdim(sum(ccg_counts, 1), 1);
            ccg_tensor = single(ccg_counts) ./ bin_opportunities ./ geom_mean_rates;
            
            % smooth with blackman window
            smoothWidthBins = ceil(2*smoothWidthMs / tool.binWidthMs);
                
            if smoothWidthBins > 1
                window = blackman(smoothWidthBins + 2);
                smooth_filter = window(2:end-1);
                smooth_filter = smooth_filter ./ sum(smooth_filter); % must be column vector

                ccg_tensor = convn(ccg_tensor, smooth_filter, 'same');
            end
        end
        
        function times_cell = interval_jitter_times(tool, times_cell, varargin)
            % must not modify the order of times within each cell since clusters will not be correspondingly updated!
            
            p = inputParser();
            p.addParameter('interval_width_ms', 25, @isscalar);
            p.parse(varargin{:});
            
            interval_width_ms = p.Results.interval_width_ms;
            if isnan(interval_width_ms)
                interval_width_ms = 25;
            end
            interval_width_samples = interval_width_ms * tool.fs / 1000;
            
%             trial_durations_samples = double(p.Results.trial_durations_samples);
%             if isempty(trial_durations_samples)
%                 trial_durations_samples = cellfun(@max, times_cell);
%             end
            
            if ~iscell(times_cell)
                times_cell = {times_cell};
            end
            nTrials = numel(times_cell);
            
            for iTrial = 1:nTrials
                times = times_cell{iTrial};
                max_time = double(max(times));
                nIntervals = ceil(max_time / interval_width_samples);
                width_last_interval = rem(max_time, interval_width_samples);
                
                which_interval = floor(times / interval_width_samples);
                orig_offset_within_interval = rem(times, interval_width_samples);      
                new_offset_within_interval = rand(size(times));
                new_offset_within_interval(which_interval < nIntervals) = new_offset_within_interval(which_interval < nIntervals) * interval_width_samples;
                new_offset_within_interval(which_interval == nIntervals) = new_offset_within_interval(which_interval == nIntervals) * width_last_interval;
                    
                updated_times = times - orig_offset_within_interval + cast(new_offset_within_interval, 'like', times);
                times_cell{iTrial} = updated_times;
            end   
        end
        
        function [ccg_tensor, ccg_raw, ccg_jitter, pair_counts] = compute_ccg_normalized_jitter_corrrected(tool, times_cell, clusters_cell, varargin)
            % computes an interval-jitter corrected CCG in which ccg_jitter is computed by resampling spikes within unifor
            % intervals of interval_width_ms 
            p = inputParser();
            p.addParameter('interval_width_ms', NaN, @isscalar);
            p.addParameter('trial_durations_samples', [], @isvector);
            p.addParameter('jitter_reps', 1, @isscalar);
            p.addParameter('nClusters', NaN, @isscalar);
            p.addParameter('smoothWidthMs', NaN, @isscalar);
            p.addParameter('showProgressOverTrials', false, @islogical);
            p.parse(varargin{:});
            
            jitter_reps = p.Results.jitter_reps;
            show_progress_over_trials = p.Results.showProgressOverTrials;
            
            if ~iscell(times_cell)
                times_cell = {times_cell};
                clusters_cell = {clusters_cell};
            end

            if ~show_progress_over_trials
                debug('Computing raw CCG\n');
            end
            compute_ccg_fn = @(times_cell, varargin) tool.compute_ccg_normalized(times_cell, clusters_cell, 'trial_durations_samples', p.Results.trial_durations_samples, ...
                'nClusters', p.Results.nClusters, 'smoothWidthMs', p.Results.smoothWidthMs, 'showProgress', show_progress_over_trials, varargin{:});
            [ccg_raw, pair_counts] = compute_ccg_fn(times_cell);
            
            rng(0);
            ccg_jitter = zeros(size(ccg_raw));

            if ~show_progress_over_trials
                prog = ProgressBar(jitter_reps, 'Computing jittered CCGs\n');
            end
            for i = 1:jitter_reps
                jitter_times_cell = tool.interval_jitter_times(times_cell, 'interval_width_ms', p.Results.interval_width_ms); 
                ccg_jitter = ccg_jitter + compute_ccg_fn(jitter_times_cell, 'showProgress', show_progress_over_trials, ...
                    'progressMessage', 'Computing jittered multi-CCG over trials');
                if ~show_progress_over_trials
                    prog.update(i);
                end
            end
            if ~show_progress_over_trials
                prog.finish();
            end
            
            ccg_jitter = ccg_jitter ./ jitter_reps;
                
            ccg_tensor = ccg_raw - ccg_jitter;    
        end
        
        function [has_peak, peak_latency, peak_mag] = find_ccg_peaks(tool, ccg_tensor, varargin)
            % searches the ccg_tensor to find pairs with a short latency peak. 
            % this requires that the peak of the CCG (positive or negative) exceeds the std of the flank values (outside +/- flank_ms)
            % by thresh_std times (e.g. 7x) and occurs in the interval [min_lag_ms, max_lag_ms]. 
            % % ccg_tensor(tbin, i1, i2) counts (ignoring normalization) spikes of i2 that occured at tbin latency relative to i1 as the reference spike
            % if the peak occurs at 0, then we indicate the peak for i1 <= i2 only.
            %
            % has_peak(i1, i2) indicates that i1 spikes influence the probability of i2 spikes in the future (peak_latency(i1,i2)
            % peak_mag is the value of the ccg at that latency offset normalized by the flank_std.
            % if we postulate this as a mono-synaptic connection, then has_peak(i1, i2) suggests i1 projects to i2, such that has_peak(i1, :) indicates it's post-synaptic partners
            p = inputParser();
            p.addParameter('thresh_std', NaN, @isscalar)
            p.addParameter('min_lag_ms', NaN, @isscalar);
            p.addParameter('max_lag_ms', NaN, @isscalar);
            p.addParameter('flank_ms', NaN, @isscalar);
            p.addParameter('min_pair_count', NaN, @isscalar);
            p.parse(varargin{:});
            
            thresh_std = p.Results.thresh_std;
            if isnan(thresh_std)
                thresh_std = 7;
            end
            max_lag_ms = p.Results.max_lag_ms;
            if isnan(max_lag_ms)
                max_lag_ms = 10;
            end
            min_lag_ms = p.Results.min_lag_ms;
            if isnan(min_lag_ms)
                min_lag_ms = 0;
            end
            assert(min_lag_ms >= 0);
            
            flank_ms = p.Results.flank_ms;
            if isnan(flank_ms)
                flank_ms = 50;
            end
            
            bins = tool.binsMs;
            mask_peak_valid = bins >= min_lag_ms & bins <= max_lag_ms;
%             ind_search_offset = find(mask_bin_search, 1) - 1;
            mask_bin_flank = abs(bins) > flank_ms;
            
            if ~any(mask_bin_flank)
                error('tool.windowMs=%f is not large enough to include flank_ms=%f', tool.windowMs, flank_ms);
            end
            
            N = size(ccg_tensor, 2);
            has_peak = false(N, N);
            peak_latency = nan(N, N);
            peak_mag = nan(N, N);
            for i1 = 1:N
                for i2 = 1:N
                    if i1 == i2
                        continue;
                    end     
                    this_ccg = ccg_tensor(:, i1, i2);
                    
                    % determine std of flanking points
                    mask_bin_flank_this = mask_bin_flank & this_ccg ~= 0;
                    this_flank_std = std(this_ccg(mask_bin_flank_this));
                    if this_flank_std == 0
                        peak_mag(i1, i2) = NaN;
                        continue;
                    end
                    
                    % find absolute largest ccg point
                    [~, ind] = max(abs(this_ccg));
                    val = this_ccg(ind);
                    val_norm_std = val / this_flank_std;

                    % the peak must exceed the threshold and lie in the valid range of latencies
                    if abs(val_norm_std) >= thresh_std && mask_peak_valid(ind)
                        has_peak(i1, i2) = true;
                    end
                    
                    % regardless of whether peak is eligible, note the largest value regardless
                    peak_latency(i1, i2) = bins(ind);
                    peak_mag(i1, i2) = val_norm_std;
                    
                    % in the event that the peak eligibility window overlaps 0, we may have both pairs
                    % detect the same peak for i1,i2 and i2,i1. 
                    % Ensure we only detect it for i1 <= i2 (upper triangle)
                    if i1 > i2 && has_peak(i2, i1)
                        has_peak(i1, i2) = false;
                    end
                    
                end
            end      
        end
        
        function [ccg_tensor, cluster_counts] = compute_ccg_batchwise(tool, times, cluster_inds, batch_inds, varargin)
            % compute raw CCG, normalized by triangle function T - |tau|
            % that compensates for edges where some bins are undercounted
            % c.f. https://www.jneurosci.org/content/25/14/3661#sec-2
            %
            % times, clusters, batches are all nSpikes x 1 with sorted spike times / cluster inds within them as vectors
            % batches(i) gives the batch number of the raw spike
            %
            % ccg_tensor is nbins x nBatches x nClusters x nClusters

            p = inputParser();
            p.addParameter('nClusters', NaN, @isscalar);
            p.addParameter('showProgress', true, @islogical);
            p.addParameter('progressMessage', 'Computing multi-CCG over batches', @ischar);
            p.parse(varargin{:});
            
            nClusters = p.Results.nClusters;
            if isnan(nClusters)
                nClusters = max(cluster_inds);
            end
            nBatches = max(batch_inds);
            showProgress = p.Results.showProgress && nBatches > 1;
            
            % in samples
            nbins = tool.nbins;
            cluster_counts = zeros(nClusters, nBatches); % accumulate the total per-cluster spike counts for rate normalization, must be column
            
            % we loop through the spikes once with two counters. ihigh is the leading spike, ilow is the 
            if showProgress
                prog = ProgressBar(nBatches, p.Results.progressMessage);
            end
            
            nbinsfromzero = (tool.nbins - 1) / 2;
            binWidth_samples = tool.binWidthSamples;
            dt_samples = (tool.binsMs(end) + tool.binWidthMs/2) * tool.fs / 1000; % max offset in samples to include in ccg
            ind_bin_zero = tool.ind_bin_zero; %#ok<*PROPLC>
            
            bin_indices = single(-nbinsfromzero : nbinsfromzero)'; % must be column for division below to work
            assert(numel(bin_indices) == nbins);
            
            %ccg_tensor = zeros(tool.nbins, nBatches, nClusters, nClusters, 'single');
            
            % segment by batch
            batchwise_times = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(times, 1, nBatches, batch_inds);
            batchwise_cluster_inds = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(cluster_inds, 1, nBatches, batch_inds);
            nbins = tool.nbins;
            ccg_tensor_batches = cell(nBatches, 1);
           for iBatch = 1:nBatches
                times_this = batchwise_times{iBatch};
                clusters_this = batchwise_cluster_inds{iBatch};
                batch_dur_nbins = (max(times_this) - min(times_this) + ones(1,1, 'like', times_this)) / binWidth_samples;
                bin_opportunities = single(batch_dur_nbins) - abs(bin_indices);
                cluster_counts(:, iBatch) = accumarray(clusters_this, 1, [nClusters, 1]);
                
                % use range search for all times against all other times
                idx = rangesearch(single(times_this), single(times_this), single(dt_samples));
                [i2, i1] = TensorUtils.catWhich(2, idx{:});
                i2 = i2';
                
                ibin = ind_bin_zero + round((double(times_this(i2)) - double(times_this(i1)))/binWidth_samples);
                c1 = clusters_this(i1);
                c2 = clusters_this(i2);
                % range search occasionally returns times with offsets == dt_samples which would land in the immediately
                % flanking 1 or 2 bins, but anything further out than that is problematic
                if any(ibin < -1 | ibin > nbins+2)
                    disp(iBatch)
                    error('Issue with bins');
                end
                mask_keep = i1 ~= i2 & ibin >= 1 & ibin <= nbins;
               
                % build contribution to ccg_tensor directly using accumarray
                subs = [ibin(mask_keep), c1(mask_keep), c2(mask_keep)];
                unnormalized_counts = accumarray(subs, 1, [nbins, nClusters, nClusters]);
                ccg_tensor_batches{iBatch} = reshape(single(unnormalized_counts), [nbins, 1, nClusters, nClusters]) ./ bin_opportunities;
                if showProgress, prog.update(iBatch); end %#ok<PFBNS>
            end
            if showProgress, prog.finish(); end
            
            ccg_tensor = cat(2, ccg_tensor_batches{:});
        end
           
    end
    
end