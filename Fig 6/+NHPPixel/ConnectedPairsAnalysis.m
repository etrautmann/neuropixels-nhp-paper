classdef ConnectedPairsAnalysis < handle
    % use CCGTool to perform cross-correlation analysis, by condition, averaging the results, and finding peaks
    
    properties
        datasetName (1, 1) string
        cluster_ids (:, 1)
        binsMs
        binWidthMs
        fs = 30000;
        
        % hyperparams
        mode string
        smoothWidthMs
        jitter_reps
        thresh_std = 7;
        min_lag_ms = 0;
        max_lag_ms = 10;
        flank_ms = 50;
        window_ms = 100;
        
        % these are unsmoothed versions of the ccgs
        ccg_corrected single % nBins x nClusters x nClusters
        ccg_raw single
        ccg_jittered single
        
        cluster_pair_spike_counts % nClusters x 1
        
        connected (:, :) logical % nClusters x nClusters
        peak_latency (:, :)
        peak_mag (:, :)
    end
    
    properties(Transient)
        ccgtool
        
        % generate post load with smooth_ccgs()
        ccg_corrected_smoothed single % nBins x nClusters x nClusters
        ccg_raw_smoothed single
        ccg_jittered_smoothed single
    end
    
    properties(Dependent)
        nPairsConnected
        nPairsTotal
        connected_pair_ids % nPairs x 2
        connected_pair_latency % nPairs x 1
        connected_pair_mag % nPairs x 1
    end
    
    methods
        function cpa = ConnectedPairsAnalysis(ks_or_seg_by_trials, td, varargin)
            p = inputParser();
            p.addParameter('mode', "", @isstringlike);
            p.addParameter('jitter_reps', 5, @isscalar);
            p.addParameter('cluster_ids', ks_or_seg_by_trials.cluster_ids, @isvector);
            p.addParameter('binWidthMs', 0.5, @isscalar);
            p.addParameter('smoothWidthMs', 1.5, @isscalar);
            p.KeepUnmatched = false;
            p.parse(varargin{:});
            jitter_reps = p.Results.jitter_reps;
            cpa.binWidthMs = p.Results.binWidthMs;
            
            ks_obj = ks_or_seg_by_trials;
            cpa.fs = ks_obj.fsAP;
            
            if isa(ks_obj, 'Neuropixel.KilosortDataset')
                % not using trials, using the whole time period
                ks = ks_obj;
                is_ks = true;
                
            elseif isa(ks_obj, 'Neuropixel.KilosortTrialSegmentedDataset')
                seg_by_trials = ks_obj;
                assert(~seg_by_trials.is_segmented_by_clusters);
                is_ks = false;
                
                assert(~isempty(td));
            end
                
            ccgtool = cpa.build_ccgtool();
            cpa.mode = string(p.Results.mode);
            cpa.binsMs = ccgtool.binsMs;
            cpa.binWidthMs = ccgtool.binWidthMs;
            cpa.cluster_ids = p.Results.cluster_ids;
            cpa.jitter_reps = jitter_reps;
            
            nbins = numel(cpa.binsMs);
            nClusters = numel(cpa.cluster_ids);
            
            if is_ks
                times = ks.spike_times;
                [tf, cluster_inds] = ismember(ks.spike_clusters, cpa.cluster_ids);
%                 assert(all(tf));
                
                % don't smooth here
                debug('Computing raw CCG for entire KS dataset\n');
                [cpa.ccg_corrected, cpa.ccg_raw, cpa.ccg_jittered, cpa.cluster_pair_spike_counts] = ...
                    ccgtool.compute_ccg_normalized_jitter_corrrected(times, cluster_inds, ...
                    'jitter_reps', jitter_reps, 'smoothWidthMs', 0, 'nClusters', nClusters, 'showProgressOverTrials', true);
            else
                % split by trials with the same condition idx and then average
                [ccg_corrected, ccg_raw, ccg_jitter] = deal(zeros(nbins, nClusters, nClusters));
                cc_counts = zeros(nClusters, nClusters);
                
                if td.nConditions == 1
                    warning('TrialDataConditionAlign may be grouped to average condition-wise CCGs');
                end

                trial_durations_samples = seg_by_trials.trial_stop - seg_by_trials.trial_start + ones(1, 'like', seg_by_trials.trial_start);

                unique_conditions = unique(td.conditionIdx(td.valid));
                nUC = numel(unique_conditions);

                nUC_with_trials = 0;
                for iC = 1:nUC
                    debug('Computing jitter-corrected CCG on condition %d / %d\n', iC, nUC);
                    td_mask = td.conditionIdx == unique_conditions(iC);

                    if ~any(td_mask)
                        continue;
                    end
                    
                    nUC_with_trials = nUC_with_trials + 1;
                    times = seg_by_trials.spike_times(td_mask);
                    clusters = seg_by_trials.spike_cluster_inds(td_mask);

                    % don't smooth here
                    [ccg_this, ccg_raw_this, ccg_jitter_this, cc_counts_this] = ccgtool.compute_ccg_normalized_jitter_corrrected(times, clusters, 'jitter_reps', jitter_reps, ...
                        'trial_durations_samples', trial_durations_samples, 'nClusters', seg_by_trials.nClusters, 'smoothWidthMs', 0, 'showProgressOverTrials', true);
                    ccg_corrected = ccg_corrected + ccg_this;
                    ccg_raw = ccg_raw + ccg_raw_this;
                    ccg_jitter = ccg_jitter + ccg_jitter_this;
                    cc_counts = cc_counts + cc_counts_this;
                end

                cpa.ccg_corrected = ccg_corrected ./ nUC_with_trials;
                cpa.ccg_raw = ccg_raw ./ nUC_with_trials;
                cpa.ccg_jittered = ccg_jitter ./ nUC_with_trials;
                cpa.cluster_pair_spike_counts = cc_counts;
            end
            
            % this will set cpa.smoothWidthMs and the ccg_*_smoothed properties
            cpa.smoothWidthMs = p.Results.smoothWidthMs;
            cpa.computeSmoothedCCGs('smoothWidthMs', cpa.smoothWidthMs);
        end
        
        function ccgtool = build_ccgtool(cpa)
            cpa.ccgtool = NHPPixel.CCGTool('fs', cpa.fs, 'windowMs', cpa.window_ms, 'binWidthMs', cpa.binWidthMs);
            ccgtool = cpa.ccgtool;
        end
        
        function computeSmoothedCCGs(cpa, varargin)
            p = inputParser();
            p.addParameter('smoothWidthMs', cpa.smoothWidthMs, @isscalar);
            p.KeepUnmatched = false;
            p.parse(varargin{:});
            assert(p.Results.smoothWidthMs > 0);
            cpa.smoothWidthMs = p.Results.smoothWidthMs;
            
            % smooth with blackman window
            smoothWidthBins = ceil(2*cpa.smoothWidthMs / cpa.binWidthMs);
                
            window = blackman(smoothWidthBins + 2);
            smooth_filter = window(2:end-1);
            smooth_filter = smooth_filter ./ sum(smooth_filter); % must be column vector
            smooth_fn = @(ccg_tensor) convn(ccg_tensor, smooth_filter, 'same');
            
            cpa.ccg_raw_smoothed = smooth_fn(cpa.ccg_raw);
            cpa.ccg_jittered_smoothed = smooth_fn(cpa.ccg_jittered);
            cpa.ccg_corrected_smoothed = cpa.ccg_raw_smoothed - cpa.ccg_jittered_smoothed;
        end
        
        function findConnectedPairs(cpa)
            assert(~isempty(cpa.ccg_corrected_smoothed), 'Call computeSmoothedCCGs() first');
            if isempty(cpa.ccgtool)
                cpa.build_ccgtool();
            end
            [cpa.connected, cpa.peak_latency, cpa.peak_mag] = cpa.ccgtool.find_ccg_peaks(cpa.ccg_corrected_smoothed, ...
                'thresh_std', cpa.thresh_std, 'min_lag_ms', cpa.min_lag_ms, 'max_lag_ms', cpa.max_lag_ms, 'flank_ms', cpa.flank_ms);
        end
        
        function [clusterInds, cluster_ids] = lookup_clusterIds(cpa, cluster_ids)
            if islogical(cluster_ids)
                cluster_ids = cpa.cluster_ids(cluster_ids);
             end
            [tf, clusterInds] = ismember(cluster_ids, cpa.cluster_ids);
            assert(all(tf), 'Some cluster ids were not found in cpa.clusterids');
        end
        
        function nPairs = get.nPairsConnected(cpa)
            nPairs = nnz(cpa.connected);
        end
        
        function nPairs = get.nPairsTotal(cpa)
            nClu = size(cpa.connected, 1);
            nPairs = nClu * (nClu-1);
        end
        
        function pairs = get.connected_pair_ids(cpa)
            [r, c] = find(cpa.connected);
            pairs = [cpa.cluster_ids(r), cpa.cluster_ids(c)];
        end
        
        function latencies = get.connected_pair_latency(cpa)
            latencies = cpa.peak_latency(cpa.connected);
        end
        
        function mags = get.connected_pair_mag(cpa)
            mags = cpa.peak_mag(cpa.connected);
        end
       
        function [ccg, binsMs] = retrieveCCG(cpa, cid1, cid2, varargin)
            p = inputParser();
            p.addParameter('kind', 'corrected', @isstringlike);
            p.parse(varargin{:});
            
            kind = string(p.Results.kind);
            switch kind
                case "corrected"
                    ccgfull = cpa.ccg_corrected;
                case "raw"
                    ccgfull = cpa.ccg_raw;
                case "jittered"
                    ccgfull = cpa.ccg_jittered;
                    
                case "corrected_smoothed"
                    ccgfull = cpa.ccg_corrected_smoothed;
                case "raw_smoothed"
                    ccgfull = cpa.ccg_raw_smoothed;
                case "jittered_smoothed"
                    ccgfull = cpa.ccg_jittered_smoothed;
                    
                otherwise
                    error('Unknown kind %s', kind);
            end
            
            cind1 = cpa.lookup_clusterIds(cid1);
            cind2 = cpa.lookup_clusterIds(cid2);
            ccg = ccgfull(:, cind1, cind2);
            binsMs = cpa.binsMs;
        end
        
        function h = plotCCG(cpa, cid1, cid2, varargin)
            p = inputParser();
            p.addParameter('kind', 'corrected', @isstringlike);
            p.addParameter('showGuides', true, @islogical);
            p.addParameter('style', 'stem', @isstringlike); % stem, line
            p.parse(varargin{:});

            kind = string(p.Results.kind);
            style = string(p.Results.style);
            [ccg, binsMs] = cpa.retrieveCCG(cid1, cid2, 'kind', kind);
            
            switch style
                case "stem"
                    h = stem(binsMs, ccg, '-', 'Marker', 'none', 'LineWidth', 1.5, 'Color', [0.3 0.3 0.3]); %#ok<*PROPLC
                case "line"
                    h = plot(binsMs, ccg, '-', 'Marker', 'none', 'LineWidth', 1, 'Color', [0.3 0.3 0.3]); %#ok<*PROPLC
                case 'stairs'
                    h = stairs(binsMs, ccg, '-', 'Marker', 'none', 'LineWidth', 1, 'Color', [0.3 0.3 0.3]); %#ok<*PROPLC
                case 'area'
                    h = area(binsMs, ccg, 'FaceColor',[.3 .3 .3], 'EdgeColor','none'); %#ok<*PROPLC
                otherwise
                    error('Unknown style %s', style);
            end
            
            if kind == "corrected" && p.Results.showGuides
                % plot thresholding for peak detection
                hold on;
                orange = [255 204 153] / 255;

                cind1 = cpa.lookup_clusterIds(cid1);
                cind2 = cpa.lookup_clusterIds(cid2);

                hx(1) = xline(cpa.min_lag_ms - cpa.binWidthMs/2);
                hx(2) = xline(cpa.max_lag_ms + cpa.binWidthMs/2);
                
                hx(3) = xline(-cpa.flank_ms + cpa.binWidthMs/2);
                hx(4) = xline( cpa.flank_ms - cpa.binWidthMs/2);
                set(hx, 'Color', orange);
                
                mask_bin_flank = abs(binsMs) > cpa.flank_ms & ccg ~= 0;
                if any(mask_bin_flank)
                    flank_std = std(ccg(mask_bin_flank));
                    hy = yline(flank_std * cpa.thresh_std, '-', sprintf('%dx STD', cpa.thresh_std));
                    hy.Color = orange;

                    if cpa.connected(cind1, cind2)
                        plot(cpa.peak_latency(cind1, cind2), cpa.peak_mag(cind1, cind2) * flank_std, 'ro', 'MarkerFaceColor', 'r');
                    end
                end
            end
            
            hold on;
            
            hx = xline(0);
            hx.Color = 'r';
            hx.LineStyle = '--';
            
            hold off;
        end
        
        function h = plotCCG_rawVsJittered(cpa, cid1, cid2, varargin)
            p = inputParser();
            p.addParameter('smoothed', false, @islogical);
            p.addParameter('style', 'stem', @isstringlike); % stem, line
            p.parse(varargin{:});
            style = string(p.Results.style);
            smoothed = p.Results.smoothed;
            
            if smoothed
                [ccg_raw, binsMs] = cpa.retrieveCCG(cid1, cid2, 'kind', 'raw_smoothed');
            else
                [ccg_raw, binsMs] = cpa.retrieveCCG(cid1, cid2, 'kind', 'raw');
            end
            
            switch style
                case "stem"
                    h.raw = stem(binsMs, ccg_raw, '-', 'Marker', 'none', 'LineWidth', 1, 'Color', [0.3 0.3 0.3]); %#ok<*PROPLC
                case "line"
                    h.raw = plot(binsMs, ccg_raw, '-', 'Marker', 'none', 'LineWidth', 1, 'Color', [0.3 0.3 0.3]);
                otherwise
                    error('Unknonw style %s', style);
            end
            hold on;
            
            if smoothed
                [ccg_jitter, binsMs] = cpa.retrieveCCG(cid1, cid2, 'kind', 'jittered_smoothed');
            else
                [ccg_jitter, binsMs] = cpa.retrieveCCG(cid1, cid2, 'kind', 'jittered');
            end
            
            switch style
                case "stem"
                    h.jittered = stem(binsMs, ccg_jitter, '-', 'Marker', 'none', 'LineWidth', 1, 'Color', [1 0.3 0.3]);
                case "line"
                    h.jittered = plot(binsMs, ccg_jitter, '-', 'Marker', 'none', 'LineWidth', 1, 'Color', [1 0.3 0.3]);
            end
            
%             xline(0);
            
%             [ccg, binsMs] = cpa.retrieveCCG(cid1, cid2, 'kind', 'corrected');
%             h.corrected = stem(binsMs, ccg, '-', 'Marker', 'none', 'LineWidth', 1, 'Color', 'k'); %#ok<*PROPLC
            
            hold off;
        end
        
        function h = plotConnectedPairCCG(cpa, pair_ind, varargin)
            pair_ids = cpa.connected_pair_ids(pair_ind, :);
            h = cpa.plotCCG(pair_ids(1), pair_ids(2), varargin{:});
        end 
        
        function h = plotConnectedPairCCG_rawVsJittered(cpa, pair_ind, varargin)
            pair_ids = cpa.connected_pair_ids(pair_ind, :);
            h = cpa.plotCCG_rawVsJittered(pair_ids(1), pair_ids(2), varargin{:});
        end 
    end

end

