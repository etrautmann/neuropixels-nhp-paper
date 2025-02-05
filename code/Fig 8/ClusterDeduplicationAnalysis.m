classdef ClusterDeduplicationAnalysis < handle
    properties(Transient)
        ks 
    end
    
    properties
        cluster_ids (:, 1) uint32
        
        % ccg_frac_spikes(c1, c2) now indicates the fraction of c1 spikes co-occurring with a reliably timed c2 spike
        % considering for each pair only those batches where both clusters fired at least min_spikes_per_batch times
        ccg_frac_spikes (:, :) single 
        
        remove_duplicate_cluster_ids (:, 1) uint32 % list of cluster_ids to be removed
        duplicate_of_cluster_id (:, 1) uint32 % because they were determined to be duplicates of these clusters
        duplicate_spike_frac (:, 1) single % the corresponding entry in ccg_frac_spikes
    end
    
    methods
        function cda = ClusterDeduplicationAnalysis(ks_non_dedup, cluster_ids_retained)
            assert(~ks_non_dedup.is_deduplicated, 'ClusterDeduplicationAnalysis should be performed on non-deduplicated KilosortDataset');
            assert(ks_non_dedup.isLoadedBatchwise, 'KilosortDataset must have batchwise info loaded');
            cda.ks = ks_non_dedup;
            cda.cluster_ids = cluster_ids_retained;
        end
       
        function detect_duplicates(cda, varargin)
            p = inputParser();
            p.addParameter('min_spikes_per_batch', 1, @isscalar);
            p.addParameter('thresh_frac_spikes', 0.5, @isscalar);
            p.parse(varargin{:});
            
            ks = cda.ks; %#ok<*PROPLC>
            min_spikes_per_batch = p.Results.min_spikes_per_batch;
            thresh_frac_spikes = p.Results.thresh_frac_spikes;
            
            batch_inds = ks.compute_which_batch(ks.spike_times);
            [tf, cluster_inds] = ismember(ks.spike_clusters, cda.cluster_ids); % Critical to use our cluster ids since only a subset of ks.cluster_ids will be retained.
            assert(all(tf), 'Some spikes were retained from dropped cluster_ids');
            
            [ccg_counts_batch, cluster_counts_batch] = cda.compute_ccg_batchwise_overlapping(ks.spike_times, cluster_inds, batch_inds);
            
            nClusters = numel(cda.cluster_ids);
            
            cluster_total_spikes = sum(cluster_counts_batch, 2);
            
            % cluster_counts_batch is nClusters x nBatches. Want to create a cluster-pair by batch occupancy grid nBatch x nClusters x nClusters
            % so we dot multiply nBatch x nClusters by nBatch x 1 x nClusters
            suff_spikes = (cluster_counts_batch > min_spikes_per_batch)'; % nBatches x nClusters
            suff_spikes_pairwise = suff_spikes & permute(suff_spikes, [1 3 2]); % nBatches x nClusters x nClusters
            
            % ccg_counts_batch(iBin, iBatch, ci, cj)
            % 1. dot multiply counts with suff_spikes_pairwise
            % 2. sum across batches
            % 3. max along bins axis
            % 4. drop leading 2 singleton dims --> nClusters x nClusters
            ccg_counts_total = shiftdim(max(sum(ccg_counts_batch .* shiftdim(suff_spikes_pairwise, -1), 2), [], 1), 2);; % nBatches x nClusters x nClusters --> nClusters x nClusters
            
            % then construct a tensor of number of c1 spikes per batch, zeroing entries where (c1,c2) did not have suff_spikes_pairwise
            % 1. replicate c1 counts and dot multiply with suff_spikes_pairwise
            % 2. sum across batches
            % 3. drop leading singleton dim --> nClusters x nClusters
            c1_counts_pairwise = shiftdim(sum(repmat(cluster_counts_batch', [1 1 nClusters]) .* suff_spikes_pairwise, 1), 1); % nBatches x nClusters x nClusters --> nClusters x nClusters
            
            % ccg_frac_spikes(c1, c2) now indicates the fraction of c1 spikes co-occurring with a reliably timed c2 spike
            % considering for each pair only those batches where both clusters fired at least min_spikes_per_batch times
            % we threshold this to find potential duplicates, and then pick which cluster to drop. If both ccg_frac_spikes(c1, c2)
            % and ccg_frac_spikes(c2, c1) are above threshold, we drop the cluster with fewer spikes overall
            ccg_frac_spikes = ccg_counts_total ./ c1_counts_pairwise;
            cda.ccg_frac_spikes = ccg_frac_spikes;
            
            % threshold above thresh_frac_spikes to find duplicate pairs
            ccg_frac_spikes_above_thresh = ccg_frac_spikes > thresh_frac_spikes;
            [ic1v, ic2v] = find(ccg_frac_spikes_above_thresh);
            [~, sortIdx] = sort(ccg_frac_spikes(ccg_frac_spikes_above_thresh), 'descend');
            ic1v = ic1v(sortIdx);
            ic2v = ic2v(sortIdx);
            
            mask_cluster_as_dup = false(nClusters, 1);
            as_duplicate_of = zeros(nClusters, 1, 'uint32');
            dup_frac_spikes = nan(nClusters, 1, 'single');
            for iPair = 1:numel(ic1v)
                ic1 = ic1v(iPair);
                ic2 = ic2v(iPair);
                if mask_cluster_as_dup(ic1) || mask_cluster_as_dup(ic2)
                    % one of these is already a duplicate, we can ignore this pair
                    continue;
                end
                
                if ~ccg_frac_spikes_above_thresh(ic2, ic1)
                    % only c1's spikes are reliably accompanied by c2 spikes. so we remove c1 since c2 has additional 
                    % spikes that are not duplicated by c1
                    mask_cluster_as_dup(ic1) = true;
                    as_duplicate_of(ic1) = cda.cluster_ids(ic2);
                    dup_frac_spikes(ic1) = ccg_frac_spikes(ic1, ic2);
                else
                    % consider both directions, since both spikes are reliably time-locked with each other
                    % to decide, remove the cluster with fewer spikes overall
                    if cluster_total_spikes(ic1) > cluster_total_spikes(ic2)
                        mask_cluster_as_dup(ic2) = true;
                        as_duplicate_of(ic2) = cda.cluster_ids(ic1);
                        dup_frac_spikes(ic2) = ccg_frac_spikes(ic2, ic1);
                    else
                        mask_cluster_as_dup(ic1) = true;
                        as_duplicate_of(ic1) = cda.cluster_ids(ic2);
                        dup_frac_spikes(ic1) = ccg_frac_spikes(ic1, ic2);
                    end
                end
            end
            
            cda.remove_duplicate_cluster_ids = cda.cluster_ids(mask_cluster_as_dup);
            cda.duplicate_of_cluster_id = as_duplicate_of(mask_cluster_as_dup);
            cda.duplicate_spike_frac = dup_frac_spikes(mask_cluster_as_dup);
            
            debug('Removing %d clusters as potential duplicates with min frac %.2f\n', nnz(mask_cluster_as_dup), min(cda.duplicate_spike_frac));
        end
        
        function [ccg_counts_batch, cluster_counts_batch] = compute_ccg_batchwise_overlapping(cda, times, cluster_inds, batch_inds, varargin)
            % ccg_counts_batch(iBin, iBatch, ci, cj) = number of spikes fired by cj that occur within binWindow(iBin) of spikes from ci in batch iBatch
            % cluster_counts_batch(ci, iBatch) = total number of spikes fired by ci in batch iBatch
            %
            % in batch ibatch that co-occur just after a spike from cluster cj (<= within_samples after)
            %
            % here we set up a situation where cl 2 is always detected after cl 1, with a lag of 0 or 1
            % time | 0 | 1 | 2 | 3 | 4 | 5
            % cl 1 | x |   | x |   | x | . 
            % cl 2 |   | x |   | x | x |  
            %
            % direct ccg:
            %   time | -2| -1| 0 | 1 | 2 |
            % c 1->2 | 0 | 2 | 1 | 2 | 1 |  / 3 total
            % c 2->1 | 1 | 2 | 1 | 2 | 0 |  / 3 total
            %
            % overlapping ccg (where each bin is +/-1 sample, 
            % but only incremented by 1 if any c2 spike occurred in that offset range)
            % min    | -3| -2| -1| 0 | 1 |
            % max    | -1|  0| 1 | 2 | 3 |
            % c 1->2 | 2 | 2 | 3 | 3 | 2 |  / 3 total 
            % c 2->1 | 2 | 3 | 3 | 2 | 2 |  / 3 total

            p = inputParser();
            p.addParameter('nClusters', numel(cda.cluster_ids), @isscalar);
            p.addParameter('within_samples', 10, @isscalar);
            p.addParameter('bin_overlap', [-1 1], @isvector); % slop in the bins
            %p.addParameter('batch_min_spikes', 10, @isscalar);
            p.parse(varargin{:});
            
            nClusters = p.Results.nClusters;
            nBatches = max(batch_inds);
            within_samples = p.Results.within_samples;
            bin_overlap = p.Results.bin_overlap;
            bin_overlaps = bin_overlap(1):bin_overlap(2);
            bins = -within_samples:within_samples;
            nBins = numel(bins);
            ind_bin_zero = find(bins == 0, 1);
            
            batchwise_times = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(times, 1, nBatches, batch_inds);
            batchwise_cluster_inds = Neuropixel.Utils.TensorUtils.splitAlongDimensionBySubscripts(cluster_inds, 1, nBatches, batch_inds);

            cluster_counts_batch = zeros(nClusters, nBatches, 'single');
            ccg_counts_batch = zeros(nBins, nBatches, nClusters, nClusters);
            
            prog = ProgressBar(nBatches, 'Processing batchwise ccg');
            prog.enableParallel();
            parfor iBatch = 1:nBatches
                times_this = batchwise_times{iBatch};
                clusters_this = batchwise_cluster_inds{iBatch};
                cluster_counts_batch(:, iBatch) = accumarray(clusters_this, 1, [nClusters 1]);
                
                % use range search for all times against all other times
                % i1, i2 are indices into times_this, forming a pair
                idx = rangesearch(single(times_this), single(times_this), single(within_samples));
                if isempty(idx)
                    % no nearby spikes this batch
                    continue;
                end
                [i2, i1] = TensorUtils.catWhich(2, idx{:});
                i2 = i2';
                
                ibin = uint64(int64(ind_bin_zero) + int64(times_this(i2)) - int64(times_this(i1)));
                c1 = clusters_this(i1);
                c2 = clusters_this(i2);
                % range search occasionally returns times with offsets == dt_samples which would land in the immediately
                % flanking bins, but anything further out than that is problematic
%                 if any(ibin < -3 | ibin > nBins+4)
%                     error('Issue with bins on batch %d', iBatch);
%                 end
                mask_keep = i1 ~= i2 & ibin >= 1 & ibin <= nBins;
               
                if ~any(mask_keep)
                    continue;
                end
                % now we construct a tbl that will later be accumulated into the batch ccg
                % now we add duplicate entries for each offset allowed by bin_overlap 
                % and then uniquify by the c1 spike, c2 identity, and ccg bin assignment so that each 
                % source spike can contribute at most 1 pair to each bin even if multiple c2 spikes are eligible
                tbl = cell(numel(bin_overlap), 1);
                for ibin_overlap = 1:numel(bin_overlaps) 
                    offset = bin_overlaps(ibin_overlap);
                    tbl{ibin_overlap} = [i1(mask_keep), c1(mask_keep), c2(mask_keep), ibin(mask_keep) + offset];
                end

                tbl = unique(cat(1, tbl{:}), 'rows', 'stable');
                tbl = tbl(tbl(:, 4) >= 1 & tbl(:, 4) <= nBins, :);
                
                % i1, c1, c2, ibin --> ibin, c1, c2
                subs = tbl(:, [4 2 3]);
                ccg_counts_batch(:, iBatch, :, :) = reshape(accumarray(subs, 1, [nBins, nClusters, nClusters]), [nBins, 1 nClusters, nClusters]);
                prog.update(iBatch); %#ok<PFBNS>
            end
            prog.finish();
        end
            
            
    end
end