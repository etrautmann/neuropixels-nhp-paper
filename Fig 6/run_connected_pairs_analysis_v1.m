function [cpa, ks, resultsfile] = run_connected_pairs_analysis_v1(ks_path, mapfile, varargin)
%   Wraps Dan O'Shea's NHPPixel.ConnectedPairsAnalysis
% 
%   inputs: 
%      ks_path:     path to kilosort output
%      save_path:   path for s
% 
% EMT 2022-10-21
% 2022-11-07: adding cluster de-duplication analysis preprocessor

p = inputParser();
p.addRequired('ks_path', @isstringlike);
p.addRequired('mapfile', @isstringlike);
p.addParameter('save_results', true, @islogical);
p.addParameter('use_cached', true, @islogical);     % if analysis has already run with saved results, load these instead of recomputing
p.addParameter('jitter_reps',5 , @isnumeric);
p.addParameter('plotConnectedPairsPanel',false,@islogical);
p.addParameter('plotLatencyHist',false,@islogical);
p.addParameter('plotConnectionMat',false,@islogical);
p.addParameter('plotConnectionDist',false,@islogical);



p.KeepUnmatched = false;
p.parse(ks_path, mapfile, varargin{:});


assert(isfolder(p.Results.ks_path))
assert(isfile(p.Results.mapfile))

% initialize KS object, set some default parameters here
% ks_args.deduplicate_spikes = true;
% ks_args.deduplicate_cutoff_spikes = true;

mapfile = p.Results.mapfile;
channelMap = Neuropixel.ChannelMap(mapfile);
ks = Neuropixel.KilosortDataset(ks_path, 'channelMap', channelMap, deduplicate_spikes=false, deduplicate_cutoff_spikes=false);
disp(sprintf('Running connected pairs analysis on dataset: %s\n',ks.pathLeaf))
resultsfile = fullfile(ks.path, sprintf('%s__connected_pairs.mat',ks.pathLeaf));



if isfile(resultsfile) && p.Results.use_cached
    disp(sprintf('CPA analysis already run, loading cached results from: \n%s\n',ks_path))
    load(resultsfile)
else
    % ks.deduplicate_within_samples = 5; % old, use ClusterDeduplicationAnalysis instead
    % ks.deduplicate_within_distance = 50; % old, use ClusterDeduplicationAnalysis instead

    ks.load(loadFeatures=false, loadBatchwise=true);
    ks.mask_clusters(ks.clusters_good);
    if numel(ks.clusters_good) < 50
        warning('only %d good clusters found in dataset: \n%s\nHas manual curation been performed?',numel(ks.clusters_good), ks.pathLeaf)
        return
    end

    % 1) cluster de-duplication preprocessor
    ntPerBatch = 65536;
    nBatches = ceil(double(max(ks.spike_times + uint64(100))) / ntPerBatch);
    ks.batch_sort_order = (1:nBatches)';
    ks.batch_starts = (uint64(1) : ntPerBatch : (ntPerBatch*uint64(ks.nBatches-1) + uint64(1)))';

    cda = ClusterDeduplicationAnalysis(ks, ks.clusters_good)

    delete(gcp('nocreate'))
    parpool(4); % or some reasonable number of cores, typically I don't use more than 8.
    cda.detect_duplicates

    % mask out duplicated clusters
    keep_clusters = ks.clusters_good(~ismember(ks.clusters_good, cda.remove_duplicate_cluster_ids));
    ks.mask_clusters(keep_clusters);

    % 2) run Connected pairs analysis
    cpa = NHPPixel.ConnectedPairsAnalysis(ks, [],'jitter_reps', p.Results.jitter_reps);

    % Compute metrics
    cpa.computeSmoothedCCGs();
    cpa.findConnectedPairs;


    if p.Results.save_results
        try save(resultsfile, 'cpa', 'ks', '-v7.3')
            disp(sprintf('saving output to: \n%s\n',resultsfile))
        catch
            disp(sprintf('cannot save output to: \n%s\n',resultsfile))
        end
    end
end


%% plot results

if p.Results.plotConnectedPairsPanel


end

return


