% Script to run connected pairs functional connectivity analysis -
% originally follwing code from Dan O'Shea
%
% EMT 2022-11-07
% EMT 2023-01-14

% todo: refactor to push plotting functions into
% run_connected_pairs_analysis_v1.m 

clc
clear all
saveFigs = 1;

% required libraries:
%   - djoshea/trial-data
%   - djoshea/neuropixel-utils
%   - djoshea/matlab-utils (maybe?)
%   - cortex-lab/spikes
%   - nhp-pixel (zipped)

 %% machine specific paths
[~,name]= unix('hostname');
name = strtrim(name);

if strcmp(name, 'Erics-MBP-2')  % for local MBP
    coderoot = '/Users/erictrautmann/Dropbox/shenoy_lab/code';
    setenv('NHP_PIXEL_COMPUTED_ROOT', '/Volumes/emt_ssd_6/data/pacman-task');
    setenv('DATA_ROOT', '/Volumes/emt_ssd_6/data/pacman-task')
    setenv('FIG_ROOT','/Users/erictrautmann/Dropbox/columbia/manuscripts - columbia/Neuropixels NHP paper/figures/Fig. 7 - Recording quality metrics/assets/')
    setenv('NEUROPIXEL_MAP_FILE','neuropixNHP_kilosortChanMap_v1.mat')

elseif strcmp(name, 'User-PC') % for RigJ recording machine
    coderoot = 'C:\Users\User\code\';
    setenv('NHP_PIXEL_COMPUTED_ROOT', 'H:\data\pacman-task\');
    setenv('NHP_PIXEL_DATA_ROOT', 'H:\data\pacman-task\')
    setenv('FIG_ROOT','C:\Users\User\Dropbox\rigJ_figures')
    setenv('NEUROPIXEL_MAP_FILE','neuropixNHP_kilosortChanMap_v1.mat')
  
end

% add libraries
addpath(genpath(fullfile(coderoot,'analysis','neuropixelNHP','matlab','connected_pairs')))
addpath(genpath(fullfile(coderoot,'matlab-utils')))
addpath(genpath(fullfile(coderoot,'trial-data')))
addpath(genpath(fullfile(coderoot,'neuropixel-utils')))
addpath(genpath(fullfile(coderoot,'spikes')))

% suppress warnings
warning('off', 'MATLAB:MKDIR:DirectoryExists');

figPath = fullfile(getenv('FIG_ROOT'),'connected pairs analysis');
mkdir(figPath);
disp(sprintf('Saving figures to: \n%s\n',figPath))


%% Create the kilosort dataset

mapfile = fullfile(coderoot,'analysis','runKilosort','neuropixelMapFiles','neuropixNHP_kilosortChanMap_v1.mat');
assert(isfile(mapfile))


%% get list of all datasets

[datasets] = NeuropixelsNhpAnalysisDatasetList()

%% Iterate over selected datasets

datesList = {'2021-03-18','2021-03-24','2021-03-26','2021-03-31','2021-05-07', '2021-05-10', '2021-05-18','2021-05-20','2021-05-27','2021-08-09','2021-08-11','2021-08-14','2021-08-15'}
% datesList = {'2021-03-18','2021-05-18','2021-08-09'}
% datesList = {'2021-03-24'}
% datesList = {'2021-03-26','2021-03-31'}
% datesList = {'2021-03-24','2021-05-07','2021-05-20','2021-05-27', '2021-08-14','2021-08-15'}
% datesList = {'2021-05-10','2021-08-11'}

inds = ismember(datasets.date, datesList)
datasetsSelect = datasets(inds,:)

for iSet = 1:size(datasetsSelect,1)
    ks_path = datasetsSelect.ksOutputPath{iSet}
    [cpa, ks, resultsfile] = run_connected_pairs_analysis_v1(ks_path, mapfile,'save_results', true, 'use_cached', false);
end


%% 

return





%% ========================================================================
% Analysis for single dataset 
% ========================================================================= 


figPath = fullfile(getenv('FIG_ROOT'),datasets.subject{ii}, datasets.date{ii},'/');
    mkdir(figPath)


% figure(1); clf;
% 
% pairNum = 125;
% cpa.plotConnectedPairCCG(pairNum,'style','stem');
% 
% pair = cpa.connected_pair_ids(pairNum,:);
% % 
% % figure(2); clf;
% % cpa.plotConnectedPairCCG_rawVsJittered(pairNum)
% 
% %%
% 
% xlabel('latency')
% title(sprintf('clusters %d-%d.pdf', pair(1), pair(2)))
% 
% 
% figNameFull = fullfile(figPath,sprintf('%s - Connected pairs panel - clusters %d-%d',ks.pathLeaf, pair(1), pair(2)))
% print(figNameFull,'-dpng')





%% ========================================================================
% Plot connected pairs panel 
% ========================================================================= 


clc
nX = 4;
nY = 6;

nPlots = ceil(cpa.nPairsConnected/(nX*nY));
% nPlots = 1
for kk = 1:nPlots
    ind = 1;
    
    figh = figure(1); clf;
    set(figh,'color','w');
    set(figh,'WindowStyle','normal')
    set(figh,'position',[0 0 1200 1800])
    
    for ii = 1:nY
        for jj = 1:nX
            pairInd = (kk-1)*nX*nY + ind
            
            if pairInd <= cpa.nPairsConnected
                subplot(nY,nX,ind)
                cpa.plotConnectedPairCCG(pairInd,'style','area');
                ind = ind+1;
                xlim([-25 25])
                title(sprintf('%s',num2str(cpa.connected_pair_ids(pairInd,:))))
                box off
%                 makepretty
                
            end
        end
    end
    
    sprintf('New panel\n\n')
 
    if saveFigs
        
        figNameFull = fullfile(figPath,sprintf('%s - Connected pairs panel - %d.pdf',ks.pathLeaf,kk))
        print(figNameFull,'-dpdf','-fillpage','-painters')
    end
end



%% Plot latency histogram

figh = figure(3); clf;

hi = histogram(cpa.connected_pair_latency,0:.5:10)
set(hi,'edgecolor','none')
box off

title('Connected pair spiking timing latency')
xlabel('Time (ms)')
ylabel('count')
makePrettyAxis

figNameFull = fullfile(figPath,sprintf('%s - Connected pair spike timing latency.pdf',ks.pathLeaf))
print(figNameFull,'-dpdf','-bestfit','-painters')

%% Plot connection matrix

metrics = ks.computeMetrics();
mask = ismember(cpa.cluster_ids, ks.clusters_good)

connected = cpa.connected(mask,mask)
depth = metrics.cluster_centroid(mask,2)
[b,ix] = sort(depth,'ascend')

figh = figure(4); clf;

imagesc(~connected(ix,ix))
axis square
axis on
box off
colormap gray

set(figh,'windowstyle','normal')
set(figh,'position',[0 0 500 500])

figNameFull = fullfile(figPath,sprintf('%s - Connection matrix.pdf',ks.pathLeaf))
print(figNameFull,'-dpdf','-bestfit','-painters')


%% Plot connection distance (version 1 - too verbose)

% cluster_list = ks.clusters_good;
% pair_ids = cpa.connected_pair_ids;
% centroids = metrics.cluster_centroid(mask,:);
% 
% distances = [];
% 
% for iPair = 1:size(pair_ids,1)
%     c1 = pair_ids(iPair,1);
%     c2 = pair_ids(iPair,2);
%     
%     centroid1 = centroids(find(cluster_list==c1),:);
%     centroid2 = centroids(find(cluster_list==c2),:);
%     
%     distances(iPair) = sqrt((centroid1(1) - centroid2(1))^2 + (centroid1(2) - centroid2(2))^2);
% end
% 
% figure(4); clf;
% hi = histogram(distances,50);
% set(hi,'edgecolor','none')
% xlabel('Distance (um)')
% ylabel('# connected pairs')
% makepretty


%% Plot connection distance, version 2

cc_dist = metrics.computeClusterDistanceMatrix();
pair_ids = cpa.connected_pair_ids;

% iterate over pairs and pull out connection distance, 
distances = [];
for iPair = 1:size(pair_ids,1)
    inds = metrics.lookup_clusterIds(pair_ids(iPair,:));
    distances(iPair) = cc_dist(inds(1), inds(2));
end

% plot results
figh = figure(5); clf;
hi = histogram(distances,50);
set(hi,'edgecolor','none')
xlabel('Distance (um)')
ylabel('# connected pairs')
makepretty

set(figh,'windowstyle','normal')
set(figh,'position',[0 100 500 350])


figNameFull = fullfile(figPath,sprintf('%s - hist connections vs. distance.pdf',ks.pathLeaf));
print(figNameFull,'-dpdf','-bestfit','-painters')


%% Plot positions and connection of connected pairs.

cluster_list = ks.clusters_good;
pair_ids = cpa.connected_pair_ids;
centroids = metrics.cluster_centroid(mask,:);

figh = figure(6); clf;
for iPair = 1:size(pair_ids,1)
    c1 = pair_ids(iPair,1);
    c2 = pair_ids(iPair,2);
    
    centroid1 = centroids(find(cluster_list==c1),:);
    centroid2 = centroids(find(cluster_list==c2),:);
    
    plot([centroid1(1) centroid2(1)], [centroid1(2) centroid2(2)],'.-','color',[.3 .3 .3 .3],'markersize',25)
    hold on
end

xlabel('X (um)')
ylabel('Y (um)')
title(sprintf('connected pairs - %s',ks.pathLeaf),'interpreter','none')
makePrettyAxis

set(figh,'windowstyle','normal')
set(figh,'position',[0 0 400 1000])

figNameFull = fullfile(figPath,sprintf('%s - connected pair positions.pdf',ks.pathLeaf))
print(figNameFull,'-dpdf','-bestfit','-painters')




%%




%%
return

%% ========================================================================
%   
% =========================================================================



%% ========================================================================
% Sandbox
% % =========================================================================

%% useful functions

% lookup index from cluster (note that this is index within all clusters,
% not only good units

metrics.lookup_clusterIds(ks.clusters_good);


% 
% 
% ks_args.deduplicate_spikes = true;
% ks_args.deduplicate_cutoff_spikes = true;
% ks.deduplicate_within_samples = 5;
% ks.deduplicate_within_distance = 50;
% ks = Neuropixel.KilosortDataset(ks_path, 'channelMap', channelMap, ks_args);
% 
% ks.load(loadFeatures=false);
% 
% %% Test connected pairs analysis on small subset of data
% cpa2 = NHPPixel.ConnectedPairsAnalysis(ks, [],'jitter_reps', 1);
% 
% 
% %% method 1 - specify clusters (currently throws errors)
% 
% 
% cpa2 = NHPPixel.ConnectedPairsAnalysis(ks, [],'jitter_reps', 1, 'cluster_ids',ks.clusters_good(1:10));
% 
% 
% %% method 2 - mask neurons in ks object
% 
% clusterList = ks.clusters_good
% ks.mask_clusters(clusterList(1:10))
% 
% cpa2 = NHPPixel.ConnectedPairsAnalysis(ks, [],'jitter_reps', 1);
% 
% 
% %%
% cpa2 = NHPPixel.ConnectedPairsAnalysis(ks, [],'jitter_reps', 1);
% 
% 
% 
% 
% 
% 

%%  Old analysis below
% 
% % return



%% Run analysis on single dataset (or load results if already run)
% clc
% 
% [cpa, ks, resultsfile] = run_connected_pairs_analysis_v1(ks_path, mapfile,'save_results', true, 'use_cached', false);
% 



%% 
% ks_args.deduplicate_spikes = true;
% ks_args.deduplicate_cutoff_spikes = true;
% ks = Neuropixel.KilosortDataset(ks_path, 'channelMap', channelMap, ks_args);
% ks.deduplicate_within_samples = 5;
% ks.deduplicate_within_distance = 50;
% ks.load(loadFeatures=false);
% ks.mask_clusters(ks.clusters_good);
% 
% %% Run the CCG tool (this takes a long time, like 4 hours or so).
% 
% tic
% cpa = NHPPixel.ConnectedPairsAnalysis(ks, [],'jitter_reps', 5);
% toc
% 
% % You can save cpa to disk and it should have everything you need in its properties
% 
% filename = fullfile(ks.path, sprintf('%s__connected_pairs.mat',ks.pathLeaf))
% 
% %%
% save(filename, 'cpa','-v7.3')
% 
% %%
% % cpa = load(filename)
% 
% % load('E:\emt_working\sandbox\pacman-task_c_220703_neu_stitching_bank1_g0_imec0_cpa.mat')
% 
% load(filename)
% %%
% 
% cpa.computeSmoothedCCGs()
% cpa.findConnectedPairs
% 
% cpa.connected_pair_ids
% cpa.connected_pair_latency
% cpa.connected_pair_mag





%% ========================================================================
% dev - ClusterDeduplicationAnalysis 
% ========================================================================= 


% Dan's example code:
% close any parallel pool if it's open
delete(gcp('nocreate'))

channelMap = Neuropixel.ChannelMap(mapfile);

ks = Neuropixel.KilosortDataset(ks_path, channelMap=channelMap, deduplicate_spikes=false, deduplicate_cutoff_spikes=false);
ks.load(loadFeatures=false, loadBatchwise=true);
ks.mask_clusters(ks.clusters_good);

%ntPerBatch = uint64(ops.NT - ops.ntbuff)
ntPerBatch = 65536;
nBatches = ceil(double(max(ks.spike_times + uint64(100))) / ntPerBatch);
ks.batch_sort_order = (1:nBatches)';
ks.batch_starts = (uint64(1) : ntPerBatch : (ntPerBatch*uint64(ks.nBatches-1) + uint64(1)))';

cda = ClusterDeduplicationAnalysis(ks, ks.clusters_good)
parpool(2); % or some reasonable number of cores, typically I don't use more than 8.
cda.detect_duplicates


