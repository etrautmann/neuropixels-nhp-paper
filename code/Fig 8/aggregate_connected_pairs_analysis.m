%% EMT 2023-01-14
% script to aggregate the output from walthru_connected_pairs_emt applied
% to individual datasets. 
%
% Todo: somewhat large variance in % connected pairs across datasets, would
% be interesting to look into


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

%% list of all datasets

[datasets] = NeuropixelsNhpAnalysisDatasetList()

%% ========================================================================
% Aggregated analyses across all datasets
% ========================================================================= 


%% Load datasets


datesList = {'2021-03-18','2021-03-24','2021-03-26','2021-03-31','2021-05-07', '2021-05-10', '2021-05-18','2021-05-20','2021-05-27','2021-08-09','2021-08-11','2021-08-14','2021-08-15'}
% datesList = {'2021-05-25'}


inds = ismember(datasets.date, datesList);
datasetsSelect = datasets(inds,:)

cpas = {};
kss = {};

for iSet = 1:size(datasetsSelect,1)
    ks_path = datasetsSelect.ksOutputPath{iSet}
    [cpa, ks, ~] = run_connected_pairs_analysis_v1(ks_path, mapfile,'save_results', true, 'use_cached', true);
    
    cpas{end+1} = cpa;
    kss{end+1} = ks;

end


%% Calculate % of connected pairs for each dataset

nConnected = [];
nGoodClusters = [];
nTotals = [];
connectedPairProbs = [];

for iSet = 1:size(datasetsSelect,1)
    
    cpa = cpas{iSet};
    ks = kss{iSet};

    nGoodCluster = numel(ks.clusters_good);
    nPossiblePairs = (nGoodCluster * (nGoodCluster-1))/2;
    connectedPairProb = cpa.nPairsConnected/nPossiblePairs;

    % pack outputs
    nConnected(iSet) = cpa.nPairsConnected;
    nGoodClusters(iSet) = nGoodCluster;
    nTotals(iSet) = nPossiblePairs;
    connectedPairProbs(iSet) = connectedPairProb;

end


nGoodClusters
nTotals
connectedPairProbs

disp(sprintf('mean connected pair probability: %2.2f percent',100*mean(connectedPairProbs)))
disp(sprintf('std connected pair probability: %2.2f percent',100*std(connectedPairProbs)))
disp(sprintf('# datasets: %d',height(datasetsSelect)))

%%

figure(1); clf;
histogram(connectedPairProbs,25)

figure(2); clf;

plot(nGoodClusters, nConnected,'o')


%%