%% for RigJ recording machine


clc
clear all
saveFigs = 1;

% machine specific paths
setenv('NHP_PIXEL_COMPUTED_ROOT', 'H:\data\pacman-task\');
setenv('NHP_PIXEL_DATA_ROOT', 'H:\data\pacman-task\')
setenv('FIG_ROOT','C:\Users\User\Dropbox\rigJ_figures')
setenv('NEUROPIXEL_MAP_FILE','neuropixNHP_kilosortChanMap_v1.mat')
coderoot = 'C:\Users\User\code\';

% add libraries
addpath(genpath(fullfile(coderoot,'analysis','neuropixelNHP','matlab','connected_pairs')))
addpath(genpath(fullfile(coderoot,'matlab-utils')))
addpath(genpath(fullfile(coderoot,'trial-data')))
addpath(genpath(fullfile(coderoot,'neuropixel-utils')))
addpath(genpath(fullfile(coderoot,'spikes')))

% suppress warnings
warning('off', 'MATLAB:MKDIR:DirectoryExists');


%% Create the kilosort dataset

% ks_path = 'H:\data\pacman-task\cousteau\2021-03-18\neuropixels\pacman-task_c_210318_neu_g0\pacman-task_c_210318_neu_g0_imec0'
% ks_path = 'H:\data\pacman-task\cousteau\2021-08-09\neuropixels\pacman-task_c_210809_neu_g0\pacman-task_c_210809_neu_g0_imec0'
ks_path = 'H:\data\pacman-task\cousteau\2021-05-25\neuropixels\pacman-task_c_210525_neu_g0\pacman-task_c_210525_neu_g0_imec0'

% ks_path = 'H:\data\pacman-task\cousteau\2021-03-24\neuropixels\pacman-task_c_210324_neu_g0\pacman-task_c_210324_neu_g0_imec0'
% ks_path = 'H:\data\pacman-task\cousteau\2021-03-26\neuropixels\pacman-task_c_210326_neu_g0\pacman-task_c_210326_neu_g0_imec0';
% ks_path = 'H:\data\pacman-task\cousteau\2021-08-11\neuropixels\pacman-task_c_210811_neu_g0\pacman-task_c_210811_neu_g0_imec0';
% ks_path = 'H:\data\pacman-task\cousteau\2021-08-14\neuropixels\pacman-task_c_210814_neu_g0\pacman-task_c_210814_neu_g0_imec0';
% ks_path = '';

% having an imec dataset is optional
% imec = Neuropixel.ImecDataset(paths.npixLfpPath, channelMap=channelMap); % or npixApPath
% ks_args.imecDataset = imec;

mapfile = fullfile(coderoot,'analysis','runKilosort','neuropixelMapFiles','neuropixNHP_kilosortChanMap_v1.mat');
% channelMap = Neuropixel.ChannelMap(mapfile);

%% Run analysis (or load results if already run)
clc

[cpa, ks, resultsfile] = run_connected_pairs_analysis_v1(ks_path, mapfile,'save_results', true, 'use_cached', true);


%%


cpi = cpa.connected_pair_ids


%% Find some local patch of densely connected neurons mannually (could automate this)


metrics = ks.computeMetrics();
mask = ismember(cpa.cluster_ids, ks.clusters_good)

connected = cpa.connected(mask,mask)
depth = metrics.cluster_centroid(mask,2)
[b,ix] = sort(depth,'ascend')

figh = figure(4); clf;

imagesc(~connected(ix,ix)- eye(size(connected))*.5)
axis square
axis on
colormap gray



%% find cluster inds from looking for an interesting dense patch on the connected pairs matrix

cluster_inds = 130:132;
cg = ks.clusters_good;
pair_ids = cpa.connected_pair_ids;


figure(3); clf;
clc

nClust = length(cluster_inds)
ind = 1
for ii = 1:nClust
    for kk = 1:nClust
        subplot(nClust, nClust,ind)
        ind = ind+1;
        
        c1 = cg(cluster_inds(ii));
        c2 = cg(cluster_inds(kk));
        
        cpa.plotCCG(c1,c2)
        
    end        
end






