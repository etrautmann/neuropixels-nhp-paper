addpath(genpath('/Users/erictrautmann/Dropbox/shenoy_lab/code/matlab-utils'))
addpath(genpath('/Users/erictrautmann/Dropbox/shenoy_lab/code/analysis/neuropixelNHP'))

% dataPath = '/Users/erictrautmann/data/pacman-task/cousteau_rez_mats/';

dataPath = '/Users/erictrautmann/data/Neuropixels_NHP_data/rez_mats';
derivedDataPath = '/Users/erictrautmann/data/pacman-task/cousteau_drift_data/';
figPath = '/Users/erictrautmann/Dropbox/Columbia/figures/neuropixelsNHP/drift_maps/';

%%
% Get a list of all .mat files in the directory
matFiles = dir(fullfile(dataPath, '*.mat'));

for ii = 1:length(matFiles)

    fullPath = fullfile(dataPath, matFiles(ii).name);
    load(fullPath);

    [~,session_name,~] = fileparts(matFiles(ii).name)

    dshift = rez.dshift;
    nBatch = rez.ops.Nbatch;
    secPerBatch = (rez.ops.sampsToRead / rez.ops.fs) / nBatch;

    fname = fullfile(derivedDataPath, ['drift_data_' session_name '.mat']);
    save(fname, 'dshift','nBatch','secPerBatch')

end


%%

matFiles = dir(fullfile(derivedDataPath, '*.mat'));

for ii = 1:length(matFiles)
    
    fullPath = fullfile(derivedDataPath, matFiles(ii).name);
    [~,session_name,~] = fileparts(matFiles(ii).name);
    
    rez = load(fullPath);

    tvec = rez.secPerBatch * (1:rez.nBatch);
    
    figh = figure(ii); clf;

    % remove top and bottom of probe - since drift estimates there are less
    % accurate
    dshift = rez.dshift(:,2:end-1);

    % dshift = mean(dshift,2);
    plot_drift_traces2(dshift, tvec)
    
    fname = fullfile(figPath, ['drift_map_' session_name]);
    print(figh, fname,'-dpdf','-painters','-bestfit')

end