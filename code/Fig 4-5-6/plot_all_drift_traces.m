addpath(genpath(pwd))



dataPath = '/data/';
derivedDataPath = '/data/drift/';
figPath = '/figs/';


%% Plot drift traces for all files

matFiles = dir(fullfile(derivedDataPath, '*.mat'));

for ii = 1:length(matFiles)
    
    fullPath = fullfile(derivedDataPath, matFiles(ii).name);
    [~,session_name,~] = fileparts(matFiles(ii).name);
    
    rez = load(fullPath);

    tvec = rez.secPerBatch * (1:rez.nBatch);
    
    figh = figure(ii); clf;

    % remove top and bottom of probe - since drift estimates there are less
    % accurate
    dshift = rez.dshift(:,1:end);

    % dshift = mean(dshift,2);
    plot_drift_traces3(dshift, tvec, 'Estimator','')
    
    fname = fullfile(figPath, ['drift_map_' session_name]);
    print(figh, fname,'-dpdf','-painters','-bestfit')

end