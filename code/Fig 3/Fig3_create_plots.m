datapath = '/Users/erictrautmann/Dropbox/Columbia/manuscripts - Columbia/paper - Neuropixels NHP/code and data/data/Fig 3/pacman-task_c_210318_taskdata_cond_5.mat'

load(datapath)


%% plot force


figure(1); clf;

plot(data.tvec, data.targetForce)


%% plot PSTHs

figure(2); clf

unitList = [43, 231, 232, 112, 125, 146, 197, 288];

plot(data.tvec, data.psths(:, unitList))


%% 