function [ units ] = kilosort_read_spikes_for_marius( data_folder, NumChannels )
% [ spike_time, unit_id ] = load_spike_kilosort( data_folder )
% load spike time from output of Kilosort2 and PhyGUI.

if nargin<1
    data_folder = [];
end

%% files to read
spike_time_file = fullfile(data_folder, 'spike_times.npy');

if ~exist(spike_time_file, 'file')
    units = [];
    return;
end

spike_id_file = fullfile(data_folder, 'spike_clusters.npy'); % updated by Phy after manual refine;
cluster_group_file = fullfile(data_folder, 'cluster_group.tsv');
params_file = fullfile(data_folder, 'params.py');
template_file = fullfile(data_folder, 'templates.npy');
spike_template_file = fullfile(data_folder, 'spike_templates.npy');
channel_map_file = fullfile(data_folder, 'channel_map.npy');

fileName = fullfile(data_folder,'*.bin');
filenamestruct = dir(fileName);
if ~isempty(filenamestruct) && ~isnan(NumChannels)
    fileName = fullfile(data_folder,filenamestruct(1).name);
    filenamestruct = dir(fileName);
    dataTypeNBytes = numel(typecast(cast(0, 'int16'), 'uint8')); % determine number of bytes per sample
    nSamp = filenamestruct.bytes/(NumChannels*dataTypeNBytes);  % Number of samples per channel
    mmf = memmapfile(fileName, 'Format', {'int16', [NumChannels nSamp], 'x'});
    use_templates = false;
else
    use_templates = true;
end
%% load data
sampling_rate = kilosort_read_param(params_file, 'sample_rate', '%f');
spike_time_ind = readNPY(spike_time_file);
spike_id = readNPY(spike_id_file);

spike_time = double(spike_time_ind)/sampling_rate; % sec

%% load cluster group
if exist(cluster_group_file, 'file')
    [unit_id, unit_quality, is_noise] = kilosort_read_cluster_group(cluster_group_file);
    has_spikes =  ismember(unit_id,unique(spike_id));
    if ~all(has_spikes)
        fprintf('Warning: Some units do not have spikes somehow and are excluded.')
    end
    has_spikes =  ismember(unit_id,unique(spike_id));
    b_has_1000_spikes = false(size(unit_id));
    for unit_index = 1:length(unit_id)
        b_has_1000_spikes(unit_index) = sum(spike_id==unit_id(unit_index))>1000;
    end
    
    unit_id = unit_id(b_has_1000_spikes & has_spikes & ~is_noise);
    unit_quality = unit_quality(b_has_1000_spikes & has_spikes & ~is_noise);
else
    error('Cannot find cluster group file');
end
n_unit = length(unit_id);

%% find the maximum channel of each units
template = readNPY(template_file); %[nTemplates, nTime, nChannels]
spike_template = readNPY(spike_template_file);
channel_map = readNPY(channel_map_file);

% get max ch of templates
amplitude = sum(template.^2, 2);
[~, ind] = max(amplitude, [], 3);
template_channel = channel_map(ind);

% find template of first spike of each unit and thus the max ch
unit_ch = zeros(n_unit,1);
unit_template = zeros(n_unit,size(template,2)); %templates
for i=1:n_unit
    sti = find(spike_id == unit_id(i), 1, 'first' ); % spike time index
    ti = spike_template(sti) + 1; % template index, zero based
    unit_ch(i) = template_channel(ti) + 1;
    unit_template(i,:) = template(ti,:,ind(ti));
end

%% construct data structure of unit info

units = struct('unit_id', [], 'spike_time', [], 'quality', [], 'channel', [], 'waveforms', []);

for i = 1:n_unit
    units(i).unit_id = unit_id(i);
    units(i).spike_time = spike_time( spike_id == unit_id(i));
    units(i).quality = unit_quality{i};
    units(i).channel = unit_ch(i);

end


end