session_name = '210803_133520_Alfie';

imec = 0; %Probe index
% sorted_units_folder = ['Z:\Data' filesep session_name '_g0' filesep session_name '_g0_imec' num2str(imec) filesep 'RAW']; %Automatically define based on session_name.
sorted_units_folder = ['.' filesep 'data' filesep session_name '_g0' filesep session_name '_g0_imec' num2str(imec) filesep]; %Automatically define based on session_name.

addpath(genpath(pwd))
paradigm_data = load(['.' filesep 'data' filesep 'fob_paradigm_data_' session_name '_imec' num2str(imec) '.mat']);
peri_stimulus_range = -200:400; %Range to extract for raster

[ units ] = kilosort_read_spikes_minimal( sorted_units_folder, nan); %function to read out kilosort results
units = units(~isnan(cat(1,units.unit_id))); %Ignore any nan units
n_units = length(units);
n_times = ceil(max(cat(1,units.spike_time))*1000);

%Create boolean time course of entire experiment from spike times
timecourse = false(n_units,n_times);
for unit_index = 1:n_units
    timecourse(unit_index,round(units(unit_index).spike_time*1000)) = true;
end

% %Plot time course of firing rates for all units
% figure(1);imagesc(zscore(smoothdata(timecourse','movmean',10000),[],1)')

%Compute rasters from timecourse
fob_raster = computeRaster(round(paradigm_data.stimulus_on_times_in_kilosort_time*1e3),timecourse,-peri_stimulus_range(1),peri_stimulus_range(end));

time_window_response = [50, paradigm_data.on_duration + 100];
time_window_baseline = [0, 50];
start_index_response = find(peri_stimulus_range==time_window_response(1));
end_index_response = find(peri_stimulus_range==time_window_response(2));
start_index_baseline = find(peri_stimulus_range==time_window_baseline(1));
end_index_baseline = find(peri_stimulus_range==time_window_baseline(2));
responses_avg_across_time = squeeze(mean(fob_raster(:,start_index_response:end_index_response,:),2)); %n_trials x n_units
baselines_avg_across_time = squeeze(mean(fob_raster(:,start_index_baseline:end_index_baseline,:),2)); %n_trials x n_units
unit_baselines = mean(baselines_avg_across_time,1); %n_trials x n_units
nstimuli = max(paradigm_data.stimulus_index_valid);

%Compute orange matrices (average responses with rows as cells and stimuli as columns)
orange_matrix = nan(size(responses_avg_across_time,2),nstimuli);

for stimulus_index = 1:nstimuli
    stimulus_index_trials = find(paradigm_data.stimulus_index_valid == stimulus_index);
    orange_matrix(:,stimulus_index)= mean(responses_avg_across_time(stimulus_index_trials,:),1);
end


orange_matrix_minus_baseline_normalized = subtract_baselines_and_normalize(orange_matrix,unit_baselines');


%Plot results
figure(3);
orange_image = plot_orange_plot(orange_matrix_minus_baseline_normalized);
