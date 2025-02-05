%% Calculate CCG

%number of sessions
params.date = '20210928';

%ccg calculation variables
params.start_time = -0.6;
params.end_time = 0;
params.max_lag = 100;
params.min_lag = -100;
params.jit_window = 25;

%postprocess variables
params.max_pt_lag = 10;
params.min_pt_lag = -10;
params.split_layers = 1;
params.large_output = 1;

%sig criteria
params.sig_num_stds = 7; % must be 7 stds above noise mean
params.sig_max_lag = 10; % must have peak w/in 10 ms of 0
params.sig_min_std = 0; % must have a non-zero noise std

%cluster variables
params.rng_seed = 1;
params.tsne_dims = 3;
params.max_k = 10;
params.kmeans_rep = 50;
params.kmeans_maxiter = 100;
params.fix_k = true;
params.k = 4;

%where to save
params.filename = "";
params.output_dir = "output/" + params.date + "/";
params.data_dir = "data/";
params.ccg_output_file = "ccg_output";

params.process_cluster = ["", "mean_center", "zscore", "l2norm"];
params.process_params = 3;

params.large_output = 1;

group_idx = 1;
params.figure_dir  = "figures/";

if ~exist(params.output_dir, 'dir')
    mkdir(params.output_dir);
end


% CCG is computed using the package published by Trepka et al. eLife 2022 (https://github.com/et22/trepka_etal_elife_2022)
% ccg_output = analyze_ccgs(data, params);
% save(params.output_dir + params.ccg_output_file+".mat", 'ccg_output', '-v7.3');

%
load(params.output_dir + params.ccg_output_file+".mat");
ccgs = ccg_output.ccg_control(:,91:111);
pair_ids = ccg_output.neuron_id_pairs;

maxlag = params.max_pt_lag;

[peaks, peak_lag] = max(ccgs,[],2);
[troughs, trough_lag] = min(ccgs, [], 2);
[lpeaks, lpeak_lag] = max(ccg_output.ccg_control,[],2);
[ltroughs, ltrough_lag] = min(ccg_output.ccg_control,[],2);
ccg_output.config.cont_vars = ["peaks", "troughs", "peak_lag","trough_lag", "area", "peak_trough_diff", "peak_trough_lag_diff","pair_distance","peak_width", "trough_width"];
ccg_output.config.cat_vars = ["pre_id", "post_id", "pre_ct", "post_ct", "pre_cl", "post_cl", "pre_depth", "post_depth", "clust", "syn_peak", "syn_trough"];
ccg_output.config.params = params;    
ccg_output.peaks = peaks;
ccg_output.troughs = troughs;


ccg_output.peak_lag = 100-lpeak_lag+1;
ccg_output.trough_lag = 100-ltrough_lag+1;
    
ccg_output.peak_trough_diff =  ccg_output.peaks- ccg_output.troughs;
ccg_output.peak_trough_lag_diff = ccg_output.peak_lag-ccg_output.trough_lag;   
ccg_output.syn_peak = ccg_output.peak_lag == 0; 
ccg_output.syn_trough = ccg_output.trough_lag == 0; 

% computing r-orientation, r-eye, and r-sf for each pair of neurons
pre_id = pair_ids(:,1);
post_id = pair_ids(:,2);

pair_positions = ccg_output.data.unit_info.y(pair_ids);
ccg_output.pair_distance  = abs(pair_positions(:,1)-pair_positions(:,2));

ccg_output.pre_id = pre_id;
ccg_output.post_id = post_id;

% ccg_exclusions
noise_distribution2 = [ccg_output.ccg_control(:, 1:50), ccg_output.ccg_control(:, 151:201)];
ccg_output.area = nansum(ccgs,2);
ccg_output.noise_std2 = std(noise_distribution2,[],2, 'omitnan');
ccg_output.noise_mean2 = nanmean(noise_distribution2,2);

% ccg peak and trough width
ccg_output.peak_width = zeros(size(ccgs,1),1);
ccg_output.trough_width = zeros(size(ccgs, 1),1);

for i = 1:size(ccgs, 1)
    ccg = ccgs(i,:);
    % peak width
    ccg_left = ccg(1:peak_lag(i));
    
    ccg_right = ccg(peak_lag(i):end);
    left_ind = find((ccg_left<peaks(i)-abs(peaks(i)/2)), 1, 'last')+1;
    right_ind = find((ccg_right<peaks(i)-abs(peaks(i)/2)), 1, 'first')+(peak_lag(i)-1)-1;
    if isempty(left_ind)
        left_ind = 1;
    end
    if isempty(right_ind)
        right_ind = maxlag+1;
    end
    ccg_output.peak_width(i) = right_ind-left_ind;
    
    % trough width
    ccg_left = ccg(1:trough_lag(i));
    left_ind = find((ccg_left>troughs(i)+abs(troughs(i)/2)), 1, 'last')+1;
            
    ccg_right = ccg(trough_lag(i):end);
    right_ind = find((ccg_right>troughs(i) + abs(troughs(i)/2)), 1, 'first')-1+(trough_lag(i)-1);
    
    if isempty(left_ind)
        left_ind = 1;
    end
    if isempty(right_ind)
        right_ind = maxlag+1;
    end
    ccg_output.trough_width(i) = right_ind-left_ind;
end

if params.large_output
else
    ccg_output.ccg_control = [];
end
ccg_output.ccgs = ccgs;

% Get significant CCGs
[ccg_sig, sig_idx] = get_significant_ccgs(ccg_output, params);
disp("total number of pairs: " + length(sig_idx));
disp("number of significant pairs: " + sum(sig_idx));

%% Compute signal correlation vs CCG
RF_mat_selected = RF_mat(:,:,selected_ids_sorted);
n_sig = 0; n_all = 0;
for i=1:(size(ccg_output.ccgs,1)/2)
    i_pair = i*2-1;
    if any(selected_ids==ccg_output.data.unit_info.cid(ccg_output.pre_id(i_pair))+1) && any(selected_ids==ccg_output.data.unit_info.cid(ccg_output.post_id(i_pair))+1)
        RF_pre = RF_mat_selected(:,:,ccg_output.pre_id(i_pair));
        RF_post = RF_mat_selected(:,:,ccg_output.post_id(i_pair));
        corr = corrcoef(RF_pre(:), RF_post(:));
        n_all = n_all + 1;
        if any(ccg_sig.pre_id==ccg_output.pre_id(i_pair) & ccg_sig.post_id==ccg_output.post_id(i_pair))
            n_sig = n_sig + 1;
            ccg_result.significant_signal_corr = [ccg_result.significant_signal_corr; corr(1,2)];
            ccg_result.significant_peak_efficacy = [ccg_result.significant_peak_efficacy; ccg_output.peaks(i_pair)];
            ccg_result.significant_peak_lag = [ccg_result.significant_peak_lag; ccg_output.peak_lag(i_pair)];
            ccg_result.significant_bankid = [ccg_result.significant_bankid; 0];
        else            
            ccg_result.nonsignificant_signal_corr = [ccg_result.nonsignificant_signal_corr; corr(1,2)];
        end
    end
end

%% Plot signal correlation vs CCG
figure;
ax1 = subplot(9,1,1:4);
hold on
X_non = ccg_result.nonsignificant_signal_corr;
histogram(X_non, 'FaceColor', [0.5, 0.5, 0.5], 'DisplayName',['Non-significant pairs', newline, '(N=',num2str(length(X_non)),')'])
ylabel('N(non-significant pairs)')
yyaxis right
ylabel('N(significant pairs)')
set(gca, 'ylim', [0, 150])
X_syn = ccg_result.significant_signal_corr(abs(ccg_result.significant_peak_lag)<=1);
histogram(X_syn, 50, 'FaceColor', [0.8500 0.3250 0.0980], 'DisplayName',['Significant pairs', newline, '(synchronous, N=',num2str(length(X_syn)),')'])
X_asyn = ccg_result.significant_signal_corr(abs(ccg_result.significant_peak_lag)>1);
histogram(X_asyn, 50, 'FaceColor', [0 0.4470 0.7410], 'DisplayName',['Significant pairs', newline, '(asynchronous, N=',num2str(length(X_asyn)),')'])
legend()
set(gca, 'xlim', [-0.4, 1])


ax2 = subplot(9,1,6:9);
hold on
X1 = ccg_result.significant_signal_corr(abs(ccg_result.significant_peak_lag)<=1)'; y1 = ccg_result.significant_peak_efficacy(abs(ccg_result.significant_peak_lag)<=1)';
[b, S1] = polyfit(X1,y1,1);
p1 = scatter(X1, y1, 10, 'o', 'MarkerEdgeColor', [0.8500 0.3250 0.0980], 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Synchronous pairs');
plot([-0.5,1],[-0.5*b(1)+b(2),sum(b)],'-', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5)
X2 = ccg_result.significant_signal_corr(abs(ccg_result.significant_peak_lag)>1)'; y2 =ccg_result.significant_peak_efficacy(abs(ccg_result.significant_peak_lag)>1)';
[b, S2] = polyfit(X2,y2,1);
p2 = scatter(X2, y2, 10, 'o', 'MarkerEdgeColor', [0 0.4470 0.7410], 'MarkerEdgeAlpha', 0.5, 'DisplayName', 'Asynchronous pairs');
plot([-0.5,1],[-0.5*b(1)+b(2),sum(b)],'-', 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5)
legend([p1, p2],'Location','northwest')
set(gca, 'xlim', [-0.4, 1])
set(gca, 'ylim', [0, 0.1])
ylabel('CCG peak')

linkaxes([ax1,ax2],'x')
xlabel('Signal correlation')