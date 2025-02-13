classdef CurrentSourceDensityAnalysis < handle
    % here we use the sourceDatsets LF band data to compute CSDs aligned to each trial data trial

    properties(Transient)
        imec
        tvec_lfp (:, 1) double

        ss % raw unfiltered snippets
        ss_filt % good trials, filtered LFP interpolated over bad channels
        tvec_lfp_filt (:, 1) double
    end

    properties
        datasetName string
        tsi (:, 1) Neuropixel.TrialSegmentationInfo % in lfp sampling rate
        conductivity_S_per_m = 0.4; % conductivity of primate gray matter, in Siemens / m, value from Logothetis et al., 2007

%         column_via_kriging = true;
        window_align_ms =  [-1000 6500];

        erp_smooth_um = 900;
        csd_smooth_um = 450;
        
        highpass_corner_freq = [];
        lowpass_corner_freq = 25;
        butter_filter_order = 5; % will be doubled effectively by filtfilt

        filt_window_padding_ms = [50 50]; % in samples
        y_spacing_erp (1, 1) = 10;

        mask_good_ch (:, 1) logical
        mask_good_trials (:, 1) logical

        % Combined into lfp channel group averages
        tvec_erp (:, 1) single
        ypos_erp (:, 1) single % nChGroups (um), yposition for erp variables below
        erp (:, :) single % nChGroups x time (uV), smoothed version of erp_chGroupAvg over channels
        erp_current (:, :) single  % nChGroups x time (nA/mm^2)
        erp_csd (:, :) single % nChGroups x time (nA/mm^3)
    end

    properties(Dependent)
        nChGroups
        nGoodTrials
        nGoodChannels

        yspacing_erp
        erp_csd_zscore
    end

    methods
        function csd = CurrentSourceDensityAnalysis(tsi, imec, varargin)
            p = inputParser();
            p.parse(varargin{:});
            
            if nargin < 1
                return;
            end

            if tsi.fs ~= imec.fsLF
                % tsi built from AP sync line, convert to LF sample rate
                tsi = tsi.convertToDifferentSampleRate(imec.fsLF);
            end

            csd.tsi = tsi;
            csd.imec = imec;
            assert(csd.imec.hasLF, 'ImecDataset must have LF datasets');
        end

        function compute_all(csd)
            csd.extract_lfp_snippets();
            csd.filter_lfp();
            csd.compute_erp();
            csd.compute_csd();
        end

        function n = get.nGoodTrials(csd)
            n = nnz(csd.mask_good_trials);
        end

        function n = get.nGoodChannels(csd)
            n = nnz(csd.mask_good_ch);
        end
        
        function dy = get.yspacing_erp(csd)
            dy = abs(csd.ypos_erp(2) - csd.ypos_erp(1));
        end
        
        function erp_csd_zsc = get.erp_csd_zscore(csd)
            % zscore by baseline at each height on probe
            mask_pre = csd.tvec_erp < -200;
            erp_csd_zsc = (csd.erp_csd - mean(csd.erp_csd(:, mask_pre), 'all')) ./ std(csd.erp_csd(:, mask_pre), 0, 'all');
        end

        function ss = extract_windows(csd, window_ms, chIds)
            starts = makecol(csd.tsi.idxStart);
            window = round(window_ms * csd.imec.fsLF / 1000);

            % we extract lfp snippets from source datasets
            ss = csd.imec.readLFSnippetSet(starts, window, 'channel_ids', chIds, ...
                 'center', false, 'car', false, 'fromSourceDatasets', false);
        end            

        function extract_lfp_snippets(csd)
            imec = csd.imec; %#ok<*PROP,*PROPLC> 
            chmap = imec.channelMap;
            chIds = chmap.channelIdsMapped;
            
            windowMs = csd.window_align_ms;

            % extract baseline data
            windowBaselineMs = [windowMs(1)-300, 0];
            ss_baseline = csd.extract_windows(windowBaselineMs, chIds);
            baseline = ss_baseline.data;
            baseline_level = int16(median(baseline, 2)); %, 'native'); % C x 1 x tRials

            % extract LFP windows for move aligned trials
            windowMs_padded = [windowMs(1) + csd.filt_window_padding_ms(1), windowMs(2) + csd.filt_window_padding_ms(2)];
            csd.ss = csd.extract_windows(windowMs_padded, chIds);

            % subtract off baseline
            csd.ss.data = csd.ss.data - baseline_level;
             
            % compute time vector
            dt_ms = 1000 / imec.fsLF;
            csd.tvec_lfp = lindelta(windowMs(1), dt_ms, size(csd.ss.data, 2))';
        end

        function filter_lfp(csd, args)
            arguments
                csd
                args.show_figures = false;
            end

            % filters in time, detects bad trials and channels, interpolates over bad channels
            imec = csd.imec;
            chmap = imec.channelMap;
            chIds = chmap.channelIdsMapped;
            nCh = numel(chIds);

            assert(~isempty(csd.ss), 'Run extract_lfp_snippets()');
            lfp = csd.ss.data;  

            lfp = filterLFPInTime(lfp, imec.fsLF);

            % time shift in frequency domain to account for channels-within-ADC sampling offsets
            sampleShifts =  imec.channelMap.adcSampleShift;
            lfp = NHPPixel.fdomain_timeshift(lfp, sampleShifts);

            % slice off filter padding windows
            dt_ms = 1000 / imec.fsLF;
            filt_window_samples = floor(csd.filt_window_padding_ms ./ dt_ms);
            maskTime = true(size(lfp, 2), 1);
            maskTime(1:filt_window_samples(1)) = false;
            maskTime(end-filt_window_samples(2)+1:end) = false;
            lfp = lfp(:, maskTime, :);

            % figure out which channels and trials are good
            [ch_mask, trials_mask, ~] = determineGoodChannelsTrials(lfp, imec.goodChannelInds);
            trials_mask = squeeze(trials_mask);
            csd.mask_good_ch = ch_mask;
            csd.mask_good_trials = trials_mask;

            
            % interpolate across bad channels via kriging
            lfp = krigingInterpolationBadChannels(lfp, chmap.xcoords, chmap.ycoords, ch_mask);

            % interpolate bad channels? or just mask them
            lfp = lfp(:, :, trials_mask);
            
%             lfp = lfp - median(lfp, 1);

            

            % compute time vector
            dt_ms = 1000 / imec.fsLF;
            csd.tvec_lfp_filt = lindelta(csd.window_align_ms(1), dt_ms, nnz(maskTime))';

            % reconstruct cleaned snippet set
            maskChannels = csd.mask_good_ch;
            ss_filt = csd.ss.selectData(maskChannels=maskChannels, maskSnippets=csd.mask_good_trials, maskTime=maskTime);
            ss_filt.data = lfp;
            csd.ss_filt = ss_filt;
            
            function lfp = filterLFPInTime(lfp, fs)
                % C x T x R --> T x C*R
                temp = TensorUtils.reshapeByConcatenatingDims(double(lfp), {2, [1 3]});
                
                % lowpass (or bandpass filter
                if isempty(csd.highpass_corner_freq)
                    debug('Low-pass filtering LFP snippets\n');
                    Wn = csd.lowpass_corner_freq / (fs/2);
                    [z, p, k] = butter(csd.butter_filter_order, Wn, 'low');
                else
                    debug('Band-pass filtering LFP snippets\n');
                    Wn = [csd.highpass_corner_freq csd.lowpass_corner_freq] / (fs/2);
                    [z, p, k] = butter(csd.butter_filter_order, Wn, 'bandpass');
                end
                [sos, g] = zp2sos(z, p, k);
                temp = filtfilt(sos, g, temp);
                
                % T x C*R --> C x T x R
                lfp = TensorUtils.undoReshapeByConcatenatingDims(temp, {2, [1 3]}, size(lfp));
            end

            function [ch_mask, trials_mask, var_lfp] = determineGoodChannelsTrials(lfp, goodChannelInds)
                debug('Determining good channels and trials\n');
                ch_mask = ismember((1:nCh)', goodChannelInds);
                vlfp = var(single(lfp), 1, 2); % variance over time, ch x 1 x trials
                var_lfp = squeeze(vlfp);

                % determine channels based on their variance
                vlfp_by_ch = median(vlfp, 3);
                median_vlfp = median(vlfp_by_ch(ch_mask));
                mad_vlfp = mad(vlfp_by_ch(ch_mask));

                if args.show_figures
                    bad = setdiff(1:384, goodChannelInds);
                    if ~isempty(bad)
                        clf
                        stem(1:384, vlfp_by_ch, 'k-');
                        hold on;
                        yline(median_vlfp + 3*mad_vlfp);
                        yline(median_vlfp - 3*mad_vlfp);
                        stem(bad, vlfp_by_ch(bad), 'r-', 'filled');
                        hold off;
                    end
                end

                % remove channels whose variance is outside of this band
                nChPreVarThresh = nnz(ch_mask);
                ch_mask = ch_mask & abs(vlfp_by_ch - median_vlfp) <= 3*mad_vlfp;
                debug('Removing %d / %d channels based on variance\n', nChPreVarThresh - nnz(ch_mask), nChPreVarThresh);

                % determine trials based on their variance pattern over channels
                sWarn = warning('off', 'stats:pca:ColRankDefX');
                [~, score] = TensorUtils.pcaAlongDim(vlfp(ch_mask, :, :), 1, 'NumComponents', 3);
                warning(sWarn);
                
                trials_mask = any(abs(zscore(score, 0, 3)) < 6, 1);
                debug('Removing %d / %d trials based on variance pattern over channels\n', nnz(~trials_mask), numel(trials_mask));

                % and again based on each good channels variance over trials
                med_over_trials = median(vlfp(ch_mask, :, trials_mask), 3);
                mad_over_trials = mad(vlfp(ch_mask, :, trials_mask), 1, 3);
                bad_trials = mean(abs(vlfp(ch_mask, :, :) - med_over_trials) > 6 * mad_over_trials, 1) > 0.1;
                debug('Removing %d / %d trials based on channel-wise variance\n', nnz(bad_trials), nnz(trials_mask));
                trials_mask = trials_mask & ~bad_trials;

                % and now that everything is clean, check channels one last time
                % vlfp tends to look mostly like a rank 1 matrix with variation over channels, and then over trials
                % divide that out and look for outliers
%                 norm_vlfp = vlfp ./ median(vlfp, 1) ./ median(vlfp, 3);
%                 grand_med = median(norm_vlfp(ch_mask, :, trials_mask), 'all');
%                 grand_mad = mad(norm_vlfp(ch_mask, :, trials_mask), 1, 'all');
%                 norm_vlfp_score = (norm_vlfp - grand_med) ./ grand_mad;
%                 bad_ch = mean(abs(norm_vlfp_score) > 12, 3) > 0.05;
%                 debug('Removing %d / %d channels based on outlier variance again\n', nnz(bad_ch & ch_mask), nnz(ch_mask));
%                 ch_mask = ch_mask & ~bad_ch;
                
                if args.show_figures
                    clf;
                    pmat(vlfp(ch_mask, :, trials_mask));
                    caxis([0 500]);
                end
            end

            function lfp = krigingInterpolationBadChannels(lfp, xcoords, ycoords, ch_mask_good)
                debug('Interpolating over bad channels\n')
                % lfp is channels x time

                if all(ch_mask_good)
                    return;
                end

                xc_good = xcoords(ch_mask_good);
                yc_good = ycoords(ch_mask_good);

                xc_bad = xcoords(~ch_mask_good);
                yc_bad = ycoords(~ch_mask_good);

                % dists is bad x good
                kriging_distance_um = 20;
                p = 1.3;
                dists = sqrt((xc_bad - xc_good').^2 + (yc_bad - yc_good').^2);
                weights = exp(-(dists ./ kriging_distance_um) .^ p);
                weights(weights < 0.005) = 0;
                weights = weights ./ sum(weights, 2);
        
                % only take good channels with non-zero weights for the multiplication
                mask_good_take = any(weights > 0, 1);
                lfp(~ch_mask_good, :) = weights(:, mask_good_take) * lfp(mask_good_take, :); %  (chBad x chGood) * (chGood x time)
            end
        end

        function compute_erp(csd)
            % compute ERP by channel
            lfp = csd.ss_filt.data; % good ch x time x good trials
            tvec_full = csd.tvec_lfp_filt;
            mask_pre = tvec_full < 200;
            lfp = lfp - mean(lfp(:, mask_pre, :), 2);

            erp_by_channel = mean(lfp(:, :, :), 3) .* csd.imec.lfScaleToUv;
            erp_by_channel = erp_by_channel - mean(erp_by_channel(:, mask_pre, :), 2); 

            
            yspacing_erp = csd.y_spacing_erp;
            chmap = csd.imec.channelMap;

%             ch_mask = csd.mask_good_ch; % we interpolate over the bad channels inside 
            ch_mask = true(csd.imec.channelMap.nChannelsMapped, 1);
            
%             if csd.column_via_kriging
            debug('Kriging ERP values to linear depth samples\n');
            [erp_columnar, y_columnar] = krigingToColumn(erp_by_channel, chmap.xcoords(ch_mask), chmap.ycoords(ch_mask), yspacing_erp);
%             else
%                 debug('Grid-interpolating ERP values to linear depth samples\n');
%                 [erp_columnar, y_columnar] = interpolateToColumn(erp_by_channel, chmap.xcoords(ch_mask), chmap.ycoords(ch_mask), yspacing_erp);
%             end
    
            if chmap.invertChannelsY
                erp_columnar = flipud(erp_columnar);
                y_columnar = flipud(y_columnar);
            end
            
            smoothWidth = ceil(csd.erp_smooth_um / yspacing_erp);

            csd.erp = smoothdata(erp_columnar, 1, 'gaussian', smoothWidth) ;
            csd.tvec_erp = tvec_full;
            csd.ypos_erp = y_columnar;

            function [erp_columnar, yout] = krigingToColumn(erp, xcoords, ycoords, yspacing)
                yout = (min(ycoords) : yspacing : max(ycoords))';
                nChOut = numel(yout);
                xout = repmat(mean(unique(xcoords)), nChOut, 1);
                
                % weights is nCh_columnar x nCh_in
                kriging_distance_um = 20;
                p = 1.3;
                dists = sqrt((xout - xcoords').^2 + (yout - ycoords').^2);
                weights = exp(-(dists ./ kriging_distance_um) .^ p);
                weights(weights < 0.005) = 0;
                weights = weights ./ sum(weights, 2);

                erp_columnar = pagemtimes(weights, erp);
            end

            function [erp_columnar, yout] = interpolateToColumn(erp, xcoords, ycoords, yspacing)
                % the progress bar was for when lfp was over trials, now we pass in erp so it's all one trial
                % lfp is now all good data in form ch x time x trials
                % we want to interpolate to even yspacing values
                
                yout = (min(ycoords) : yspacing : max(ycoords))';
                nChOut = numel(yout);
                xout = repmat(mean(unique(xcoords)), nChOut, 1);
                
                [~, nTime, nTrials] = size(erp);
                erp_columnar = nan(nChOut, nTime, nTrials);
                
%                 prog = ProgressBar(nTrials, 'Interpolating LFP to column');
%                 prog.enableParallel();
                for r = 1:nTrials
%                     for t = 1:nTime
%                         lfp_columnar(:, t, r) = griddata(xcoords, ycoords, lfp(:, t, r), xout, yout, 'cubic'); %#ok<GRIDD>
%                     end
%                     prog.update(r);
                    F = scatteredInterpolant(xcoords, ycoords, zeros(numel(xcoords), 1), 'natural');
                    for t = 1:nTime
                         F.Values = erp(:, t, r);
                        erp_columnar(:, t, r) = F(xout, yout);
                    end
%                     prog.update(r); %#ok<PFBNS>
                end
%                 prog.finish();
                
                % edge channels often come out nan
                if isnan(erp_columnar(1, 1, 1))
                    erp_columnar(1, :, :) = erp_columnar(2, :, :);
                end
                if isnan(erp_columnar(end, 1, 1))
                    erp_columnar(end, :, :) = erp_columnar(end-1, :, :);
                end
            end
        end

    
%             mask_pre = csd.tvec_lfp_filt < -200;
%             lfp_col = lfp_col - mean(lfp_col(:, mask_pre, :), 2);

%             
% 
%             erp = mean(lfp, 3) .* csd.imec.lfScaleToUv;
%             yspacing_erp = 10;

            
%                     debug('Smoothing down each column\n');
%                     for pad = 1:4
%                         this_pad_ch = pad:4:size(lfp, 1);
%                         this_pad_ch_good = ch_mask(this_pad_ch);
%                         this_pad_ch = this_pad_ch(this_pad_ch_good);
%                         lfp(this_pad_ch, :, :) = smoothdata(lfp(this_pad_ch, :, :), 1, 'gaussian', 10, 'SamplePoints', find(this_pad_ch_good));
%                     end

%                     erp_by_channel = mean(lfp(ch_mask, :, trials_mask), 3) .* imec.lfScaleToUv;
%                     erp_by_channel = erp_by_channel - mean(erp_by_channel(:, mask_pre, :), 2); 
%                 
%                     debug('Grid-interpolating LFP values to linear depth samples\n');
%                     [erp_columnar, y_columnar] = interpolateLFPToColumn(erp_by_channel, chmap.xcoords(ch_mask), chmap.ycoords(ch_mask), yspacing_erp);
%             
%                     if chmap.invertChannelsY
%                         erp_columnar = flipud(erp_columnar);
%                         y_columnar = flipud(y_columnar);
%                     end
%                     
%                     smoothWidth = ceil(csd.erp_smooth_um / yspacing_erp);
% 
%                     csd.erp = smoothdata(erp_columnar, 1, 'gaussian', smoothWidth) ;
%                     csd.tvec_erp = tvec_full;
%                     csd.ypos_erp = y_columnar;
% 
%             erp_columnar = mean(lfp_col, 3) .* csd.imec.lfScaleToUv;
%             erp_columnar = erp_columnar - mean(erp_columnar(:, mask_pre, :), 2); 
%                 
%             if chmap.invertChannelsY
%                 erp_columnar = flipud(erp_columnar);
%                 yc_col = flipud(yc_col);
%             end
%             
%             % smooth the columnar data in Y
%             smoothWidth = ceil(csd.erp_smooth_um / csd.y_spacing_erp);
%             csd.erp = smoothdata(erp_columnar, 1, 'gaussian', smoothWidth) ;
%             csd.tvec_erp = csd.tvec_lfp_filt;
%             csd.ypos_erp = yc_col;
% 

        function compute_csd(csd)
            %% Compute CSD from ERP
            % output of savitzkyGolayFilt needs to be inverted for CSD (negative second derivative), but is scaled to mV / mm^2
            % to scale_to_nA/mm^3, we multiply by the negative conductivity in nA / (mV*mm)
            spacing_mm = csd.y_spacing_erp / 1000;
            scale_uV_to_mV = 1/1000;
            conductivity_nA_per_mV_mm = csd.conductivity_S_per_m * 1000;

            frameSize = ceil(csd.csd_smooth_um / csd.y_spacing_erp);
            if round(frameSize / 2) == frameSize/2 % ensure odd
                frameSize = frameSize + 1;
            end
            
            % replicate edges so that CSD fades away at edges
            padding = frameSize;
            erp_padded = padarray(csd.erp, [padding 0 0], 'replicate', 'both');
            
            % current is first derivative
            erp_current = conductivity_nA_per_mV_mm .* scale_uV_to_mV .* TrialDataUtilities.Data.savitzkyGolayFilt(erp_padded, 'differentiationOrder', 1, 'frameSize', frameSize, 'dim', 1, 'samplingInterval', spacing_mm);
            csd.erp_current = erp_current(padding+1:end-padding, :, :);
            
            % csd is neg second derivative
            erp_csd = -conductivity_nA_per_mV_mm .* scale_uV_to_mV .* TrialDataUtilities.Data.savitzkyGolayFilt(erp_padded, 'differentiationOrder', 2, 'frameSize', frameSize, 'dim', 1, 'samplingInterval', spacing_mm);
            csd.erp_csd = erp_csd(padding+1:end-padding, :, :);
        end

        function extract_erp_and_compute_csd(csd, args)
            arguments
                csd
                args.show_figures = false;
%                 args.windowMs = [-1000 8500];
            end
%             show_figures = args.show_figures;

            imec = csd.imec; %#ok<*PROPLC> 
            chmap = imec.channelMap;
            chIds = chmap.channelIdsMapped;
            nCh = numel(chIds);
            
            windowMs = csd.window_align_ms;
            
            % extract baseline data
            windowBaselineMs = [windowMs(1)-300, 0];
            ss_baseline = csd.extract_windows(windowBaselineMs, chIds);
            baseline = ss_baseline.data;
            baseline_level = int16(median(baseline, 2)); %, 'native'); % C x 1 x tRials

%             % extract LFP windows for move aligned trials
%             csd.ss = csd.extract_windows(windowMs, chIds);
% 
%             % subtract off baseline
%             csd.ss.data = csd.ss.data - baseline_level;
%             lfp = csd.ss.data;            
%             
%             debug('Low-pass filtering LFP snippets\n');
%             lfp = filterLFP(lfp, imec.fsLF);
            
            % compute time vector
            dt_ms = 1000 / imec.fsLF;
            tvec_full = lindelta(windowMs(1), dt_ms, csd.ss_filt.nTimepoints)';

            lfp = csd.ss_filt.data;
            %% convert grid sampled LFP to linear column and average over trials
            
            strategy = "grid_interpolate";
            switch strategy
                case "grid_interpolate"
                    yspacing_erp = 10;

                    mask_pre = tvec_full < 200;
                    lfp = lfp - mean(lfp(:, mask_pre, :), 2);

%                     debug('Smoothing down each column\n');
%                     for pad = 1:4
%                         this_pad_ch = pad:4:size(lfp, 1);
%                         this_pad_ch_good = ch_mask(this_pad_ch);
%                         this_pad_ch = this_pad_ch(this_pad_ch_good);
%                         lfp(this_pad_ch, :, :) = smoothdata(lfp(this_pad_ch, :, :), 1, 'gaussian', 10, 'SamplePoints', find(this_pad_ch_good));
%                     end

                    erp_by_channel = mean(lfp, 3) .* imec.lfScaleToUv;
                    erp_by_channel = erp_by_channel - mean(erp_by_channel(:, mask_pre, :), 2); 
                
                    debug('Grid-interpolating LFP values to linear depth samples\n');
                    [erp_columnar, y_columnar] = interpolateLFPToColumn(erp_by_channel, chmap.xcoords(csd.mask_good_ch), chmap.ycoords(csd.mask_good_ch), yspacing_erp);
            
                    if chmap.invertChannelsY
                        erp_columnar = flipud(erp_columnar);
                        y_columnar = flipud(y_columnar);
                    end
                    
                    smoothWidth = ceil(csd.erp_smooth_um / yspacing_erp);

                    csd.erp = smoothdata(erp_columnar, 1, 'gaussian', smoothWidth) ;
                    csd.tvec_erp = tvec_full;
                    csd.ypos_erp = y_columnar;
                    
                case "group_average" % this is the old method that is pretty much the same but less elegant
                    
                    % average over channel groups of 4 (old strategy)
                    debug('Performing averaging over groups of channels\n');
                    [W, ypos_by_group, mask_group_inpainted] = computeGroupAveragingWeights(chmap, ch_mask);
                    
                    erp_by_channel = mean(lfp(ch_mask, :, trials_mask), 3) .* imec.lfScaleToUv;
                    
                    mask_pre = tvec_full < -400;
                    erp_by_channel = erp_by_channel - mean(erp_by_channel(:, mask_pre, :), 2); 
                    
                    erp_by_group = W' * erp_by_channel;

                    if any(mask_group_inpainted)
                        erp_by_group(mask_group_inpainted, :) = NaN;
                        erp_by_group = fillmissing(erp_by_group, 'pchip', 1);
                    end 

                   yspacing_erp = double(abs(median(diff(ypos_by_group)))); % expected to be 60
                   smoothWidth = ceil(csd.erp_smooth_um / yspacing_erp); % expected to be 19
                   csd.erp = smoothdata(erp_by_group, 1, 'gaussian', smoothWidth);
                   csd.tvec_erp = tvec_full;
                   csd.ypos_erp = ypos_by_group;
                   
                otherwise
                    error('Unknown interpolation strategy');
            end

            %% Compute CSD from ERP
            % output of savitzkyGolayFilt needs to be inverted for CSD (negative second derivative), but is scaled to mV / mm^2
            % to scale_to_nA/mm^3, we multiply by the negative conductivity in nA / (mV*mm)
            spacing_mm = yspacing_erp / 1000;
            scale_uV_to_mV = 1/1000;
            conductivity_nA_per_mV_mm = csd.conductivity_S_per_m * 1000;

            frameSize = ceil(csd.csd_smooth_um / yspacing_erp);
            if round(frameSize / 2) == frameSize/2 % ensure odd
                frameSize = frameSize + 1;
            end
            
            % replicate edges so that CSD fades away at edges
            padding = frameSize;
            erp_padded = padarray(csd.erp, [padding 0 0], 'replicate', 'both');
            
            % current is first derivative
            erp_current = conductivity_nA_per_mV_mm .* scale_uV_to_mV .* TrialDataUtilities.Data.savitzkyGolayFilt(erp_padded, 'differentiationOrder', 1, 'frameSize', frameSize, 'dim', 1, 'samplingInterval', spacing_mm);
            csd.erp_current = erp_current(padding+1:end-padding, :, :);
            
            % csd is neg second derivative
            erp_csd = -conductivity_nA_per_mV_mm .* scale_uV_to_mV .* TrialDataUtilities.Data.savitzkyGolayFilt(erp_padded, 'differentiationOrder', 2, 'frameSize', frameSize, 'dim', 1, 'samplingInterval', spacing_mm);
            csd.erp_csd = erp_csd(padding+1:end-padding, :, :);
            
            csd.pmatbal_csd();


            % old style
%             function lfp = filterLFP(lfp, fs) 
%                 % C x T x R --> T x C*R
%                 temp = TensorUtils.reshapeByConcatenatingDims(double(lfp), {2, [1 3]});
%                 dfilt = designfilt('lowpassiir', 'SampleRate', fs, 'FilterOrder', 10, 'PassbandFrequency', csd.lowpass_corner_freq, 'PassbandRipple', 0.01, 'StopbandAttenuation', 80); %#ok<*PROPLC>
%                 temp = filtfilt(dfilt, temp);
%                 lfp = TensorUtils.undoReshapeByConcatenatingDims(temp, {2, [1 3]}, size(lfp));
%             end
%             
            
            
            function [lfp_columnar, yout] = interpolateLFPToColumn(lfp, xcoords, ycoords, yspacing)
                % the progress bar was for when lfp was over trials, now we pass in erp so it's all one trial
                % lfp is now all good data in form ch x time x trials
                % we want to interpolate to even yspacing values
                
                yout = (min(ycoords) : yspacing : max(ycoords))';
                nChOut = numel(yout);
                xout = repmat(mean(unique(xcoords)), nChOut, 1);
                
                [~, nTime, nTrials] = size(lfp);
                lfp_columnar = nan(nChOut, nTime, nTrials);
                
%                 prog = ProgressBar(nTrials, 'Interpolating LFP to column');
%                 prog.enableParallel();
                for r = 1:nTrials
%                     for t = 1:nTime
%                         lfp_columnar(:, t, r) = griddata(xcoords, ycoords, lfp(:, t, r), xout, yout, 'cubic'); %#ok<GRIDD>
%                     end
%                     prog.update(r);
                    F = scatteredInterpolant(xcoords, ycoords, zeros(numel(xcoords), 1), 'natural');
                    for t = 1:nTime
                         F.Values = lfp(:, t, r);
                        lfp_columnar(:, t, r) = F(xout, yout);
                    end
%                     prog.update(r); %#ok<PFBNS>
                end
%                 prog.finish();
                
                % edge channels often come out nan
                if isnan(lfp_columnar(1, 1, 1))
                    lfp_columnar(1, :, :) = lfp_columnar(2, :, :);
                end
                if isnan(lfp_columnar(end, 1, 1))
                    lfp_columnar(end, :, :) = lfp_columnar(end-1, :, :);
                end
            end
            
            function [W, ypos_by_group, mask_group_inpainted] = computeGroupAveragingWeights(chmap, mask_good_ch)
                % erp should already be only good channels
                ch_y = chmap.ycoords;
                ch_row = ch_y / chmap.yspacing;
                
                % we will end up only keeping the fully populated rows accomodating drift, but we track all of them
                row_inds_with_offsets = min(ch_row) : max(ch_row);
                nRowsWithOffset = numel(row_inds_with_offsets);
                nRowsPerGroup = 2;
                nGroups = floor(nRowsWithOffset / nRowsPerGroup);

                nCh = numel(chmap.channelIdsMapped);
                % compute the weighting matrix that will do drift correciton and groups of 4 averaging
                which_group = nan(nCh, 'single');
                W = zeros(nCh, nGroups, 'single');
                ypos_by_group = zeros(nGroups, 1, 'single');
                for iG = 1:nGroups
                    row_inds_this_group = row_inds_with_offsets((iG-1)*nRowsPerGroup + (1:nRowsPerGroup))';
                    ch_rows_with_current_offset = ch_row;
                    mask_all_this_group = ismember(ch_rows_with_current_offset , row_inds_this_group);
                    mask_this_group = mask_all_this_group & mask_good_ch;
                    nGood = nnz(mask_this_group);
                    W(mask_this_group, iG) = 1 / nGood;
                    which_group(mask_this_group) = iG;

                    % use all channels to comute the y position so that they are evenly spaced
                    ypos_by_group(iG) = mean(ch_y(mask_all_this_group));
                end

                % flip weighting matrix so that averaged data comes out with correct orientation
                if chmap.invertChannelsY
                    ypos_by_group = flipud(ypos_by_group);
                    W = fliplr(W);
                end

                % use only middle band of averaged-groups which include channels at all times
                channelsPerGroup = sum(W > 0, 1)';
                maxChannelsPerGroup = max(channelsPerGroup, [], 'all');
                mask_group_has_min_channels = all(sum(W > 0, 1) >= maxChannelsPerGroup / 2, 3);
                mask_groups = false(nGroups, 1);
                mask_groups(find(mask_group_has_min_channels, 1, 'first') : find(mask_group_has_min_channels, 1, 'last')) = true;

                mask_group_inpainted = mask_groups & channelsPerGroup == 0;
                if any(mask_group_inpainted)
                    warning('%d / %d channel groups have no valid channels, using interpolating infill', nnz(mask_group_inpainted), nnz(mask_groups));
                end

                % how to weight across channels into groups at each time
                W = W(mask_good_ch, mask_groups);
                mask_group_inpainted = mask_group_inpainted(mask_groups);
                ypos_by_group = ypos_by_group(mask_groups);
            end
        end
        
    end

    methods % Plotting
        function h = pmat_var_lfp(csd)
            h = pmat(csd.var_lfp(csd.mask_good_ch, csd.mask_good_trials));
        end

        function h = ptstack_lfp(csd, trialInd, varargin)
            if nargin < 2
                trialInd = 1;
            end

            h = ptstack(2, 1, csd.tvec_erp, csd.lfp_by_trial(csd.mask_good_ch, :, trialInd), 'gain', 4, 'normalize', false, varargin{:});
        end

        function h = ptstack_lfp_filt(csd, trialInd, varargin)
            if nargin < 2
                trialInd = 1;
            end

            h = ptstack(2, 1, csd.tvec_erp, csd.lfp_filt_by_trial(csd.mask_good_ch, :, trialInd), 'gain', 4, 'normalize', false, varargin{:});
        end

        function h = ptstack_lfp_filt_compare(csd, trialInd, varargin)
            if nargin < 2
                trialInd = 1;
            end

            cmap = [0 0 0; 0.898 0.447 0.690]; % black, blue;
            data = cat(3, csd.lfp_by_trial(csd.mask_good_ch, :, trialInd), csd.lfp_filt_by_trial(csd.mask_good_ch, :, trialInd));
            h = ptstack(2, 1, csd.tvec_erp, data, 'gain', 4, 'normalize', false, 'colormap', cmap, varargin{:});
        end

        function pmatbal_erp(csd)
            pmatbal(csd.erp, 'x', csd.tvec_erp, 'y', csd.ypos_erp, 'colorAxisLabel', '{\mu}V', 'colormap', 'default');
            xline(0);
            axis xy;
            set(gca, 'XAxisLocation', 'bottom');
            box off;
            xlabel('Time from Move (ms)');
            ylabel('Probe position (um)');
        end

        function ptstack_erp(csd)
            ptstack(2, 1, csd.tvec_erp, csd.erp, 'gain', 10, 'normalize', false);
        end

        function pmatbal_current(csd)
            pmatbal(csd.erp_current, 'x', csd.tvec_erp, 'y', csd.ypos_erp, 'colorAxisLabel', 'nA/mm^2', 'colormap', 'default');
            xline(0);
            axis xy;
            set(gca, 'XAxisLocation', 'bottom');
            box off;
            shading interp
            xlabel('Time from Move (ms)');
            ylabel('Probe position (um)');
        end

        function pmatbal_csd(csd)
            pmatbal(csd.erp_csd, 'x', csd.tvec_erp, 'y', csd.ypos_erp, 'colorAxisLabel', 'nA/mm^3', 'colormap', 'default');
            xline(0);
            axis xy;
            set(gca, 'XAxisLocation', 'bottom');
            box off;
            shading interp
            xlabel('Time from Move (ms)');
            ylabel('Probe position (um)');
        end
        
        function pmatbal_csd_zscore(csd)
            pmatbal(csd.erp_csd_zscore, 'x', csd.tvec_erp, 'y', csd.ypos_erp, 'colorAxisLabel', 'zscore', 'colormap', 'default');
            xline(0);
            axis xy;
            set(gca, 'XAxisLocation', 'bottom');
            box off;
            shading interp
            xlabel('Time from Move (ms)');
            ylabel('Probe position (um)');
        end
    end

    methods(Static) % simple manual construction
        function csd = build_csd_directly(erp_csd, tvec_erp, ypos_erp)
            csd = NeuropixelExpt.CurrentSourceDensityAnalysisRevised();
            csd.erp_csd = erp_csd;
            csd.tvec_erp = tvec_erp;
            csd.ypos_erp = ypos_erp;
        end
    end
    
    methods % locating specific sinks sources
        function [ypos, tpos, source_max, row, col] = find_superficial_source(csd, varargin)
            p = inputParser();
            p.addParameter('timeWindow', [-100 100], @isvector);
            p.addParameter('belowSuperficialSink', false, @(x) islogical(x) || isstringlike(x));
            p.addParameter('depthWindowRelSink', [0 800], @isvector);
            p.addParameter('depthWindowRelTop', [0 1000], @isvector);
            p.addParameter('showPlot', false, @islogical);
            p.parse(varargin{:});

            timeWindow = p.Results.timeWindow;
            depthWindowRelTop = p.Results.depthWindowRelTop;
            belowSuperficialSink = p.Results.belowSuperficialSink;

            mask_time = csd.tvec_erp >= timeWindow(1) & csd.tvec_erp <= timeWindow(2);
            top = max(csd.ypos_erp);
            depthFromY = @(y) -y + top;
            depthRelTop = depthFromY(csd.ypos_erp);
            mask_depth = depthRelTop >= depthWindowRelTop(1) & depthRelTop <= depthWindowRelTop(2);

            if belowSuperficialSink
                % shift top of depth window to lie below the superficial sink
                [ypos_sink, tpos_sink] = csd.find_superficial_sink('showPlot', false, 'belowSuperficialSource', false);
                depthSink = depthFromY(ypos_sink);
                depthWindowSink = depthSink + p.Results.depthWindowRelSink;
                mask_depth = mask_depth & depthRelTop >= depthWindowSink(1) & depthRelTop <= depthWindowSink(2);
            end
            
            strip = -inf(size(csd.erp_csd));
            strip(mask_depth, mask_time) = csd.erp_csd(mask_depth, mask_time);
            [source_max, ind] = max(strip, [], 'all', 'linear');
            [row, col] = ind2sub(size(strip), ind);

            ypos = csd.ypos_erp(row);
            tpos = csd.tvec_erp(col);

            if p.Results.showPlot
                clf;
                csd.pmatbal_csd();
                hold on;
                
                if belowSuperficialSink
                    plot(tpos_sink, ypos_sink, '+', 'Color', [0.2 0.2 0.2], 'MarkerSize', 25, 'LineWidth', 5);
                    plot(tpos_sink, ypos_sink, '+', 'Color', [0.4 1 0.4] , 'MarkerSize', 20, 'LineWidth', 1.5);
                end
                plot(tpos, ypos, '+', 'MarkerSize', 25, 'Color', [0.2 0.2 0.2], 'LineWidth', 5);
                plot(tpos, ypos, '+', 'MarkerSize', 20, 'Color', [0.4 1 0.4], 'LineWidth', 1.5);
                hold off;
            end
        end

        function [ypos, tpos, sink_min, row, col] = find_superficial_sink(csd, varargin)
            p = inputParser();
            p.addParameter('timeWindow', [-150 50], @isvector);
            p.addParameter('belowSuperficialSource', true, @(x) islogical(x) || isstringlike(x));
            p.addParameter('depthWindowRelSource', [0 800], @isvector);
            p.addParameter('depthWindowRelTop', [0 1500], @isvector);
            p.addParameter('showPlot', false, @islogical);
            p.parse(varargin{:});

            timeWindow = p.Results.timeWindow;
            depthWindowRelTop = p.Results.depthWindowRelTop;

            mask_time = csd.tvec_erp >= timeWindow(1) & csd.tvec_erp <= timeWindow(2);
            top = max(csd.ypos_erp);
            depthFromY = @(y) -y + top;
            depthRelTop = depthFromY(csd.ypos_erp);
            mask_depth = depthRelTop >= depthWindowRelTop(1) & depthRelTop <= depthWindowRelTop(2);

            belowSuperficialSource  = p.Results.belowSuperficialSource;
%             if isstringlike(belowSuperficialSource)
%                 belowSuperficialSource = string(belowSuperficialSource);
%                 switch belowSuperficialSource
%                     case 'auto'
%                         % do belowSuperficialSource if close to surface and large
%                         
                    
            if belowSuperficialSource
                % shift top of depth window to lie below the superficial source
                [ypos_source, tpos_source] = csd.find_superficial_source('showPlot', false);
                depthSource = depthFromY(ypos_source);
                depthWindowSource = depthSource + p.Results.depthWindowRelSource;
                mask_depth = mask_depth & depthRelTop >= depthWindowSource(1) & depthRelTop <= depthWindowSource(2);
            end

            strip = inf(size(csd.erp_csd));
            strip(mask_depth, mask_time) = csd.erp_csd(mask_depth, mask_time);
            [sink_min, ind] = min(strip, [], 'all', 'linear');
            [row, col] = ind2sub(size(strip), ind);

            ypos = csd.ypos_erp(row);
            tpos = csd.tvec_erp(col);

            if p.Results.showPlot
                clf;
                csd.pmatbal_csd();
                hold on;
                if belowSuperficialSource
                    plot(tpos_source, ypos_source, '+', 'Color', [0.2 0.2 0.2], 'MarkerSize', 25, 'LineWidth', 5);
                    plot(tpos_source, ypos_source, '+', 'Color', [0.4 1 0.4] , 'MarkerSize', 20, 'LineWidth', 1.5);
                end

                plot(tpos, ypos, '+', 'MarkerSize', 25, 'Color', [0.2 0.2 0.2], 'LineWidth', 5);
                plot(tpos, ypos, '+', 'MarkerSize', 20, 'Color', [0.4 1 0.4], 'LineWidth', 1.5);
                hold off;
            end
        end
        
        
        function [ypos, tpos, sink_min, row, col] = find_sink_near_reference(csd, ref, varargin)
            % find the largest sink within the vicinity of a reference position
            p = inputParser();
            p.addParameter('timeWindow', [-150 50], @isvector);
            p.addParameter('depthWindowRelRef', [-500 500], @isvector); % depth, + means deeper, - means superficial
            
            % we search above and below the reference for the sink, and by default take the minimum. However, if there 
            % are sinks below and above the ref, you can prefer the more superficial / deeper if it is below the value
            % preferBelowThresh, even if it is not the minimal value of the two.
            p.addParameter('preferSuperficial', false, @isvector);
            p.addParameter('preferDeep', false, @isvector);
            p.addParameter('preferBelowThresh', -30, @isvector);  % max value of sink to count if prefer* is active
            p.addParameter('showPlot', false, @islogical);
            p.parse(varargin{:});

            timeWindow = p.Results.timeWindow;
            depthWindowRelRef = p.Results.depthWindowRelRef;
            
            preferSuperficial = depthWindowRelRef(1) < 0 && p.Results.preferSuperficial;
            preferDeep = ~preferSuperficial && depthWindowRelRef(2) > 0 && p.Results.preferDeep;
            preferThresh = p.Results.preferBelowThresh;
            windowSuperficial = [depthWindowRelRef(1) 0];
            windowDeep = [0 depthWindowRelRef(2)];
            
            mask_time = csd.tvec_erp >= timeWindow(1) & csd.tvec_erp <= timeWindow(2);
            
            function [val, row, col] = findSink(mask_depth)
                strip = inf(size(csd.erp_csd));
                strip(mask_depth, mask_time) = csd.erp_csd(mask_depth, mask_time);
                [val, sink_ind] = min(strip, [], 'all', 'linear');
                [row, col] = ind2sub(size(strip), sink_ind);
            end

            depthFromRef = -csd.ypos_erp + ref;
            mask_superficial = depthFromRef >= windowSuperficial(1) & depthFromRef <= windowSuperficial(2);
            mask_deep = depthFromRef >= windowDeep(1) & depthFromRef <= windowDeep(2);
            [sink_s, row_s, col_s] = findSink(mask_superficial);
            [sink_d, row_d, col_d] = findSink(mask_deep);
            
            if preferSuperficial && sink_s <= preferThresh
                % take the superficial sink
                row = row_s;
                col = col_s;
                sink_min = sink_s;
            elseif preferDeep && sink_d <= preferThresh
                row = row_d;
                col = col_d;
                sink_min = sink_d;
            elseif sink_s <= sink_d
                row = row_s;
                col = col_s;
                sink_min = sink_s;
            else
                row = row_d;
                col = col_d;
                sink_min = sink_d;
            end
                
            ypos = csd.ypos_chGroups(row);
            tpos = csd.tvec_erp(col);

            if p.Results.showPlot
                clf;
                csd.pmatbal_csd();
                hold on;
                yline(ref, '-', 'Color', 'k');
                plot(tpos, ypos, '+', 'MarkerSize', 25, 'Color', [0.2 0.2 0.2], 'LineWidth', 5);
                plot(tpos, ypos, '+', 'MarkerSize', 20, 'Color', [0.4 1 0.4], 'LineWidth', 1.5);
                hold off;
            end
        end
    end

    methods(Static)
        function [offsets_um] = alignment_optimization(csdSet, varargin)
            p = inputParser();
            p.addParameter('origin', [], @(x) isvector(x) || isempty(x)); % initial offsets for each probe. if [] or where NaN, the superficial sink will be used
            p.addParameter('edgeRowSkip', csdSet(1).diff_smooth_nChGroups - 1, @isscalar);
            p.addParameter('maxShift', 300, @isscalar);
            p.addParameter('tWindow', [-200 200], @isscalar);
            p.addParameter('normalizeEach', true, @islogical);
            p.addParameter('extractWindowRelOrigin', [-1000 2250], @isvector);
            p.addParameter('showInitialAlignment', true, @islogical);
            p.addParameter('optimoptions', {}, @iscell);
            p.parse(varargin{:});
            showInitialAlignment = p.Results.showInitialAlignment;
            
            datasetNames = cat(1, csdSet.datasetName);

            tMin = p.Results.tWindow(1);
            tMax = p.Results.tWindow(2);
            normalizeEach = p.Results.normalizeEach;
            extractWindowRelOrigin = p.Results.extractWindowRelOrigin;

            depth_spacing = unique(arrayfun(@(csd) csd.yspacing_csd, csdSet));
            assert(numel(depth_spacing) == 1, 'CSDs have differing yspacing_csd');
            maxShiftRows = floor(p.Results.maxShift / depth_spacing);

            nC = numel(csdSet);
            skip = p.Results.edgeRowSkip;
            
            origin = p.Results.origin;
            if isempty(origin)
                origin = nan(nC, 1);
            end
            assert(numel(origin) == nC);

            % loop through and extract the matchable portion of each CSD, relative to the "origin" (e.g. the superficial sink)
            [row_origin, t_origin] = nanvec(nC);
            [ind_lim, depth_lim] = deal(nan(nC, 2));
            for iC = 1:nC
                csd = csdSet(iC);
                if isnan(origin(iC))
                    [origin(iC), t_origin(iC), ~, row_origin(iC)] = csd.find_superficial_sink('showPlot', false, 'belowSuperficialSource', true, 'depthWindowRelTop', [200 1000]);
                else
                    t_origin = 0;
                end
                
                depth_from_origin = -(csd.ypos_erp - origin(iC));

                ymask = depth_from_origin >= extractWindowRelOrigin(1) & depth_from_origin <= extractWindowRelOrigin(2);
                ymask(1:skip) = false;
                ymask(end-skip+1:end) = false;

                ind_lim(iC, 1) = find(ymask, 1, 'first');
                ind_lim(iC, 2) = find(ymask, 1, 'last');
                depth_lim(iC, 1) = depth_from_origin(ind_lim(iC, 1));
                depth_lim(iC, 2) = depth_from_origin(ind_lim(iC, 2));
            end

            min_depth = min(depth_lim(:, 1));
            max_depth = max(depth_lim(:, 2));
            row_depth = min_depth : depth_spacing : max_depth;
            nRows = numel(row_depth);

            tmask = csdSet(1).tvec_erp >= tMin & csdSet(1).tvec_erp <= tMax;
            nT = nnz(tmask);
            tvec_mask = csdSet(1).tvec_erp(tmask);
            col_sink_mask = nan(nC, 1);
            for iC = 1:nC
                [~, col_sink_mask(iC)] = min(abs(t_origin(iC) - tvec_mask));
            end

            data = nan(nRows, nT, nC, 'single');
            data_valid = false(nRows, nC);
            for iC = 1:nC
                take = ind_lim(iC, 1) : ind_lim(iC, 2);
                [~, insert_first] = min(abs(row_depth - depth_lim(iC, 1)));
                insert = insert_first : insert_first + numel(take) - 1;
                data(insert, :, iC) = csdSet(iC).erp_csd(take, tmask);
                data_valid(insert, iC) = true;

                if normalizeEach
                    normalizer = quantile(abs(data(insert, :, iC)), 0.98, 'all');
                    data(insert, :, iC) = data(insert, :, iC) ./ normalizer;
                end
            end
            
            if showInitialAlignment
                nanstrip = nan(nRows, 1, nC);
                dmat = TensorUtils.reshapeByConcatenatingDims(cat(2, data, nanstrip), {1, [2 3]});
                xvec = linspace(-0.3, 0.3, nT+1);
                xmat = xvec' + (0:nC-1);
                xvals = xmat(:);
                
                figure(1); clf;
                pmatbal(dmat, 'y', row_depth, 'x', xvals, 'colormap', 'reverse');
                
                hold on;
                xsink = xvec(col_sink_mask)' + (0:nC-1)';
                ysink = zeros(nC, 1);
                scatter(xsink, ysink, 100, 'k', 'filled', 'MarkerEdgeColor', 'w');
                
                set(gca, 'XAxisLocation', 'bottom', 'XTick', 0:nC-1, 'XTickLabels', datasetNames, 'XTickLabelRotation', 45);
                
                clim = quantile(abs(dmat(:)), 0.99);
                caxis([-clim clim]);
                hold off;
            end
            
            % pre-compute all pairwise distances so that the algorithm can proceed quickly
            % we want to compute costs(i, j, k) which is the penalty for mismatch between
            % csdSet(i) and csdSet(j) when j is shifted by shifts(k) down relative to i.
            shift_vals = -maxShiftRows : maxShiftRows;
            nShifts = numel(shift_vals);
            
            function D2 = distfun(ZI, ZJ)
                nn = sum(~isnan(ZI) & ~isnan(ZJ), 2);
                D2 = sum((ZI - ZJ).^2, 2, 'omitnan') ./ nn;
            end
            
            % first we want to generate shifted versions of all the csds in data
            prog = ProgressBar(nShifts, 'Precomputing costs of each pairwise CSD shift');
            shift_costs = nan(nC, nC, nShifts);
            for iS = 1:nShifts
                shift = shift_vals(iS);
                if shift >= 0
                    % session 2 above session 1
                    data_flat1 = reshape(data(1:(nRows-shift), :, :), [], nC)';
                    data_flat2 = reshape(data((shift+1):nRows, :, :), [], nC)';
                else
                    % session 1 above session 1
                    data_flat1 = reshape(data((-shift+1):nRows, :, :), [], nC)';
                    data_flat2 = reshape(data( 1:(nRows+shift), :, :), [], nC)';
                end

                shift_costs(:, :, iS) = pdist2(data_flat1, data_flat2, @distfun);
                prog.update(iS);
            end
            prog.finish();

            nvars = nC;
            lower_bound = repmat(double(fix(-maxShiftRows/2)), nC, 1);  % Lower bound
            upper_bound = repmat(double(fix( maxShiftRows/2)), nC, 1);  % Upper bound
            rng(1); % set random number state to predictable constant
            options = optimoptions('ga', 'PlotFcn', @plot_fun, 'UseParallel', false, 'PopulationSize', 50, 'Display', 'iter', 'MaxGenerations', 1000, p.Results.optimoptions{:});

            % 1:nC arg specifies that all nC vars are integer variables
            [row_offsets, ~] = ga(@obj_fun, nvars, [], [], [], [], lower_bound, upper_bound, [], 1:nC, options);
            row_offsets_um = row_offsets' * depth_spacing; 
            offsets_um = origin + row_offsets_um;

            function cost = obj_fun(shifts)
                % shifts is a vector of offsets with (positive being that this session should be translated upwards relative to the others)
                % nRows is over the population used in the ga
                nPop = size(shifts, 1);
                cost = zeros(nPop, 1, 'double');
                
                for iR = 1:nPop
                    [tf, relShiftInds] = ismember(shifts - shifts', shift_vals);
                    assert(all(tf(:)));
                    for iC1 = 1:nC
                        for iC2 = iC1+1:nC
                            cost(iR) = cost(iR) + shift_costs(iC1, iC2, relShiftInds(iC1, iC2));
                        end
                    end
                end
            end

%             function cost = obj_fun(shifts)
%                 % shifts is a vector of offsets with (positive being that this session should be translated upwards relative to the others)
%                 % nRows is over the population used in the ga
%                 nPop = size(shifts, 1);
%                 cost = zeros(nPop, 1, 'double');
%                 for iR = 1:nPop
%                     for iC = 1:nC
%                         for iC2 = iC1+1:nC
%                             offset = shifts(iR, iC2) - shifts(iR, iC1);
%                             if offset >= 0
%                                 % session 2 above session 1
%                                 ch1 = 1 : (nRows-offset);
%                                 ch2 = (offset+1) : nRows;
%                             else
%                                 % session 1 above session 2
%                                 ch2 = 1 : (nRows + offset);
%                                 ch1 = (-offset+1) : nRows;
%                             end
% 
%                             cost(iR) = cost(iR) + sum((data(ch1, :, iC1) - data(ch2, :, iC2)).^2, 'all', 'omitnan') ./ numel(ch1);
%                         end
%                     end
%                 end
%             end

            function state = plot_fun(options,state,flag) %#ok<INUSL>
                %GAPLOTBESTF Plots the best score and the mean score.
                %   STATE = GAPLOTBESTF(OPTIONS,STATE,FLAG) plots the best score as well
                %   as the mean of the scores.
                %
                %   Example:
                %    Create an options structure that will use GAPLOTBESTF
                %    as the plot function
                %     options = ptimoptions('ga','PlotFcn',@gaplotbestf);

                %   Copyright 2003-2016 The MathWorks, Inc.

                if size(state.Score,2) > 1
                    msg = getString(message('globaloptim:gaplotcommon:PlotFcnUnavailable','gaplotbestf'));
                    title(msg,'interp','none');
                    return;
                end

                persistent himg X0 Y0;
                switch flag
                    case 'init'
                        hold on;
                        [Y0, X0] = ndgrid(nRows:-1:1, linspace(-0.45, 0.45, nT));

                        himg = gobjects(nC, 1);
                        for iiC = 1:nC
                            X = X0 + (iiC-1);
                            himg(iiC) = pcolor(X, Y0, data(:, :, iiC));
                            himg(iiC).EdgeColor = 'none';
                        end

                        cmap = TrialDataUtilities.Color.cbrewer('div', 'RdBu', 256);
                        colormap(cmap);
                        cmax = quantile(abs(data), 0.99, 'all');
                        caxis([-cmax cmax]);
                        set(gca, 'XTick', 0:nC, 'XTickLabel', datasetNames, 'XTickLabelRotation', 45);
                        ylim([-maxShiftRows, nRows + maxShiftRows]);
                        xlim([-1 nC+1]);

                    case {'iter', 'done'}
                        [~, best] = min(state.Score);
                        offsets = state.Population(best, :);

                        for iiC = 1:nC
                            himg(iiC).YData = Y0 + offsets(iiC);
                        end
                        title(sprintf('Generation %d, Best cost %.3g', state.Generation, min(state.Best)));
                end
            end
        end

        function show_alignment_plot(csdSet, offsets_um, varargin)
            p = inputParser();
            p.addParameter('tWindow', [-200 200], @isscalar);
            p.addParameter('tSpacing', 20, @isscalar);
%             p.addParameter('skipEdges', false, @islogical);
            p.addParameter('showLabels', true, @islogical);
            p.addParameter('showGrandAverage', true, @islogical);
            p.addParameter('grandAverageFracDatasets', 0.8, @isscalar);
            p.addParameter('limits', [], @(x) true);
            p.addParameter('source', 'csd', @isstringlike);
            p.parse(varargin{:});

            tMin = p.Results.tWindow(1);
            tMax = p.Results.tWindow(2);
%             skipEdges = p.Results.skipEdges;
            showLabels = p.Results.showLabels;
            ypos = csdSet(1).ypos_chGroups;
            yspacing = csdSet(1).yspacing_chGroups;
            showGrandAverage = p.Results.showGrandAverage;
            grandAverageFracDatasets = p.Results.grandAverageFracDatasets;
            source = string(p.Results.source);

            tmask = csdSet(1).tvec_erp >= tMin & csdSet(1).tvec_erp <= tMax;
            tvec_plot = csdSet(1).tvec_erp(tmask);
            tSpacing = p.Results.tSpacing;
            nS = numel(csdSet);
            tOffsets = (tSpacing + tMax - tMin) * (0:nS-1)';

            nT = nnz(tmask);
%             if skipEdges
%                 skip = floor(csdSet(1).diff_smooth_nChGroups/2);
%             else
%                 skip = 0;
%             end
            nChGroups = csdSet(1).nChGroups - 2*skip;
            nS = numel(csdSet);

            data = zeros(nChGroups, nT, nS, 'single');
            for iS = 1:nS
                if source == "csd"
                    data(:, :, iS) = csdSet(iS).erp_csd(skip+1:end-skip, tmask, :);
                elseif source == "erp"
                    data(:, :, iS) = csdSet(iS).erp(skip+1:end-skip, tmask, :);
                else
                    error('Unkown source %s', source);
                end
            end
            ypos = ypos(skip+1:end-skip);

            ylabel = max(ypos) + yspacing;
            tlabel = (tMax + tMin) / 2;

            cla;

            if showGrandAverage
                grandAverage = mean(data, 3, 'omitnan');
                nDatasetsByRow = sum(~isnan(data(:, 1, :)), 3);
                thresh = nDatasetsByRow >= grandAverageFracDatasets * nS;
                tOffThis = -tOffsets(2);

                grandAverage(nDatasetsByRow < thresh, :, :) = NaN;
                NeuropixelExpt.CurrentSourceDensityAnalysis.pmat(grandAverage, 'x', tvec_plot + tOffThis, 'y', ypos, 'addColorbar', false);
                hold on

                plot(tOffThis, min(ypos) - 2*yspacing, '^', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

                if showLabels
                    text(tlabel + tOffThis, double(ylabel), "Average", 'FontWeight', 'bold', ...
                        'HorizontalAlign', 'center', 'VerticalAlign', 'bottom', 'YLimInclude', 'on');
                end
            end

            himg = gobjects(nS, 1);
            for iS = 1:nS
                himg(iS) = NeuropixelExpt.CurrentSourceDensityAnalysis.pmat(data(:, :, iS), 'x', tvec_plot + tOffsets(iS), 'y', ypos + offsets_um(iS), 'addColorbar', false);
                hold on;

                plot(tOffsets(iS), min(ypos) - 2*yspacing + offsets_um(iS), '^', 'MarkerSize', 5, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');

                if showLabels
                    text(tlabel + tOffsets(iS), double(ylabel + offsets_um(iS)), [ csdSet(iS).datasetName, sprintf("%+.0f {\\mu}m", offsets_um(iS)) ], ...
                        'HorizontalAlign', 'center', 'VerticalAlign', 'bottom', 'YLimInclude', 'on');
                end
            end

            cmap = flipud(TrialDataUtilities.Color.cbrewer('div', 'RdBu', 256));
            colormap(cmap);
            if isempty(p.Results.limits)
                cmax = quantile(abs(data), 0.999, 'all');
                caxis([-cmax cmax]);
            else
                limits = p.Results.limits;
                if isscalar(limits)
                    limits = [-abs(limits) abs(limits)];
                end
                caxis(limits);
            end
            set(gca, 'XTick', []);
            axis tight
            hold off;
            axis xy;
            shading interp;

            aa = AutoAxis.replaceScaleBars('xUnits', 'ms', 'yLength', 500, 'yUnits', '{\mu}m', 'xLength', 100);
            aa.addColorbar('location', AutoAxis.FullPositionSpec.outsideRightTop, 'labelCenter', 'nA/mm^3', 'labelAbove', 'source', 'labelBelow', 'sink');
            aa.axisMarginRight = 1.5;
            aa.axisMarginTop = 1;
            aa.update();
        end

        function [h, hcbar] = pmat(m, varargin)
            % visualize a matrix using pcolor

            p = inputParser();
            p.addParameter('x', [], @(x) isempty(x) ||  isvector(x));
            p.addParameter('y', [], @(x) isempty(x) || isvector(x));
            p.addParameter('addColorbar', true, @islogical);
            p.addParameter('colorAxisLabel', '', @isstringlike);
            p.parse(varargin{:});

            % cla;

            m = squeeze(m);

            if isvector(m)
                m = repmat(makerow(m), 2, 1);
            end

            if islogical(m)
                m = double(m);
            end

            if ndims(m) > 2 %#ok<ISMAT>
                warning('Selecting (:, :, 1) of tensor to display');
                m = m(:, :, 1);
            end

            % add an extra row onto m
            addRowCol = @(v) [v, v(:, end)+diff(v(:, end-1:end), 1, 2); ...
                v(end, :) + diff(v(end-1:end, :), 1, 1), 2*v(end, end)-v(end-1, end-1)];

            if isempty(p.Results.x)
                x = 0.5:size(m, 2)-0.5;
            else
                x = p.Results.x;
                dx = diff(x);
                x = x - dx([1:end end])/2;
            end
            if isempty(p.Results.y)
                y = 0.5:size(m, 1)-0.5;
            else
                y = p.Results.y;
                dy = diff(y);
                y = y - dy([1:end end])/2;
            end

            [X, Y] = meshgrid(x, y);

            % need an extra row and column because of the way that pcolor works
            m = addRowCol(m);
            X = addRowCol(X);
            Y = addRowCol(Y);

            h = pcolor(X,Y, m);

            set(h, 'EdgeColor', 'none');
            %colormap(parula);
            % colormap(flipud(cbrewer('div', 'RdYlBu', 256)));
            % colormap(pmkmp(256));
            %colormap gray;
            TrialDataUtilities.Color.cmocean('haline');

            if p.Results.addColorbar
                hcbar = colorbar;
                box(hcbar, 'off');
                set(hcbar, 'TickLength', 0);

                colorAxisLabel = string(p.Results.colorAxisLabel);
                if colorAxisLabel ~= ""
                    hcbar.YLabel.String = colorAxisLabel;
                end
            else
                hcbar = [];
            end

            box off
            axis ij
            axis on;

            axis on;
            set(gca, 'TickLength', [0 0], 'XAxisLocation', 'top');
            axis tight;
            box on;

        end

    end
end
