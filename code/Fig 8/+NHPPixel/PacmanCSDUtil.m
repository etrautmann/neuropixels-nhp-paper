classdef PacmanCSDUtil < handle

    properties
        sessinfo (:, 1) NHPPixel.PacmanSessionInfo

        %% load_trial_edges
        fsSG (1, 1) = 1000;
        fsNI (1, 1) double
        fsAP (1, 1) double

        % for all trials, before subselecting
        Tsync table
        imecTrialEdges_unfiltered (:, 2) uint64

        %% load_trial_alignment
        % after subselecting good trials
        trialInds (:, 1)
        condInd (:, 1) 
        task (:, 1) PacmanTaskCond
        
        imecTrialEdges (:, 2) uint64
        apIdx_Target (:, 1) uint64
        apIdx_TargetWithOptimalLag (:, 1) uint64
        apIdx_FirstPeak (:, 1) uint64
        apIdx_FirstPeakWithOptimalLag (:, 1) uint64

        csd (:, 1) NHPPixel.CurrentSourceDensityAnalysis
    end

    properties(Transient)
        imec (:, 1) Neuropixel.ImecDataset;
    end

    properties(Dependent)
        nCond
    end

    methods
        function pcu = PacmanCSDUtil(sessinfo)
            pcu.sessinfo = sessinfo;
        end

        function nCond = get.nCond(pcu)
            if isempty(pcu.condInd)
                nCond = NaN;
            else
                nCond = max(pcu.condInd);
            end
        end

        function csd = compute_csd(pcu)
            % this is the workhorse that calls everything below
            pcu.load_trial_edges();
            pcu.load_trial_alignment();

            tsi = pcu.build_tsi(start="Target", includeOptimalLag=true, windowMs=[0 5000]); % this window doesn't matter
            imec = pcu.load_imec_dataset();

            csd = NHPPixel.CurrentSourceDensityAnalysis(tsi, imec);
            csd.compute_all();
        end
        
        function load_trial_edges(pcu, args)
            arguments
                pcu
                args.regenerate (1, 1) logical = false;
            end

            if ~isempty(pcu.imecTrialEdges_unfiltered) && ~args.regenerate
                return;
            end
            % ultimately loads imecTrialEdges_unfiltered

            paths = pcu.sessinfo.paths;

            % 1) load behavior data            
            T = loadsession(fullfile(paths.sgDataPath, paths.prefixBehavior));
            debug('Session loaded: %d trials \n',size(T,1))

            % 2) load kilosort output
            
            % [spikeIdxMat, clusterID, clusterLabels] = ksResults2spikeMat(paths.ksResultsPath);
            % fprintf('Spike Times loaded %.1f \n', toc)
            % spikeIdxMat_orig = spikeIdxMat;
            % 
            % timing.spikeMatCreated = toc

            % 3) convert spiketimes from imec headstage sample indices into NIDAQ indices:
            % not sure if equivalent?
            % readSpikeGLXmeta = @(file) Neuropixel.readINI(file);
            niMeta = readSpikeGLXmeta(paths.nidaqMetaPath);
            pcu.fsNI = niMeta.sRateHz;
           
            apMeta = readSpikeGLXmeta(paths.npixApMetaPath);
            if isempty(apMeta)
                debug('Warning: using default fsAP\n')
                pcu.fsAP = 30000;
            else
                pcu.fsAP = apMeta.sRateHz;
            end

            pcu.fsSG = 1000; % hardcoded speedgoat (behavior) at 1 kHz
            
            % spikeIdxMatNi = convertSpikeTimeIndices(spikeIdxMat, FsImec, FsNi);
            
            % 4) sync GLX times with behavior, using AP band timescale
            [pcu.Tsync, ~, nidaqTrialEdges] = syncSpeedgoatNidaq(paths.nidaqPath, T, ...
                'SGsyncChan', pcu.sessinfo.sync_chan, 'SGsyncBit',pcu.sessinfo.sync_bit);
            
            % nidaqTrialEdges is nTrials x 2 indices into the nidaq-defined trials, but using AP band sampling
            pcu.imecTrialEdges_unfiltered = uint64(round(nidaqTrialEdges*pcu.fsAP/pcu.fsNI));
        end

        function load_trial_alignment(pcu, args)
            % this does the work that paccond_gain_switch would do, and figures out the alignment windows for each trial
            % ultimately 

            arguments
                pcu
                args.regenerate (1, 1) logical = false;
                args.build_task (1, 1) logical = false;
            end

            if ~isempty(pcu.condInd) && ~args.regenerate
                return;
            end

            T = pcu.Tsync;
            FsSg = pcu.fsSG;
            FsNeural = pcu.fsAP;
            alignState = 'InTarget';
            padDur = [0.5 0.5]; % changed per Eric's email
            errThr = 2; % changed per Eric's email
            stdThr = 3;

            %% Parse trial table
            
            % @djoshea we'll track which trials get thrown away
            trialIds = (1:height(T))'; %#ok<*PROP> 
            
            % write task states to struct
            stateNo = T.Properties.UserData.TaskStates(:,1);
            stateName = T.Properties.UserData.TaskStates(:,2);
            TaskStates = cell2struct(stateNo, stateName);
            
            % verify alignment state exists
            assert(ismember(alignState, stateName))
            
            % filter by save tags
            if ~isempty(pcu.sessinfo.save_tags)
                trialIds = trialIds(ismember(T.saveTag, pcu.sessinfo.save_tags));
                T = T(ismember(T.saveTag, pcu.sessinfo.save_tags), :);
            end

            % remove invalid trials
            valid = T.validTrial;
            T = T(valid, :);
            trialIds = trialIds(valid);
            
            % remove glitches
            glitched = cellfun(@(ts) ismember(TaskStates.Glitch,ts), T.taskState);
            %nGlitch = nnz(glitched);
            trialIds = trialIds(~glitched);
            T = T(~glitched,:);

            %% Initialize task object
            
            % condition fields
            trialFields = fieldnames(T.trialParams{1});
            
            % remove extraneous trial fields
            trialFields = trialFields(~ismember(trialFields,{'stimDelay'}));
            
            nTrialFields = length(trialFields);
            condFields = [{'id'}; trialFields];
            
            % % extract relevant trial parameters
            % if ismember('stimElectrode',trialFields)
            %     trialParams = cellfun(@(tp) [tp.condNo,tp.stimElectrode,tp.stimCurrent], T.trialParams,'uni',false);
            % else
            %     trialParams = cellfun(@(tp) tp.condNo, T.trialParams,'uni',false);
            % end
            
            % Strip some bad trials
            trialMask = ~isnan(T.saveTag);
            T = T(trialMask,:);
            trialIds = trialIds(trialMask);
            
            % extract trial parameters from all trials.
            % if perturbations exist in the dataset, add those to the conditionparameter list
            if nnz(strcmp(fieldnames(T.trialParams{1}),'pertFlag')) > 0
                trialParams = cellfun(@(tp) [tp.frcPol, tp.type, tp.offset, tp.amplitude, tp.duration, tp.frequency, tp.gain, tp.pertFlag, tp.pertAmp, tp.pertTime], T.trialParams,'uni',false);
            else
                trialParams = cellfun(@(tp) [tp.gain, tp.offset, tp.amplitude, tp.duration, tp.frequency], T.trialParams,'uni',false);
            end
            
            % trialParams = cellfun(@(tp) [tp.gain], T.trialParams,'uni',false);
            
            % determine number of unique conditions
            unqCond = unique(cell2mat(trialParams),'rows');
            nCond = size(unqCond,1);
            unqCondNo = (1:nCond)';
            
            % map conditions to unique identifier
            condID = cellfun(@(tp) unqCondNo(ismember(unqCond,tp,'rows')), trialParams);
            
            if length(unique(condID)) < nCond
                warning('All conditions are not present in the data - nCond: %d, present: %d',nCond, length(unique(condID)))
                nCond = length(unique(condID));
%                 unqCond = unqCondNo(ismember(unqCondNo,condID));
            end
            
            % initialize condition table
            Condition = cell2table(cell(nCond,1+nTrialFields), 'VariableNames', condFields);
            
            % condition ID
            Condition.id = (1:nCond)';
            
            % representative set of trial parameters
            Condition(:,trialFields) = cell(nCond,nTrialFields);
            nTrialsMatched = 0;
            for ii = 1:nCond
                trPar = T.trialParams{find(condID == unqCondNo(ii), 1)};
                nTrialsMatched = nTrialsMatched + nnz(condID == unqCondNo(ii));
                for jTF = 1:nTrialFields
                    Condition.(trialFields{jTF}){ii} = trPar.(trialFields{jTF});
                end
                Condition.type{ii} = char(Condition.type{ii});
                Condition.offset{ii} = Condition.offset{ii} * trPar.frcMax;
                Condition.amplitude{ii} = Condition.amplitude{ii} * trPar.frcMax;
            end
            
            assert(nTrialsMatched == numel(trialIds));
            
            % extract scalar entries from cells
            for iTF = 1:nTrialFields
                if all(cellfun(@isscalar,Condition.(trialFields{iTF})))
                    Condition.(trialFields{iTF}) = cell2mat(Condition.(trialFields{iTF}));
                end
            end

            if args.build_task
                Task = PacmanTaskCond(Condition); %(:,trialFields(ismember(trialFields,COND_PARAMS))));
            end
            %% Populate task object
            
            AlignStats = cell(nCond,1);

            % alignIndex is the time lag in ms (speedgoat samples) to the target event
            
            [alignIndices, targetFnOffsetFirstPeak] = deal(cell(nCond,1)); % stores per trial alignment info, for each condition, to be combined at the end
            % optimalLags is the additional lag in ms to align force with the target function
%             optimalLags = cell(nCond, 1);
%             maxNRMSE = cell(nCond,1);
            
            % loop over conditions, align each group of trials internally
            pbar = ProgressBar(nCond, 'Processing Condition data', nCond);
            trialInds_indsByCond = cell(nCond, 1);
            
            for ii = 1:nCond
            %     pbar.update(ii, sprintf('Condition %d of %d',ii,nCond) )
                pbar.update(ii)

                if strcmp(Condition.type{ii},'STA')
                    % skip static condition
                    continue;
                end

                
                % condition indices
                condIdx = find(condID==unqCondNo(ii));
                
                % alignment point per trial
                alignIdx = cellfun(@(ts) find(ts == TaskStates.(alignState),1), T{condIdx,'taskState'}, 'uni', false);
                
                % remove trials without alignment point
                emptyAlignment = cellfun(@isempty,alignIdx);
                condIdx(emptyAlignment) = [];
                alignIdx(emptyAlignment) = [];
                
                % trial indices (Speedgoat and Blackrock)
                trIdxSG = -round(FsSg*padDur(1)) : round(FsSg*(padDur(2)+Condition.duration(ii)));
                trIdxNeural = -round(FsNeural*padDur(1)) : round(FsNeural*(padDur(2)+Condition.duration(ii)));
                
                % check bounds
                isBounded = cellfun(@(ai,fr) diff([ai+trIdxSG(end),length(fr)],[],2)>=0, alignIdx, T.forceRaw(condIdx));
                isBounded = isBounded & cellfun(@(ai,ff) diff([ai+trIdxSG(end),length(ff)],[],2)>=0, alignIdx, T.forceFilt(condIdx));
                if ismember('emgSG',fieldnames(T))
                    isBounded = isBounded & cellfun(@(ai,esg) diff([ai+trIdxSG(end),size(esg,2)],[],2)>=0, alignIdx, T.emgSG(condIdx));
                end
                if ismember('emgBR',fieldnames(T))
                    isBounded = isBounded & cellfun(@(ai,ebr) diff([ai+trIdxNeural(end),size(ebr,1)],[],2)>=0, alignIdx, T.emgBR(condIdx));
                end
                condIdx = condIdx(isBounded);
                alignIdx = alignIdx(isBounded);

                if ~any(isBounded)
                    debug('Warning: Skipping condition %d, none isBounded\n', ii);
                    continue;
                end
            
                % alignment stats
                AlignStats{ii} = struct('emptyPts',nnz(emptyAlignment), 'unbounded',nnz(~isBounded));
                
                alignIndices{ii} = cat(1, alignIdx{:});

                
                colNames = Condition.Properties.VariableNames;
	            colInds = ismember(colNames,{'type','gain','blockMembership','offset','amplitude','duration','frequency','decay','power'});
                
                % phase correction for dynamic conditions. 
                % NOTE: this is not appropriate for any analyses requiring detailed
                % timing (i.e. perturbation response, divergence between sets of
                % conditions, etc.)
            %         warning('Aligning trials to behavior - impacts timing of neural responses')

               % skipping alignment, since it doesn't matter
%                 if ~strcmp(Condition.type{ii},'STA')
                    MAX_LAG = 0.2;
                    padDurTrunc = padDur-2*MAX_LAG;
%         
                    [targFn, ~, ~, ~, tFirstPeakThisCond] = NHPPixel.pacmantargfns(Condition(ii,colInds),1,'dt',1/FsSg,'padDur',padDurTrunc);
%                     tIdx = -round(FsSg*padDurTrunc(1)):round(FsSg*(padDurTrunc(2)+Condition.duration(ii)));
%                     maxLagSamp = round(MAX_LAG*FsSg);
%                     lags = -maxLagSamp:maxLagSamp;
%                     optLag = zeros(length(condIdx),1);
%                     maxNRMSE{ii} = zeros(length(condIdx),1);
%                     for trial = 1:length(condIdx)
%                         normRMSE = zeros(length(lags),1);
%                         for ll = 1:length(lags)
%         
%                             frcAlg = T.forceFilt{condIdx(trial)}(tIdx+alignIdx{trial}+lags(ll))';   
%                             if Condition.gain(ii) < 0
%                                 frcAlg = Condition.frcMax(ii) + Condition.gain(ii)*frcAlg;
%                             end
%                             normRMSE(ll) = 1 - sqrt(mean((frcAlg-targFn).^2)/var(targFn));
%                         end
%                         [~,maxIdx] = max(normRMSE);
%                         optLag(trial) = lags(maxIdx);
%                         maxNRMSE{ii}(trial) = normRMSE(maxIdx);
%                     end
%         
%                     % don't shift alignment index, just compute the optimal lag in ms
%                     optimalLags{ii} = optLag;
% %                   
%                     % and keep track of the time the first peak occurs relative to target in ms
                    targetFnOffsetFirstPeak{ii} = repmat(1000*tFirstPeakThisCond, length(condIdx), 1);
%                       alignIdx = num2cell(cell2mat(alignIdx)+optLag);
%                 end

                % check improper length trials
                if ismember('emgBR',fieldnames(T))
                    badTrials = cellfun(@(x) x*(FsBr/FsSg)+trIdxNeural(end),alignIdx) > cellfun(@(x) size(x,1),T.emgBR(condIdx));
                    condIdx(badTrials) = [];
                    alignIdx(badTrials) = [];
                end
                
                % fix trials too short
                force_lengths = cellfun(@(fr) length(fr),T.forceRaw(condIdx));
                alignedMax = cellfun(@(ai) ai+trIdxSG(end),alignIdx);
                gap = alignedMax - force_lengths;
                if any(gap > 0)
                    Condition.duration(ii) = Condition.duration(ii) - max(gap) / FsSg;
                    trIdxSG = -round(FsSg*padDur(1)) : round(FsSg*(padDur(2)+Condition.duration(ii)));
                end

                badTrials = cellfun(@(fr) length(fr),T.forceRaw(condIdx))...
                    < cellfun(@(ai) ai+trIdxSG(end),alignIdx);
                condIdx(badTrials) = [];
                alignIdx(badTrials) = [];
                

                % remove inaccurate trials
                targFn = pacmantargfns(Condition(ii,colInds),1,'dt',1/FsSg,'padDur',padDur);
                tIdx = -round(FsSg*padDur(1)) : round(FsSg*(padDur(2)+Condition.duration(ii)));
                
                % #TODO: fix this 
            %     y = Condition.frcMax(ii) * cell2mat(cellfun(@(ci,ai,lag) T.forceFilt{ci}(tIdx+ai), num2cell(condIdx), alignIdx,'uni',false))';
                
                y = cell2mat(cellfun(@(ci,ai,lag) T.forceFilt{ci}(tIdx+ai), num2cell(condIdx), alignIdx,'uni',false))';
                if Condition.gain(ii) < 0
                    y = Condition.frcMax(ii) + Condition.gain(ii)*y;
                end
             
                % 1) filter out trials based on maximum deviation from target
                absErr = max(abs(y-targFn),[],1)/max(2,max(targFn)-min(targFn));
                badTrials = absErr > errThr;
                
                condIdx(badTrials) = [];
                alignIdx(badTrials) = [];
                y(:,badTrials) = [];
                
                % 2) remove imprecise trials
                mu = mean(y,2);
                sd = std(y,[],2);
                badTrials = any(y<(mu-stdThr*sd)...
                    | y>(mu+stdThr*sd),1);
                
                condIdx(badTrials) = [];
                alignIdx(badTrials) = [];

                % @djoshea last step - align the entire condition using a threshold 
            
                trialInds_indsByCond{ii} = trialIds(condIdx);
            %     Condition.previousGain{ii} =  previousGain(condIdx);
                
                % forces
            %     forceRaw = cellfun(@(fr,ai,tp) tp.frcMax*(((MAX_FORCE_POUNDS*NEWTONS_PER_POUND)/tp.frcMax * (fr(ai+trIdxSG)/MAX_FORCE_VOLTS)) - tp.frcOff)',...
            %         T.forceRaw(condIdx), alignIdx, T.trialParams(condIdx), 'uni', false);
            %     forceFilt = cellfun(@(ff,ai,tp) Condition.frcMax(ii) * ff(ai+trIdxSG)',...
            %         T.forceFilt(condIdx), alignIdx, T.trialParams(condIdx), 'uni', false);
            
                if args.build_task
                    X_FORCE_IND = 1;
                    Y_FORCE_IND = 2;
                    Z_FORCE_IND = 3;
                    
                    % leaving these here for backwards compatibility prior to three-axis measurements
                    forceRaw = cellfun(@(fr,ai)  fr(Y_FORCE_IND, ai+trIdxSG)',...
                        T.threeAxisForceRaw(condIdx), alignIdx, 'uni', false);
                    forceFilt = cellfun(@(ff,ai) ff(Y_FORCE_IND, ai+trIdxSG)',...
                        T.threeAxisForceFilt(condIdx), alignIdx, 'uni', false);
                
%                     ForceRawX = cellfun(@(fr,ai)  fr(X_FORCE_IND, ai+trIdxSG)', T.threeAxisForceRaw(condIdx), alignIdx, 'uni', false);
                    ForceFiltX = cellfun(@(ff,ai) ff(X_FORCE_IND, ai+trIdxSG)', T.threeAxisForceFilt(condIdx), alignIdx, 'uni', false);    
                    
%                     ForceRawY = cellfun(@(fr,ai)  fr(Y_FORCE_IND, ai+trIdxSG)', T.threeAxisForceRaw(condIdx), alignIdx, 'uni', false);
                    ForceFiltY = cellfun(@(ff,ai) ff(Y_FORCE_IND, ai+trIdxSG)', T.threeAxisForceFilt(condIdx), alignIdx, 'uni', false);    
                
%                     ForceRawZ = cellfun(@(fr,ai)  fr(Z_FORCE_IND, ai+trIdxSG)', T.threeAxisForceRaw(condIdx), alignIdx, 'uni', false);
                    ForceFiltZ = cellfun(@(ff,ai) ff(Z_FORCE_IND, ai+trIdxSG)', T.threeAxisForceFilt(condIdx), alignIdx, 'uni', false);    
                
%                     forces = [{forceFilt}, {forceRaw}];
%                     forces = cellfun(@(x) permute(cell2mat(x'),[1 3 2]), forces, 'uni', false);
                    
                    Fxyz = [{ForceFiltX},{ForceFiltY},{ForceFiltZ}];
                    Fxyz = cellfun(@(x) permute(cell2mat(x'),[1 3 2]), Fxyz,'uni', false);
                    
                    Task.Force(ii).data = cell2mat(Fxyz);  % changed from 
                    Task.Force(ii).Fs = FsSg;
                    Task.Force(ii).alignIndex = 1 + padDur(1)*FsSg;
                    Task.Force(ii).variableLabels = {'Fx','Fy','Fz'};
                end
              
            end
            pbar.finish('Done');
            
            [trialInds, condInd] = TensorUtils.catWhich(1, trialInds_indsByCond{:}); %#ok<*PROPLC> 
            [trialInds, sortIdx] = sort(trialInds);
%             trialIds = trialIds(sortIdx);
            condInd = condInd(sortIdx);

            % each of these in ms

            % time relative from start of trial to InTarget event
            alignIndices = cat(1, alignIndices{:});
            alignIndices = alignIndices(sortIdx);
            
            % additional lag to add to align force with target force
%             optimalLags = cat(1, optimalLags{:});
%             optimalLags = optimalLags(sortIdx);

            % time from InTarget event to first peak of force for each condition
            targetFnOffsetFirstPeak = cat(1, targetFnOffsetFirstPeak{:});
            targetFnOffsetFirstPeak = targetFnOffsetFirstPeak(sortIdx);

            % filter down the trial edges to trials that were kept
            trialStart = pcu.imecTrialEdges_unfiltered(trialInds, 1);

            msToAP = @(x) uint64(x * pcu.fsAP / pcu.fsSG);
            
            pcu.imecTrialEdges = pcu.imecTrialEdges_unfiltered(trialInds, :);
            pcu.apIdx_Target = trialStart + msToAP(alignIndices);
%             pcu.apIdx_TargetWithOptimalLag = pcu.apIdx_Target + msToAP(optimalLags);
            pcu.apIdx_FirstPeak = trialStart + msToAP(alignIndices + targetFnOffsetFirstPeak);
%             pcu.apIdx_FirstPeakWithOptimalLag = pcu.apIdx_FirstPeak + msToAP(optimalLags);

            % swap in the modified version of the conditions table with previousGain
            % field added
            if args.build_task
                Task.Conditions = Condition;
                pcu.task = Task;
            end

            pcu.trialInds = trialInds;
            pcu.condInd = condInd;
            
        end

        function imec = load_imec_dataset(pcu, args)
            arguments
                pcu
                args.lf_rmsRange = pcu.sessinfo.lf_rmsRange;
                args.show_lf_rms_histogram (1, 1) logical = false;
            end

            channelMap = Neuropixel.ChannelMap("~/npl/nhp-pixel/neuropixNHP_kilosortChanMap_v1.mat");
            imec = Neuropixel.ImecDataset(pcu.sessinfo.paths.npixLfpPath, channelMap=channelMap);

            if ~isempty(args.lf_rmsRange)
                imec.markBadChannelsByRMS(band="lf", rmsRange=args.lf_rmsRange, printMessage=true);
            end

            if args.show_lf_rms_histogram
                rmsByChannel = imec.computeRMSByChannel(band="lf");
                histogram(rmsByChannel)
                xlabel('RMS (μV)');
                xline(args.lf_rmsRange(1));
                xline(args.lf_rmsRange(2));
                niceGrid
            end

            pcu.imec = imec;
        end

        function tsi = build_tsi(pcu, args)
            arguments
                pcu
                args.start (1, 1) string
                args.includeOptimalLag (1, 1) logical = false;
                
                args.windowMs = []; % relative to condition alignment point (first peak)
                args.condInds (:, 1) = 1:pcu.nCond;
            end

            trialIndsList = (1:numel(pcu.trialInds))';

            switch args.start
                case "TrialStart"
                    imecZero = pcu.imecTrialEdges(:, 1);
                case "Target"
                    if args.includeOptimalLag
                        imecZero = pcu.apIdx_Target;
                    else
                        imecZero = pcu.apIdx_TargetWithOptimalLag;
                    end
                case "FirstPeak"
                    if args.includeOptimalLag
                        imecZero = pcu.apIdx_FirstPeak;
                    else
                        imecZero = pcu.apIdx_FirstPeakWithOptimalLag;
                    end
                otherwise
                    error("Unknown")
            end

            window_imecAP = round(args.windowMs * pcu.fsAP / 1000);
            if isempty(args.windowMs)
                % go until end of trial
                idxStart = imecZero;
                idxStop = pcu.imecTrialEdges(:, 2);
            else
                % take window around start
                idxStart = uint64(imecZero + window_imecAP(1));
                idxStop = uint64(imecZero + window_imecAP(2));
            end

            % mask by condition
            trialIndsList = trialIndsList(ismember(pcu.condInd, args.condInds));
            nTrials = numel(trialIndsList);
            
            tsi = Neuropixel.TrialSegmentationInfo(nTrials, pcu.fsAP);
            for iiT = 1:nTrials
                iT = trialIndsList(iiT);
                tsi.trialId(iiT) = pcu.trialInds(iiT);
                tsi.conditionId(iiT) = pcu.condInd(iT);
                tsi.idxStart(iiT) = idxStart(iT);
                tsi.idxStop(iiT) = idxStop(iT);
            end 
        end
    end


end
