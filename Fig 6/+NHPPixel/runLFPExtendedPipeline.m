function [info, imec_cat] = runLFPExtendedPipeline(imec, tsi, varargin)
% likely not needed!

return

%%
isstringlike = @(x) ischar(x) || isstring(x);
p = inputParser();
p.addParameter('cleanedStem', '', isstringlike);
p.addParameter('rmsRange', [30 140], @isvector); % for LFP
p.addParameter('dryRun', false, @islogical);
p.parse(varargin{:});
dataRoot = char(p.Results.dataRoot);
dryRun = p.Results.dryRun;
estim = p.Results.estim;

if exist('writeNPY', 'file') ~= 2
    error('npy-matlab was not found on path');
end

imec.setSyncBitNames(7, "speedgoat");

% if imec.hasAP
%     sqBadChannels = NeuropixelExpt.DataClean.detectSquarewaveChannels(imec);
%     debug('Marking %d channels bad based on detected square waves\n', numel(sqBadChannels));
%     imec.markBadChannels(sqBadChannels);
% end

rmsBadChannels = imec.markBadChannelsByRMS(band="lf", rmsRange=p.Results.rmsRange);
debug('Marking %d channels bad based on rms\n', numel(rmsBadChannels));

% Load the trial segmentation info so we can excise the outside trial regions later
if tsi.fs == imec.fsAP
    % tsi built from AP sync line, convert to LF sample rate
    tsi = tsi.convertToDifferentSampleRate(imec.fsLF);
end

% figure out paths for the cleaned data and the kilosort location
if isempty(p.Results.cleanedStem)
    % use stem from raw file
    cleanedStem = char(imec.fileStem);    
else
    cleanedStem = char(p.Results.cleanedStem);
end
cleanedBinFile = [cleanedStem sprintf('.imec%d.lf.bin', imec.fileImecNumber)];

[parent, leaf] = fileparts(char(imec.pathRoot));
parent = [parent, '_lfpCleaned'];
cleanedPath =  fullfile(parent, leaf, cleanedBinFile);

% print path info so user can check
debug('Path info:\n');
debug('  Raw file: %s\n', imec.pathLF);
debug('  Cleaned : %s\n', cleanedPath);


% flush the unused sync bits to the cleaned_datasets and mark the artifact regions
% also excise the regions outside of the trials
debug('Writing sync-corrected concatenated lfp at %s\n', cleanedPath);

% constrain each file to sample window if requested
extraMeta = struct();

% allow override of syncBits connected, but default to 1-3 for estim and 1-2 if not
syncBitsConnected = p.Results.syncBitsConnected;
if isempty(syncBitsConnected)
    if estim
        syncBitsConnected = 1:3;
    else
        syncBitsConnected = 1:2;
    end
end
clearSyncBitsFn = @(varargin) PrimatePixel.DataClean.clearSpecifiedUnusedSyncBits(syncBitsConnected, varargin{:});
lf_fnList = {clearSyncBitsFn};

runKilosort = p.Results.runKilosort;
if runKilosort
    debug('Writing concatenated AP band as well as LF in order to run Kilosort for LFP spike removal\n');
    writeAP = true;

    extraMeta.commonAverageReferenced = true; % haappens inside detectAndZeroStimArtifactWindows

    ap_extraArg = struct();
    if p.Results.hpFilterPostCAR % needed for a minority of estim datasets where the post-stim transient lasts a while
        extraMeta.hpFilteredPostCAR = true;
        extraMeta.hpFilterCornerHz = p.Results.hpFilterCornerHz;
        extraMeta.hpFilterHalfOrder = p.Results.hpFilterHalfOrder;

        debug('Running CAR on AP band with HP filtering above %g Hz\n', extraMeta.hpFilterCornerHz );
        ap_extraArg.hp_filter = true;
        ap_extraArg.hp_filter_corner = extraMeta.hpFilterCornerHz;
        ap_extraArg.fs = imecList{1}.fsAP;
        ap_extraArg.hp_filter_half_order = extraMeta.hpFilterHalfOrder;
    else
        debug('Running CAR on AP band\n');
    end

    if estim
        % preserve stim bit, zero during stim, CAR inside detectAndZeroStimArtifactWindows
        extraMeta.run_detectAndZeroStimArtifactWindows = true;
        ap_fnList = {clearSyncBitsFn; @PrimatePixel.DataClean.detectAndZeroStimArtifactWindows};
    else
        % also clean stim bit, don't worry about zeroing stim, just CAR
        extraMeta.run_detectAndZeroStimArtifactWindows = false;
        ap_fnList = {clearSyncBitsFn; @Neuropixel.DataProcessFn.commonAverageReferenceEachBank};
    end
    chunkEdgeExtraSamplesAP = [10000 10000];
else
    writeAP = false;
    ap_fnList = {};
    ap_extraArg = struct();
    chunkEdgeExtraSamplesAP = [0 0];
end

% this function will automatically convert timeShiftsLF to the corresponding timeShiftsAP if writeAP == true
imec_cat = Neuropixel.ImecDataset.writeConcatenatedFileMatchGains(imecList, cleanedPath, ...
    'writeAP', writeAP, 'transformAP', ap_fnList, 'chunkEdgeExtraSamplesAP', chunkEdgeExtraSamplesAP, 'transformAPExtraArg', ap_extraArg, ...
    'writeLF', true, 'transformLF', lf_fnList, 'extraMeta', extraMeta, 'timeShiftsLF', timeShiftsLF, ...
    'dryRun', dryRun);

info.rawPaths = imecRawPaths;
info.cleanedPath = cleanedPath;

if p.Results.runKilosort
    info.ksPath = ksPath;
    if dryRun
        warning('Skipping KS step due for dryRun true');
        return;
    end

    if ~exist(ksPath, 'dir')
        mkdirRecursive(ksPath);
    end

    % sym link into ks directory
    debug('Sym-linking to %s\n', ksPath);
    imec_ks = imec_cat.symLinkAPIntoDirectory(ksPath);

    % and run KiloSort
    debug('Running PrimatePixel vanilla KiloSort2 on (extended LFP) concatenated AP data\n');
    PrimatePixel.Kilosort2.runKilosort2(imec_ks, 'export_batchwise', false, p.Results.kilosortOpts);
end

end
