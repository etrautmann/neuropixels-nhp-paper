classdef PacmanSessionInfo
    properties(Transient)
        mgr (:, 1) NHPPixel.DatasetManager
        data (1, 1) struct
    end

    properties
        id (1, 1) string
        data_root (1, 1) string
        subject (1, 1) string
        date (1, 1) datetime
        npix_region (1, 1) string
        save_tags (1, :) {mustBeInteger} = []
        gNum (1, 1) {mustBeInteger} = 0
        tNum (1, 1) {mustBeInteger} = 0
        imecNum (1, 1){mustBeInteger} = 0
        sync_chan (1, 1) {mustBeInteger} = 3
        sync_bit (1, 1) {mustBeInteger} = 8;

        lf_rmsRange = [30 150];
    end

    properties(SetAccess=protected)
        paths (1, 1) struct
    end

    methods
        function sess = PacmanSessionInfo(mgr, args)
            arguments
                mgr (:, 1) NHPPixel.DatasetManager
                args.?NHPPixel.PacmanSessionInfo
            end

            sess.mgr = mgr;
            sess.data_root = mgr.data_root;

            flds = fieldnames(args);
            for iF = 1:numel(flds)
                sess.(flds{iF}) = args.(flds{iF});
            end

            if sess.id == ""
                sess.id = sess.subject + "_" + datestr(sess.date, "yyyymmdd");
            end

            sess.paths = sess.buildPacmanPaths();
        end

        function paths = buildPacmanPaths(sess)
            % generates paths for neuropixels data files assuming the following folder
            % structure:
            %
            % - data_root
            % 	- subject
            % 		- raw
            % 			- date
            % 				- speedgoat
            % 				- neuropixels
            %                   - recording gNum folder
            %       				- probe1 folder
            %                       	- file.ap.bin
            %                       	- file.lf.bin
            %                       	- file.ap.meta
            %                       	- file.lf.meta
            %                   - nidq.bin
            %                  	- nidq.meta
            % 				- blackrock
            % 		- processed
            % 			- date
            %               - kilosort-manually-sorted
            %                   - recording gNum folder
            %                       - probe1 sort folder
            %                       - proben sort folder
            % 
            %
            %
            % EMT 2021-03-08
            
            data_root = char(sess.data_root);
            subject = char(sess.subject); %#ok<*PROP> 
            date = datestr(sess.date, 'yyyy-mm-dd');
            gNum = sess.gNum;
            tNum = sess.tNum;
            imecNum = sess.imecNum;

            % 0) *** prefixes
            subjInit = lower(subject(1));
            
            prefix = ['pacman-task_' subjInit '_' date([3,4,6,7,9,10])];
            paths.prefix = prefix;
            
            prefixBehavior = [prefix '_beh'];
            paths.prefixBehavior = prefixBehavior;
            
            prefixNpix = ['pacman-task_' subjInit '_' date([3,4,6,7,9,10]) '_neu'];
            paths.prefixNpix = prefixNpix;
            
            % 1) *** behavioral data via speedgoat
            sgDataPath = fullfile(data_root, subject,'raw', date, 'speedgoat');
            warnIfNotExist(sgDataPath)
            paths.sgDataPath = sgDataPath;
            
            
            % 2) *** Neuropixels data files (see folder structure above for documentation of expected paths)]
            gNumFolder = [prefixNpix '_g' num2str(gNum)];
            probeFolder = [prefixNpix '_g' num2str(gNum) '_imec' num2str(imecNum)];  % #TODO: add support for multiple probes here
            
            % build paths
            nPixRootPath = fullfile(data_root, subject, 'raw', date, 'neuropixels', gNumFolder); 
            warnIfNotExist(nPixRootPath)
            
            % Probe path
            nPixProbePath = fullfile(nPixRootPath, probeFolder);        % #TODO: add support for multiple probes here
            warnIfNotExist(nPixProbePath)
            
            % Raw data paths
            npixApPath = fullfile(nPixProbePath, [prefixNpix '_g' num2str(gNum) '_t' num2str(tNum) '.imec' num2str(imecNum) '.ap.bin']);
%             warnIfNotExist(npixApPath)
            paths.npixApPath = npixApPath;
            
            npixApMetaPath = fullfile(nPixProbePath, [prefixNpix '_g' num2str(gNum) '_t' num2str(tNum) '.imec' num2str(imecNum) '.ap.meta']);
%             warnIfNotExist(npixApMetaPath)
            paths.npixApMetaPath = npixApMetaPath;
            
            npixLfpPath  =  fullfile(nPixProbePath, [prefixNpix '_g' num2str(gNum) '_t' num2str(tNum) '.imec' num2str(imecNum) '.lf.bin']);
            warnIfNotExist(npixLfpPath)
            paths.npixLfpPath = npixLfpPath;
            
            npixLfpMetaPath  =  fullfile(nPixProbePath, [prefixNpix '_g' num2str(gNum) '_t' num2str(tNum) '.imec' num2str(imecNum) '.lf.meta']);
            warnIfNotExist(npixLfpMetaPath)
            paths.npixLfpMetaPath = npixLfpMetaPath;
            
            % 2.1 NIDAQ I/O
            nidaqPath = fullfile(nPixRootPath, [prefixNpix '_g' num2str(gNum) '_t' num2str(tNum) '.nidq.bin']);
            warnIfNotExist(nidaqPath)
            paths.nidaqPath = nidaqPath;
            
            nidaqMetaPath = fullfile(nPixRootPath, [prefixNpix '_g' num2str(gNum) '_t' num2str(tNum) '.nidq.meta']);
            warnIfNotExist(nidaqMetaPath)
            paths.nidaqMetaPath = nidaqMetaPath;
            
            % Manually-sorted kilosort output 
            ksResultsPath = fullfile(data_root, subject, 'processed', date, 'kilosort-manually-sorted',gNumFolder, probeFolder);
%             warnIfNotExist(ksResultsPath)
            paths.ksResultsPath = ksResultsPath;
            
            % 3) *** blackrock NSX path
            brDataPath = fullfile(data_root, subject, 'raw', date, 'blackrock', ['pacman-task_' subjInit '_' date([3,4,6,7,9,10]) '_emg_001.ns6']);
%             warnIfNotExist(brDataPath)
            paths.brDataPath = brDataPath;

            % task output
            paths.taskTableOutputPath = fullfile(data_root, subject, 'processed', date, 'mergedTaskData', [paths.prefix '_taskdata.mat']);

            function [] = warnIfNotExist(path)
                if exist(path, 'dir') == 0 && exist(path, 'file') == 0
                    debug(['Warning: ' char(path) ' doest not exist\n'])
                end
            end
        end

        function out = load_key(site, key, varargin)
            % uses DatasetManger.load_key to load something
            out = site.mgr.load_key(site, key, varargin{:});
        end
        
        function out = load_or_generate_key(site, key, varargin)
            % uses DatasetManger.load_or_generate_key to load something
            out = site.mgr.load_or_generate_key(site, key, varargin{:});
        end
        
        function path = build_computed_path(site)
            path = site.mgr.build_computed_path_for_site(site);
        end
    end

    methods % Utilities for accessing / setting data & meta keys
        function meta = getSessionMeta(stg, key, default) % for StimDynamics.StimSession
            if nargin == 1
                meta = stg.meta;
            elseif nargin == 2
                meta = stg.get_meta(key);
            else
                meta = stg.get_meta(key, default);
            end 
        end
        
        function v = get_meta(stg, fld, or)
            if isfield(stg.meta, fld)
                v = stg.meta.(fld);
            elseif nargin > 2
                v = or;
            else
                v = [];
            end
        end

        function set_meta(stg, fld, value)
            stg.meta.(fld) = value;
        end

        function tf = has_meta(stg, fld)
            tf = isfield(stg.meta, fld);
        end

        function clear_meta(stg, fld)
            if isfield(stg.meta, fld)
                stg.meta = rmfield(stg.meta, fld);
            end
        end

        function v = get_data(stg, fld, or)
            if isfield(stg.data, fld)
                v = stg.data.(fld);
            elseif nargin > 2
                v = or;
            else
                v = [];
            end
        end

        function set_data(stg, fld, value)
            stg.data.(fld) = value;
        end

        function clear_data(stg, fld)
            if nargin < 2
                stg.data = struct();
            elseif isfield(stg.data, fld)
                stg.data = rmfield(stg.data, fld);
            end
        end

        function tf = has_data(stg, fld)
            tf = isfield(stg.data, fld);
        end
    end
end