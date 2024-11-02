classdef DatasetManager < handle
    properties
        sites (:, 1) NHPPixel.PacmanSessionInfo % list of pre-loaded sites to be accessed with get_site_by_id

        common (1, 1) struct
        
        data_root (1, 1) string
        computed_root (1, 1) string % OPTOREACH_COMPUTED_ROOT
    end
    
    properties(Dependent)
        nSites
        site_ids
    end
    
    methods
        function mgr = DatasetManager()
            nhpPixelSettings;            

            mgr.computed_root = get_path_from_env('computed_root', 'NHP_PIXEL_COMPUTED_ROOT');
            mgr.data_root = get_path_from_env('data_root', 'NHP_PIXEL_DATA_ROOT');            
            
            mgr.sites = NHPPixel.buildSites(mgr);
            
            function path = get_path_from_env(name, key)
                path = getenv(key);
                if isempty(path)
                    error('Missing environment variable ''%s''', name, key);
                end
                assert(exist(path, 'dir') > 0, '%s path %s does not exist', name, path);
            end
        end


        % switch statement that maps key values to the load methods below
        function [val, should_save] = generate_object_for_key(mgr, site, key, varargin)
            key = string(key);
            should_save = true; % set to false below if it's quick to construct on the fly
            switch key
                case 'imec'
                    val = mgr.build_imec_dataset(site);
                    should_save = false;

                case 'ks'
                    val = mgr.build_kilosort_dataset(site);
                    should_save = false;

                case 'csd_util'
                    val = mgr.build_csd_util(site);
                    should_save = false;

                case 'current_source_density'
                    val = mgr.build_current_source_density(site);
                    
                otherwise
                    val = [];
                    should_save = false;
            end
        end
        
        function val = modify_post_load_for_key(mgr, site, key, val, varargin)
            key = string(key);
            
        end

         % switch statement that maps key values to the load methods below
        function [val, should_save] = generate_object_for_key_common(mgr, key, varargin)
            key = string(key);
            should_save = true; % set to false below if it's quick to construct on the fly
            switch key
                
                otherwise
                    val = [];
                    should_save = false;
            end
        end
        
        function val = modify_post_load_for_key_common(mgr, key, val, varargin)
            key = string(key);
            switch key
                    
            end
        end
        
    end
    
    methods % Specific analyses
        function imec = build_imec_dataset(mgr, site, varargin) %#ok<*INUSL>
            channelMap = Neuropixel.ChannelMap("~/npl/nhp-pixel/neuropixNHP_kilosortChanMap_v1.mat");
            imec = Neuropixel.ImecDataset(site.paths.npixLfpPath, channelMap=channelMap);
        end

        function ks = build_kilosort_dataset(mgr, site, varargin)
            p = inputParser();
            p.addParameter('deduplicate_spikes', true, @islogical);
            p.addParameter('deduplicate_cutoff_spikes', true, @islogical);
            p.addParameter('deduplicate_within_samples', 5, @isscalar);
            p.addParameter('deduplicate_within_distance', 50, @isscalar);
            p.KeepUnmatched = true; % pass load? flags to ks.load(...)
            p.parse(varargin{:});

            channelMap = Neuropixel.ChannelMap("~/npl/nhp-pixel/neuropixNHP_kilosortChanMap_v1.mat");
            ks_path = site.paths.ksResultsPath;

            imec = mgr.build_imec_dataset(site);
            ks_args.imecDataset = imec;
                
            ks = Neuropixel.KilosortDataset(ks_path, 'channelMap', channelMap, ks_args);
            
            if ~ks.is_deduplicated && (p.Results.deduplicate_spikes || p.Results.deduplicate_cutoff_spikes)
                % force a load so that the initial spike deduplication mask is loaded
                ks.load('loadFeatures', false, p.Unmatched);
            end
        end

        function pcu = build_csd_util(mgr, site, varargin)
            pcu = NHPPixel.PacmanCSDUtil(site);
        end
        
        function csd = build_current_source_density(mgr, site, varargin)
            pcu = NHPPixel.PacmanCSDUtil(site);
            csd = pcu.compute_csd();
        end
    end

    methods % Specific common analyses
        
    end
    
    methods % path building
        function path = build_computed_path_for_site(mgr, site)
            path = fullfile(mgr.computed_root, site.id);
        end
        
        function file = get_save_location_for_key(mgr, site, key)
            file = fullfile(mgr.build_computed_path_for_site(site), sprintf('%s.mat', key));
        end

        function path = build_computed_path_common(mgr)
            path = fullfile(mgr.computed_root, "common");
        end

        function file = get_save_location_for_key_common(mgr, key)
            file = fullfile(mgr.build_computed_path_common(), sprintf('%s.mat', key));
        end
    end
    
    methods % Per-site key based loading / caching infrastructure
        function n = get.nSites(mgr)
            n = numel(mgr.sites);
        end
        
        function ids = get.site_ids(mgr)
            ids = cat(1, mgr.sites.id);
        end
     
        % use this to access the canonical StimSite object for each site
        function site = get_site_by_id(mgr, id)
            id = string(id);
            ids = cat(1, mgr.sites.id);
            [tf, ind] = ismember(id, ids);
            assert(all(tf));
            
            site = mgr.sites(ind);
        end
        
        function out = load_all_objects_for_key(mgr, key)
            key = string(key);
            out = cell(numel(mgr.sites), 1);
            for iS = 1:numel(mgr.sites)
                site = mgr.sites(iS);
                out{iS} = mgr.load_key(site, key);
            end
            out = cat(1, out{:});
        end

        function val = load_key(mgr, site, key, varargin)
            val = mgr.get_or_load_keys(site, key, struct(), varargin{:});
        end
        
        function varargout = get_or_load_keys(mgr, site, keys, override, varargin)
            p = inputParser();
            p.addParameter('regenerate', false, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            regenerate = p.Results.regenerate;
            
            % either retrieve a key from site.data or load it using 
            if nargin < 4
                override = struct();
            end
            if isempty(override)
                override = struct();
            end
            keys = string(keys);
            varargout = cell(numel(keys), 1);
            for iK = 1:numel(keys)
                key = keys(iK);
                
                if contains(key, "|")
                    % accept any of the specified keys
                    options = strsplit(key, "|");
                    for iO = 1:numel(options)
                        val = search_single_key(options(iO));
                        if ~isempty(val)
                            key = options(iO); % make sure we store the found object in the right key
                            break;
                        end
                    end
                    if isempty(val)
                        % none found, load the first option as key
                        key = options(1);
                    end
                else
                    val = search_single_key(key);
                end
                
                if isempty(val)
                    % not loaded, generate using do_load_for_key lookup
                    val = mgr.load_or_generate_key(site, key, 'regenerate', p.Results.regenerate, p.Unmatched);
                    assert(~isempty(val), 'Key %s loaded as empty', key);
                    
                    % and store the object in site
                    site.set_data(key, val);
                end
                varargout{iK} = val;
            end
            
            function val = search_single_key(key)
                if regenerate
                    val = [];
                elseif isfield(override, key) && ~isempty(override.(key))
                    val = override.(key);
                    
                elseif site.has_data(key) && ~isempty(site.get_data(key))
                    val = site.get_data(key);
                else
                    val = [];
                end
            end
        end
        
        function varargout = global_get_or_load_keys(mgr, sites, keys, overrides, varargin)
            nS = numel(sites);
            nK = numel(keys);
            if nargin < 4
                overrides = repelem(struct(), nS, 1);
            end
            if isempty(overrides)
                overrides = repelem(struct(), nS, 1);
            end
            
            varargout = cell(nK, 1);
            for iK = 1:numel(keys)
                for iS = 1:numel(sites)
                    varargout{iK}(iS, 1) = mgr.get_or_load_keys(sites(iS), keys(iK), overrides(iS), varargin{:});
                end
            end
        end
        
        function val = load_or_generate_key(mgr, site, key, varargin)
            % this will either load key from disk (using get_save_location_for_key) or if not found, 
            % compute it (using generate_object_for_key) and save to disk
            
            p = inputParser();
            p.addParameter('regenerate', false, @islogical);
            p.addParameter('save_to_disk', true, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            file = mgr.get_save_location_for_key(site, key);
            if exist(file, 'file') && ~p.Results.regenerate
                % load from file in imec cleaned directory
                fprintf('Loading site %s key %s from %s\n', site.id, key, file);
                ld = load(file, key);
                val = ld.(key);
                
                % enable post load hooks that compute anything that needs to be computed on load
                val = mgr.modify_post_load_for_key(site, key, val);
                
            else
                fprintf('Generating site %s key %s\n', site.id, key);
                args = namedargs2cell(p.Unmatched);
                [val, should_save] = mgr.generate_object_for_key(site, key, args{:});
                if p.Results.save_to_disk && should_save
                    fprintf('Saving site %s key %s to %s\n', site.id, key, file);
                    mkdirRecursive(fileparts(file));
                    
                    to_save.(key) = val; %#ok<STRNU>
                    saveLarge(file, '-struct', 'to_save');
                end
            end
        end
        
        function save_modified_key_to_disk(mgr, site, key, value)
            file = mgr.get_save_location_for_key(site, key); 
            fprintf('Saving site %s key %s to %s\n', site.id, key, file);
            to_save.(key) = value;
            save(file, '-struct', 'to_save');
        end
        
        function clear_data(mgr)
            for iS = 1:mgr.nSites
                mgr.sites(iS).clear_data();
            end
        end
    end

    methods % Common data key loading infrastructure
        function val = load_key_common(mgr, key, varargin)
            val = mgr.get_or_load_keys_common(key, struct(), varargin{:});
        end

        function varargout = get_or_load_keys_common(mgr, keys, override, varargin)
            p = inputParser();
            p.addParameter('regenerate', false, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            regenerate = p.Results.regenerate;
            
            % either retrieve a key from site.data or load it using 
            if nargin < 4
                override = struct();
            end
            if isempty(override)
                override = struct();
            end
            keys = string(keys);
            varargout = cell(numel(keys), 1);
            for iK = 1:numel(keys)
                key = keys(iK);
                
                if contains(key, "|")
                    % accept any of the specified keys
                    options = strsplit(key, "|");
                    for iO = 1:numel(options)
                        val = search_single_key(options(iO));
                        if ~isempty(val)
                            key = options(iO); % make sure we store the found object in the right key
                            break;
                        end
                    end
                    if isempty(val)
                        % none found, load the first option as key
                        key = options(1);
                    end
                else
                    val = search_single_key(key);
                end
                
                if isempty(val)
                    % not loaded, generate using do_load_for_key lookup
                    val = mgr.load_or_generate_key_common(key, 'regenerate', p.Results.regenerate, p.Unmatched);
                    assert(~isempty(val), 'Key %s loaded as empty', key);
                    
                    % and store the object in site
                    mgr.common.(key) = val;
                end
                varargout{iK} = val;
            end
            
            function val = search_single_key(key)
                if regenerate
                    val = [];
                elseif isfield(override, key) && ~isempty(override.(key))
                    val = override.(key);
                    
                elseif isfield(mgr.common, key) && ~isempty(mgr.common.(key))
                    val = mgr.common.(key);
                else
                    val = [];
                end
            end
        end

        function val = load_or_generate_key_common(mgr, key, varargin)
            % this will either load key from disk (using get_save_location_for_key) or if not found, 
            % compute it (using generate_object_for_key) and save to disk
            
            p = inputParser();
            p.addParameter('regenerate', false, @islogical);
            p.addParameter('save_to_disk', true, @islogical);
            p.KeepUnmatched = true;
            p.parse(varargin{:});
            
            file = mgr.get_save_location_for_key_common(key);
            if exist(file, 'file') && ~p.Results.regenerate
                % load from file in imec cleaned directory
                fprintf('Loading common key %s from %s\n', key, file);
                ld = load(file, key);
                val = ld.(key);
                
                % enable post load hooks that compute anything that needs to be computed on load
                val = mgr.modify_post_load_for_key_common(key, val);
                
            else
                fprintf('Generating common key %s\n', key);
                args = namedargs2cell(p.Unmatched);
                [val, should_save] = mgr.generate_object_for_key_common(key, args{:});
                if p.Results.save_to_disk && should_save
                    fprintf('Saving common key %s to %s\n', key, file);
                    mkdirRecursive(fileparts(file));
                    
                    to_save.(key) = val; %#ok<STRNU>
                    saveLarge(file, '-struct', 'to_save');
                end
            end
        end    
    end
   
end