function [d] = buildCalciumImagingTable(root_directory, guide_name, ...
    experiment_type, n_movies, n_stims, n_movie_frames, n_cluster_frames, n_pcs, lr_params, not_lr_params)
% Builds a table for calcium imaging data belonging to specific neurons.
% Each row is for one neuron and holds a lot of metadata about the neuron,
% in addition to the calcium imaging traces for n_movies, each
% n_movie_frames long. The neurons should all be from the same
% experiment_type, e.g. 'intensity response' and must have the same number
% of movie frames in each movie. 

    % Set up an empty array to hold the Neuron objects
    neurons = []; 
    
    % Set up an empty array to keep track of which rows of the guide table
    % correspond to each neuron 
    guide_rows = []; 
    
    % Load the data guide spreadsheet
    guide = readtable(guide_name); 
    
    % Manual exclusion of datasets 
    guide = guide(logical(guide.include), :); 
    
    % Select data corresponding to the correct experiment type
    guide = guide(strcmp(guide.experiment_type, experiment_type), :); 

    % For each row of the data guide spreadsheet... 
    for r = 1:length(guide.data_name)
        
        % Start at the root directory 
        cd(root_directory); 
    
        % Load the Retina object 
        data_name = guide.data_name{r}; 
        folder = guide.folder{r};
        directory = guide.directory{r}; 
        cd([directory, '\', folder]);
        R = load([data_name, '.mat']);
        
        % For whatever reason, matlab loads retina objects as structs with
        % a '.R' field. Go in to the struct and get out the retina object. 
        if isstruct(R)
            field_names = fields(R);
            field_name = field_names{1};
            R = R.(field_name);
        end 
        
        % Neuron IDs for this retina piece
        ids = R.getNeuronIDList(); 
        
        % For each Neuron object...
        for i = 1:length(ids) 
        
            % Add the Neuron to the array of Neurons 
            neurons = [neurons, R.getNeuron(ids{i})];
            
            % Keep track of which row of the data guide table it came from. 
            guide_rows = [guide_rows, r];
            
        end
    end 
    
    % Create the table with pre-allocated columns 
    var_names = {'date', 'age', 'genotype', 'mouse_id', 'data_name', 'neuron_id', 'neuron_object', 'motion_artifact',... % Metadata
        'conditions', 'framerate', ... % Experimental parameters
        'gfp_pos_thresh', 'gfp_neg_thresh', 'fluorescent_marker_name', 'fluorescent_marker_intensity', ... % morphological and molecular information
        'stratification', 'immuno_marker', 'subtype', 'subtype_code', 'mouse_code'};
    var_types = {'double', 'uint8', 'string', 'double', 'string', 'double', 'Neuron', 'logical',... % Metadata
        'string', 'double', ...   % Experimental parameters
        'double', 'double', 'string', 'double', ... % Morphological and molecular information 
        'string', 'string', 'string', 'double', 'double'}; 
    sz = [length(neurons) length(var_names)];
    d = table('Size', sz, 'VariableTypes', var_types, 'VariableNames', var_names);
    
    % Pre-allocate the calcium imaging data arrays - a raw and preprocessed 'clustering' version for each possible
    % movie along with the snr and z score columns 
    for m = 1:n_movies
        
        d.(['movie_', num2str(m), '_data']) = NaN * ones(length(neurons), n_movie_frames);
        d.(['movie_', num2str(m), '_clustering_data']) = NaN * ones(length(neurons), n_cluster_frames);
        d.(['movie_', num2str(m), '_feature_vector']) = NaN * ones(length(neurons), n_pcs);
        d.(['movie_', num2str(m), '_max_responses']) = NaN * ones(length(neurons), n_stims);
        d.(['movie_', num2str(m), '_snr']) = NaN * ones(length(neurons), 1); 
        d.(['movie_', num2str(m), '_max_z_score']) = NaN * ones(length(neurons), 1);
        d.(['movie_', num2str(m), '_lr']) = zeros(length(neurons), 1); 
        d.(['movie_', num2str(m), '_not_lr']) = zeros(length(neurons), 1); 
    end
    
    % Pre-allocate the stim frame and intensity arrays
    d.stim_on_frames = NaN * ones(length(neurons), n_stims); 
    d.stim_off_frames = NaN * ones(length(neurons), n_stims);
    d.intensity = NaN * ones(length(neurons), n_stims);
 
    
    % Set up a mapping between full mouse id numbers and an integer
    % numbering of mice
    mouse_id_map = containers.Map(unique(guide.mouse_id), 1:length(unique(guide.mouse_id)));
    
    % For each Neuron, fill in a row of the table 
    for i = 1:length(neurons)
        
        n = neurons(i); 
        r = guide_rows(i); 
        
        % Metadata
        d.date(i) = guide.date(r);
        d.age(i) = guide.age(r);
        d.genotype{i} = guide.genotype{r};
        d.mouse_id(i) = guide.mouse_id(r); 
        d.mouse_code(i) = mouse_id_map(guide.mouse_id(r));  % Also assign an integer numbering to the mice 
        d.data_name{i} = guide.data_name{r}; 
        d.neuron_id(i) = 10^10 * guide.mouse_id(r) + i; % Generate a unique neuron ID
        d.neuron_object{i} = n;
        d.motion_artifact(i) = sum(d.neuron_id(i) == str2double(guide.manual_lr_exclusion{r})) > 0; % Is the id in the list of neuron ids with motion artifacts? 
        
        % Experimental parameters
        d.conditions{i} = guide.conditions{r}; 
        d.framerate(i) = guide.framerate(r); 
        d.stim_on_frames(i, :) = eval(guide.stim_on_frames{r});  
        d.stim_off_frames(i, :) = eval(guide.stim_off_frames{r});
        d.intensity(i, :) = eval(guide.intensity{r});
        
        % Morphological and molecular information 
        d.gfp_pos_thresh(i) = guide.gfp_pos_thresh(r);
        d.gfp_neg_thresh(i) = guide.gfp_neg_thresh(r);
        d.fluorescent_marker_name{i} = n.getMarkerID();
        d.fluorescent_marker_intensity(i) = n.getMarkerIntensity(); 
        d.stratification{i} = n.getStratification(); 
        d.immuno_marker{i} = n.getFixedMarker(); 
        [subtype, subtype_code] = getSubtype(n, d(i,:)); 
        d.subtype{i} = subtype;        
        d.subtype_code(i, :) = subtype_code; 
        

        % Collect calcium imaging data, snr, and max z score for each movie 
        movies = n.getMovieList();
        for m = 1:n_movies
            
            % If the movie is available, enter the data
            if m <= length(movies)
                
                trace = n.getTrace(movies{m});
                d.(['movie_', num2str(m), '_data'])(i,1:length(trace)) = trace;
                
                max_responses = n.getTraceStat(movies{m}, 'amplitude');
                d.(['movie_', num2str(m), '_max_responses'])(i,1:length(max_responses)) = max_responses;
                
                d.(['movie_', num2str(m), '_snr'])(i) = ...
                    getSNR_v2(n, movies{m}, d.stim_on_frames(i, :), d.stim_off_frames(i, :));
                d.(['movie_', num2str(m), '_max_z_score'])(i) = ...
                    getZScore(n, movies{m}, d.stim_on_frames(i,:), d.framerate(i));
                
                % Determine if cell is light-sensitive in this movie 
                d.(['movie_', num2str(m), '_lr'])(i) = (d.(['movie_', num2str(m), '_max_z_score'])(i) > lr_params.z_score_thresh) && ... 
                                                       (d.(['movie_', num2str(m), '_snr'])(i) > lr_params.snr_thresh); 
                d.(['movie_', num2str(m), '_not_lr'])(i) = (d.(['movie_', num2str(m), '_max_z_score'])(i) < not_lr_params.z_score_thresh);
                
            end % Otherwise leave the fields as NaNs 
            
        end 
        
    end 

end 

% Narrow down an ipRGC's subtype from morphology and molecular markers
function [subtype, subtype_code] = getSubtype(neuron, d_row)

stratification = neuron.getStratification();

    switch stratification
        case 'OFF'
            
            subtype = 'M1';
            subtype_code = 1;
            
        case 'ON'
            
            if strcmp(neuron.getFixedMarker(), 'SMI32')
                
                subtype = 'M4';
                subtype_code = 3;
                
            else
                
                subtype = 'M2/5/6';
                subtype_code = 4;
                
            end
            
        case 'ON/OFF'
            
            subtype = 'M3';
            subtype_code = 2;
            
        otherwise
            
            if strcmp(neuron.getFixedMarker(), 'SMI32')
                
                subtype = 'M4';
                subtype_code = 3;
                
            elseif d_row.fluorescent_marker_intensity >= d_row.gfp_pos_thresh
                
                subtype = 'GFP+';
                subtype_code = 5;
                
            elseif d_row.fluorescent_marker_intensity < d_row.gfp_neg_thresh
                
                subtype = 'GFP-';
                subtype_code = 6;
                
            else
                
                subtype = 'Unknown';
                subtype_code = 7;
                
            end
            
    end
    
end


























