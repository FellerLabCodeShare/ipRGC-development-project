% Build a table of calcium imaging data (neurons x observations) 

%% Parameters 

% Parameters for loading data  
home = false; % Is this my PC at home (vs. in lab)? 
if home
    
    root_directory = 'C:\Google Drive\data on cloud\imaging'; % Home PC 

else 
    
    root_directory = 'C:\Users\Feller Lab\Google Drive\data on cloud\imaging'; % Work PC 
    
end 
guide_name = 'Opn4 eGFP interline imaging - updated gfp thresholds.xlsx';
experiment_type = 'intensity-response';
n_movies = 2;
n_stims = 10; 
n_movie_frames = 1370;
n_cluster_frames = 924; % Number of movie frames used as input to the clustering algorithm 

% General clustering parameters 
new_pca = true; % run new sPCA? 
new_models = true; % train new GMMs? 
data_column_name = 'movie_1_'; % Specify the data column name prefix the clustering algorithm will operate on 

% Parameters for sPCA 
pca_params.Gram = []; % Gram matrix. Optional. I don't understand what this is for yet.
pca_params.k = 10; % Desired number of sparse principal components
pca_params.delta = inf; % infinity for soft thresholding, if p >> n (I do have more frames than cells)
pca_params.stop = -300; % If stop is negative, stop refers to integer # of desired nonzero variables. If stop is positive, the upper bound on L1 norm of Beta coefficients. was 300
pca_params.convergence_criterion = 1e-9; % Default
pca_params.max_steps = 300; % Default is 300 
pca_params.plot_results = true; 

% Parameters for GMM 
gmm_params.max_n_clusters = 11; 
gmm_params.n_iterations = 500; % was 500
gmm_params.plot_results = true; 

% Parameters for quality control - e.g. only do clustering on these cells 
qc_params.snr_thresh = 5;
qc_params.z_score_thresh = 6; 
qc_params.plot = true;

% Parameters for calling a cell light sensitive or definitively not light
% sensitive
lr_params.snr_thresh = 5;
lr_params.z_score_thresh = 6; 
not_lr_params.z_score_thresh = 5; % must be less than this threshold 

% Parameters for plot aesthetics 
subtype_names = {'M1', 'M3', 'M4', 'M2/5/6', 'GFP+', 'GFP-', 'Unknown'};
m1_color = [0.2235, 0.3255, 0.6353]; % RGB triplet (illustrator rgb / 255) 
m3_color = [0.4157, 0.7373, 0.2706]; 
m4_color = [0.8706, 0.2706, 0.5882];
m2_5_6_color = [0.6863, 0.5843, 0.1922];
gfp_pos_color = [0.5, 0.5, 0.5];
gfp_neg_color = [0, 0, 0]; 
unknown_color = [1, 1, 1];
subtype_colormap = [m1_color; m3_color; m4_color;  m2_5_6_color;  gfp_pos_color; gfp_neg_color; unknown_color];
subtype_colors = {m1_color, m3_color, m4_color, m2_5_6_color, gfp_pos_color,  gfp_neg_color, unknown_color};
subtype_color_lookup = containers.Map(subtype_names, subtype_colors); 
trace_plot_type = 'heatmaps'; % Can be 'heatmaps' or 'mean traces'
second_plot_type = 'subtype color code'; % Can be one of 'subtype color code', 'subtype pie', 'mouse color code', 'mean trace', 'irradiance-response'
sort_traces = true; % Sort traces according to the earliest supra-threshold light response? 
scatter_xylimits = [0 5];


%% Create the table and cluster cells 

% Build the table 
d_original = buildCalciumImagingTable(root_directory, guide_name, experiment_type, n_movies, n_stims, n_movie_frames, n_cluster_frames, pca_params.k, lr_params, not_lr_params);

% Select a subset of the data rows for clustering

% Don't include knockouts or the mouse for which there were some
% glutamatergic waves
d_original = d_original(~strcmp(d_original.genotype, 'Cx30.2 lacZ/lacZ; Opn4 eGFP') & d_original.date ~= 170814, :);

% Select data for clustering 
data_selector = true(size(d_original, 1), 1); % Include everything 

% Cluster the cells that, in control conditions, are light responsive and
% pass quality control. Can pass in existing PCA and GMM, or leave the
% arguments as empty arrays if new_models = false - the function will run a
% new PCA, then train and select a new GMM. 

% Version where I run PCA and train a new GMM 
[d, pca_results, gmm, min_BICs, logBFs] = clusterCells(d_original, data_column_name, data_selector, new_pca, new_models, [], pca_params, [], gmm_params, qc_params); 

% Version where I provide a gmm and/or pca 
%[d, pca_results, gmm, min_BICs, logBFs] = clusterCells(d_original, data_column_name, data_selector, new_pca, new_models, pca_results, pca_params, [], gmm_params, qc_params); 

% Make a color-code for the clusters according to where the majority of
% subtypes end up. Leave out the unknown cells and GFP+ cells 
cluster_colormap = makeClusterColormap(d, subtype_names([1:4, end-1]), subtype_colormap([1:4, end-1], :));

%% Plot a bar chart of the percentage of GFP+ cells that are light-responsive
if false
    
    [h, C] = gfpLRBarChart(d(strcmp(d.conditions, 'SCH'), :), 'movie_1');
        
end

%% Visualize clustered data 
mouse_colormap = distinguishable_colors(length(unique(d.mouse_code)));


% Plot clustering result on PC plot 
cluster_edge_colormap = cluster_colormap;
cluster_edge_colormap(1, :) = [0 0 0]; % Make f4b have a black ring 

if true
   
   d_cluster = d(~isnan(d.cluster_idx), :); 
   cluster_idxs = unique(d_cluster.cluster_idx); 
    
   figure 
   for i = 1:length(cluster_idxs)
        selector = d_cluster.cluster_idx == i;
%         plot3(d.movie_1_feature_vector(selector,1), d.movie_1_feature_vector(selector,2), d.movie_1_feature_vector(selector, 3), ...
%             'o', 'MarkerFaceColor', cluster_colormap(i, :), 'MarkerEdgeColor', cluster_edge_colormap(i, :), 'MarkerSize', 5);
        plot3(d.movie_1_feature_vector(selector,1), d.movie_1_feature_vector(selector,2), d.movie_1_feature_vector(selector, 3), ...
            'o', 'MarkerFaceColor', cluster_colormap(i, :), 'MarkerEdgeColor', 'k', 'MarkerSize', 5);
        hold on 
        
   end
   xlabel('PC 1');
   ylabel('PC 2'); 
   zlabel('PC 3'); 
   hold off
   set(gcf, 'Position', [190 400 230 200]);
   pbaspect([1 1 1]); 
   xlim([-0.25 0.25]); 
   ylim([-0.25 0.25]); 
   zlim([-0.25 0.25]); 
   xticks([-0.25 0 0.25]);
   yticks([-0.25 0 0.25]);
   zticks([-0.25 0 0.25]);
   grid on 
   set(gca, 'Projection','perspective')
      
end


% Plot heatmaps and subtype for...

% Full dataset in control conditions 
if true 
    
    cluster_heatmaps = plotClusterHeatmaps(d, subtype_colormap, mouse_colormap, ...
        cluster_colormap, 'movie_1_data', 'movie_1_max_responses', 'movie_1', trace_plot_type, second_plot_type, sort_traces, 'within_cluster_order', 'cluster_idx');

end 

% Cells before and after MFA 
if true 
    
    data_selector = strcmp(d.conditions, 'MFA'); % Pick neurons to include in the figure
    cluster_heatmaps = plotClusterHeatmaps(d(data_selector, :), subtype_colormap,...
        mouse_colormap, cluster_colormap, 'movie_1_data', 'movie_1_max_responses', 'movie_1', trace_plot_type, second_plot_type, sort_traces, 'within_cluster_order', 'cluster_idx');
    
    cluster_heatmaps = plotClusterHeatmaps(d(data_selector, :), subtype_colormap,...
        mouse_colormap, cluster_colormap, 'movie_2_data', 'movie_2_max_responses', 'movie_1', trace_plot_type, second_plot_type, sort_traces, 'within_cluster_order', 'cluster_idx');

    % Scatter plots of max amplitude before and after the drug 
    dscatter = d((d.movie_1_lr | d.movie_2_lr) & strcmp(d.conditions, 'MFA') & ~d.motion_artifact, :); 
    open_selector = logical(zeros(1, length(dscatter.movie_1_lr)));
    edge_colormap = cluster_colormap; 
    edge_colormap(6, :) = [0 0 0]; 
    amplitude_scatter = plotAmplitudeScatter(dscatter, cluster_colormap, edge_colormap, 'cluster_idx',...
        'movie_1_max_responses', 'movie_2_max_responses', open_selector, scatter_xylimits);
    pbaspect([1 1 1]); 
    
    % Plot of fold-change in number of light-responsive cells by mouse 
    fold_change_plot = plotFoldChangeInNumberLRCells(dscatter, 'movie_1', 'movie_2'); 
    
    % Plot percent of cells of each subtype that are light-responsive in
    % the two conditions
    percentLRBySubtype(d(strcmp(d.conditions, 'MFA') & ~d.motion_artifact, :), subtype_colormap);
    
    % Plot the percentage of GFP+ cells that are light-responsive for each
    % mouse in MFA. C is a mouse x class array where the columns are lr
    % cells, non-responsive cells, and ambiguous cells 
    [h, C] = gfpLRBarChart(d(strcmp(d.conditions, 'MFA'), :), 'movie_2');
        
end 

% Cells before and after SCH 
if true 
    
    data_selector = strcmp(d.conditions, 'SCH') & d.movie_1_lr & d.movie_2_lr; % Pick neurons to include in the figure
    cluster_heatmaps = plotClusterHeatmaps(d(data_selector, :), subtype_colormap,...
        mouse_colormap, cluster_colormap, 'movie_1_data', 'movie_1_max_responses', 'movie_1', trace_plot_type, second_plot_type, sort_traces, 'within_cluster_order', 'cluster_idx');
    
    cluster_heatmaps = plotClusterHeatmaps(d(data_selector, :), subtype_colormap,...
        mouse_colormap, cluster_colormap, 'movie_2_data', 'movie_2_max_responses', 'movie_1', trace_plot_type, second_plot_type, sort_traces, 'within_cluster_order', 'cluster_idx');

    % Scatter plots of max amplitude before and after the drug 
    dscatter = d((d.movie_1_lr | d.movie_2_lr) & strcmp(d.conditions, 'SCH') & ~d.motion_artifact, :); 
    %open_selector = ~dscatter.movie_1_lr & dscatter.movie_2_lr;
    open_selector = logical(zeros(1, length(dscatter.movie_1_lr)));
    edge_colormap = cluster_colormap; 
    edge_colormap(6, :) = [0 0 0]; 
    amplitude_scatter = plotAmplitudeScatter(dscatter, cluster_colormap, edge_colormap, 'cluster_idx',...
        'movie_1_max_responses', 'movie_2_max_responses', open_selector, scatter_xylimits);
    pbaspect([1 1 1]); 
    
    % Plot of fold-change in number of light-responsive cells by mouse 
    fold_change_plot = plotFoldChangeInNumberLRCells(dscatter, 'movie_1', 'movie_2');
    
    % Plot percent of cells of each subtype that are light-responsive in
    % the two conditions
    percentLRBySubtype(d(strcmp(d.conditions, 'SCH') & ~d.motion_artifact, :), subtype_colormap);
    
    % Plot the percentage of GFP+ cells that are light-responsive for each
    % mouse in SCH. C is a mouse x class array where the columns are lr
    % cells, non-responsive cells, and ambiguous cells 
    [h, C] = gfpLRBarChart(d(strcmp(d.conditions, 'SCH'), :), 'movie_2');      

end 

% Cells that gained light responses in SCH
if true
    
    % Pick neurons to include in the figure
    data_selector = strcmp(d.conditions, 'SCH') & d.movie_1_max_z_score < not_lr_params.z_score_thresh & d.movie_2_lr; 
    
    % Sort the rows of the table according to which cells respond earliest
    % in the movie
    d_gained_sch = sortWithinCategory(d(data_selector, :), 'subtype_code', 'within_subtype_order', 'movie_2'); 
    
    % Make the heatmaps 
    cluster_heatmaps = plotClusterHeatmaps(d_gained_sch, subtype_colormap, ...
        mouse_colormap, cluster_colormap, 'movie_1_data', 'movie_1_max_responses', ...
        'movie_2', trace_plot_type, second_plot_type, true, 'within_subtype_order', 'subtype_code');
    
    cluster_heatmaps = plotClusterHeatmaps(d_gained_sch, subtype_colormap,...
        mouse_colormap, cluster_colormap, 'movie_2_data', 'movie_2_max_responses', ...
        'movie_2', trace_plot_type, second_plot_type, true, 'within_subtype_order', 'subtype_code');

end 


%% Scatter plot of maximum response amplitudes in different conditions 
if false 
    xylim  = [0 6]; 
    dscatter = d(d.movie_1_lr | d.movie_2_lr, :); % Cells that are light sensitive in some condition 
    figure 
    hold on 
    plot(max(dscatter.movie_1_max_responses(strcmp(dscatter.conditions, 'MFA'), :), [], 2), ...
        max(dscatter.movie_2_max_responses(strcmp(dscatter.conditions, 'MFA'), :), [], 2), '.r', 'MarkerSize', 8);

    plot(max(dscatter.movie_1_max_responses(strcmp(dscatter.conditions, 'SCH'), :), [], 2), ...
        max(dscatter.movie_2_max_responses(strcmp(dscatter.conditions, 'SCH'), :), [], 2), '.g', 'MarkerSize', 8);

    plot(xylim, xylim, '-k');
    xlim(xylim);
    ylim(xylim);
    hold off

    set(gcf, 'Position', [0 760, 240, 200]);
end 










