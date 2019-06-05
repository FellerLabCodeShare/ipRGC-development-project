% Main function for assigning cells to clusters based on their physiology 
function [d_result, pca_result, gmm, min_BICs, logBFs] = clusterCells(d, movie_prefix, data_selector, ...
    new_pca, new_gmms, pca_result, pca_params, gmm, gmm_params, qc_params)

    % Add columns for keeping track of cluster indexes and
    % within-cluster-order
    d.cluster_idx = NaN * ones(size(d, 1), 1);
    d.within_cluster_order = NaN * ones(size(d, 1), 1);

    % Extract neurons in d that pass quality control and are light-sensitive
    [d_cluster, d_misc] = qualityControl(d, data_selector, movie_prefix, qc_params); 
    
    % Run PCA 
    [X, d_cluster] = preProcess(d_cluster, movie_prefix);
    if new_pca
        [pca_result, X_feature] = runPCA(X, pca_params); 
        d_cluster.([movie_prefix, 'feature_vector']) = X_feature; 
    else
        X_feature = normalize(X*pca_result.B);
    end 
    
    % Plot results of PCA 
    if pca_params.plot_results; plotPCAResults(X_feature, pca_result); end
    
    % Cluster the data with a GMM 
    if new_gmms
        [gmm, min_BICs, logBFs] = trainGMM(X_feature, gmm_params);
    else
        min_BICs = gmm.BIC;
        logBFs = NaN;
    end 
    d_cluster = clusterData(d_cluster, X_feature, gmm); 
    
    % Sort data within clusters 
    d_cluster_sort = sortWithinCluster(d_cluster, movie_prefix); 
    
    % Re-combine sorted data with the larger dataset 
    d_result = vertcat(d_cluster_sort, d_misc); 

end

%% Split dataset into portion that will be clustered and portion that won't
function [d_cluster, d_misc] = qualityControl(d, data_selector, movie_prefix, qc_params)
    
    % Reject cells with a significant motion artifact that occurred during
    % a stimulus 
    motion_artifact_selector = ~d.motion_artifact; 

    % Determine which cells have sufficiently high Z score. 
    %z_selector = strcmp(d.subtype, 'M1') | d.([movie_prefix,
    %'max_z_score']) > qc_params.z_score_thresh; %  exempts M1s 
    z_selector = d.([movie_prefix, 'max_z_score']) > qc_params.z_score_thresh; 
    
    % Determine which cells have sufficiently high SNR 
    snr_selector = d.([movie_prefix, 'snr']) > qc_params.snr_thresh;
    
    % Select cells for clustering
    cluster_selector = z_selector & snr_selector & motion_artifact_selector & data_selector;
    d_cluster = d(cluster_selector, :); 
    
    % And keep the cells that are not selected in a separate dataset while
    % clustering proceeds
    d_misc = d(~cluster_selector, :); 
    
    if qc_params.plot
        
       figure
       plot(d.([movie_prefix, 'max_z_score'])(cluster_selector, :), d.([movie_prefix, 'snr'])(cluster_selector, :), '.k');
       hold on 
       plot(d.([movie_prefix, 'max_z_score'])(~cluster_selector & data_selector, :), d.([movie_prefix, 'snr'])(~cluster_selector & data_selector, :), '.r');
       hold off 
       xlabel('Max. Z score');
       ylabel('SNR');
       set(gcf, 'Position', [190 400 290 260]);
       pbaspect([1 1 1]); 
       
    end    

end

%% Extract and preprocess the data for PCA 
function [X, d_cluster] = preProcess(d_cluster, movie_prefix) 

    % Take out the dataset that will be fed in to the clustering algorithm
    X = d_cluster.([movie_prefix, 'data']);
    
    % Highpass filter the data 
    fpass = 0.01; 
    X = (highpass(X', fpass, d_cluster.framerate(1)))';  
    
    % Select a subset of the data that consists of the stimulus and an
    % equal amount of time following the stimulus (to catch the decay) 
    stim_and_decay_frames = getStimAndDecayFrames(d_cluster);
    X = X(:, stim_and_decay_frames);
    
    % Normalize each trace to its own max response
    X = X ./ max(X, [], 2);
    
    % Plot a heatmap of when the stimulus was on for the preprocessed data 
    if false
        
        figure
        stim_on_frames = getStimFrames(d_cluster);
        imagesc(stim_on_frames(:, stim_and_decay_frames)); 
        title('stimulus frames for data used in clustering');
        
        figure
        imagesc(stim_on_frames); 
        title('stimulus frames');
    
    end
    
    % Store the preprocessed data
    d_cluster.([movie_prefix, 'clustering_data'])(:, 1:size(X,2)) = X;
    
    % Center and normalize X such that the column means are 0 and
    % the column euclidean lengths are 1.
    X = normalize(X);
    
end

%% Run Sparse Principal Components Analysis on the preprocessed data 
function [pca_result, X_feature_normalized] = runPCA(X, pca_params)
    
    % Run the sparse pca algorithm. Last argument is verbose mode. 
    [B, SD, L, D, Paths] = spca(X, pca_params.Gram, pca_params.k, ...
        pca_params.delta, pca_params.stop, pca_params.max_steps, pca_params.convergence_criterion, true); 

    % Create sPCA results struct 
    pca_result.B = B;
    pca_result.SD = SD;
    pca_result.L = L; 
    pca_result.D = D;
    pca_result.Paths = Paths; 
    
    % Make a k-dimensional feature vector for each cell. X is a (cells x
    % frames) matrix. B is a (frames x k) matrix. 
    X_feature = X*B; 

    % Standardize each feature across cells 
    X_feature_normalized = normalize(X_feature);

end

%% Select and train a Gaussian Mixture Model 
function [gmm, min_BICs, log_BFs] = trainGMM(X_feature, gmm_params)

    % Train models with different numbers of Gaussians and keep track of
    % their BICs and log likelihood. 
    Models = cell(gmm_params.max_n_clusters, gmm_params.n_iterations); 
    BICs = NaN * ones(gmm_params.max_n_clusters, gmm_params.n_iterations); 
    neg_log_Likely = NaN * ones(gmm_params.max_n_clusters, gmm_params.n_iterations); 
    
    for c = 1:gmm_params.max_n_clusters
        
        for i = 1:gmm_params.n_iterations 
            
            options = struct();
            options.MaxIter = 300; % default 100
            
            Models{c, i} = fitgmdist(X_feature, c,...
            'CovarianceType','diagonal','RegularizationValue',10^-5, 'Options', options);
        
            if Models{c, i}.Converged % If the model converged 
                
                BICs(c, i) = Models{c, i}.BIC;
                neg_log_Likely(c, i) = Models{c, i}.NegativeLogLikelihood; 
                
            end 
            
        end 
        
    end 

    % Select the model with the lowest BIC
    [~, min_lin_idx] = min(BICs(:)); 
    gmm = Models{min_lin_idx};
    
    if gmm_params.plot_results
        
        % Plot BIC as a function of # of clusters 
        figure
        hold on 
        plot(BICs, 'o', 'Color', [0.5 0.5 0.5], 'MarkerSize', 1);
        min_BICs = min(BICs, [], 2);
        
        % Normalized BIC 
        %plot(2 - min_BICs ./ min(min(min_BICs)), 'k', 'LineWidth', 2); 
        
        plot(min_BICs, 'k', 'LineWidth', 2); 
        xlabel('# of clusters'); 
        ylabel('Normalized BIC'); 
        xticks(1:11); 
        
        %title('Optimization of cluster number'); 
        
        % Add the chosen model as a red dot 
        plot(gmm.NumComponents, gmm.BIC, '.r');
        hold off
        set(gcf, 'Position', [190 400 230 200]);
        pbaspect([1 1 1]);
        
        % Calculate and plot Bayes factor for n+1 clusters 
        log_BFs = NaN * ones(1, gmm_params.max_n_clusters - 1);
        for c = 1:(gmm_params.max_n_clusters - 1)
            
            log_BFs(c) = (min_BICs(c) - min_BICs(c + 1)) / 2; 
        
        end
        
        figure
        plot(1:(gmm_params.max_n_clusters - 1), log_BFs, '.-k', 'MarkerSize', 20);
        xlabel('Cluster number'); 
        ylabel('Bayes Factor with respect to n+1 cluster model'); 
        hline(6); % Bayes Factor of >6 is strong evidence for further splittinig 
        set(gcf, 'Position', [190 400 230 200]);
        pbaspect([1 1 1]);
        
    end 

end

%% Use the GMM to assign cluster indexes and cluster probabilities to cells
function [d_cluster] = clusterData(d_cluster, X_feature, gmm)
    
    d_cluster.cluster_idx = gmm.cluster(X_feature); 
    
end

%% Sort cells within cluster by assigning them a number ranking when their first light response occurred or by subtype 
function [d_cluster] = sortWithinCluster(d_cluster, movie_prefix)
    
    cluster_idxs = unique(d_cluster.cluster_idx);

    % For each cluster number
    for i = 1:length(cluster_idxs)
        
        % Make a selection array for this cluster
        selector = d_cluster.cluster_idx == cluster_idxs(i);
    
        % Pick out the data 
        cluster_data = d_cluster.([movie_prefix, 'data'])(selector, :); 
        
        % Get the stimulus and decay frames
        stim_and_decay_frames = getStimAndDecayFrames(d_cluster);
        
        % Normalize data for each cell to maximum value
        cluster_data_norm = cluster_data ./ max(cluster_data, [], 2); 
        
        % Determine which data points are over a threshold 
        data_over_thresh = cluster_data(:, stim_and_decay_frames) > 0.3; 
        
        % Find out when the first responses occurred for each cell
        [~, first_response_idx] = max(data_over_thresh, [], 2); 
        
        % Make a sorting vector that orders cells by the time of their
        % first response (from earliest to latest) 
        [~, max_response_sort] = sort(first_response_idx); 
        
        % Make a sorting vector that orders cells by their subtype
        [~, subtype_sort] = sort(d_cluster.subtype_code(selector));
        
        % Use the sorting vector to make a vector that describes what the
        % ordering of the cell should be within the cluster 
        %ordering = orderingFromSorting(max_response_sort); % By max response timing 
        ordering = orderingFromSorting(subtype_sort); % By subtype 
        
        % Store that sorting vector in the within-cluster-order column
        d_cluster.within_cluster_order(selector) = ordering; 
        
    end 

end


% Generates several plots to illustrate the results of sparse Principal
% Components Analysis 
function [] = plotPCAResults(X_feature, pca_result) 
        
        % Display adjusted variance (e.g. percentage of total variance explained)
        % for the sparse prinicipal components
        figure
        plot(pca_result.SD/sum(pca_result.D) * 100, '.-k', 'MarkerSize', 15); 
        %title('Variance explained by sparse principal components');
        xlabel('Principal component #');
        ylabel('% of variance explained'); 
        set(gcf, 'Position', [190 400 230 200]);
        pbaspect([1 1 1]); 

        % Display the loading vectors 
        figure
        imagesc(pca_result.B'); 
        %title(['Loading vectors for ', num2str(pca_params.k), ' principal components']); 
        xlabel('Frames');
        ylabel('Principal Components'); 
        set(gcf, 'Position', [190 400 230 200]);
        pbaspect([1 1 1]);
        
        % Apply a custom colormap where zero is white (to emphasize
        % sparseness of PCA) 
        caxis([-0.2 0.2]); 
        custom_parula = parula; 
        custom_parula(33, :) = [1 1 1];
        colormap(custom_parula); 
        
        
        % Display the feature vectors for each cell 
        % first as a heatmap 
        figure
        imagesc(sortrows(X_feature));
        %title('Feature vectors for each cell'); 
        xlabel('Features');
        ylabel('Cells'); 
        set(gcf, 'Position', [190 400 230 200]);
        pbaspect([1 1 1]);
        
        % And then as a plot of first two PCs 
        figure
        plot(X_feature(:, 1), X_feature(:, 2), '.k', 'MarkerSize', 10);
        title('First two PCs (norm.)'); 
        xlabel('PC 1');
        ylabel('PC 2'); 
        set(gcf, 'Position', [190 400 230 200]);
        pbaspect([1 1 1]);
        
end 
    




