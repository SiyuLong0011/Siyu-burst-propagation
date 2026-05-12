function [sorted_pre, sorted_post, pca_results] = burst_onset_pca(ann1, ann2, n_components)
% BURST_ONSET_PCA Sort burst onsets within subjects and extract group-level patterns via PCA
%
% For each subject, sorts all channel onset times in ascending order
% (earliest to latest). Then performs PCA across subjects on the sorted
% onset distributions to extract group-level temporal patterns of burst
% onset ordering.
%
% Pre-event and post-event onsets are processed independently.
%
% Usage:
%   [sorted_pre, sorted_post, pca_results] = burst_onset_pca(ann1, ann2)
%   [sorted_pre, sorted_post, pca_results] = burst_onset_pca(ann1, ann2, n_components)
%
% Inputs:
%   ann1            - [channels x subjects] pre-event burst onset times in
%                     ms (negative values). 0 = no burst detected.
%                     (output from detect_burst_onset).
%   ann2            - [channels x subjects] post-event burst onset times in
%                     ms (positive values). 0 = no burst detected.
%   n_components    - (Optional) Number of PCA components to retain.
%                     Default: 10.
%
% Outputs:
%   sorted_pre      - [channels x subjects] pre-event onsets sorted within
%                     each subject (ascending order). NaN values (no burst)
%                     are placed at the end.
%   sorted_post     - [channels x subjects] post-event onsets sorted within
%                     each subject (ascending order). NaN values (no burst)
%                     are placed at the end.
%   pca_results     - Struct with fields:
%       .pre
%           .coeff          - Principal component coefficients (loadings).
%           .score          - PC scores for each subject.
%           .explained      - Percentage of variance explained.
%           .mu             - Mean used for centering.
%       .post               - Same structure as .pre for post-event.
%       .n_components       - Number of components retained.
%
% Notes:
%   - Onset values of 0 are treated as missing data (no burst detected)
%     and replaced with NaN before sorting.
%   - For PCA, NaN values are imputed with the column mean (mean onset at
%     that sorted rank across subjects) to allow matrix decomposition.
%   - Columns (sorted ranks) that are entirely NaN across all subjects are
%     removed before PCA and excluded from the output coefficients.
%
% Example:
%   [sorted_pre, sorted_post, pca_results] = burst_onset_pca(ann1, ann2);
%
% See also: detect_burst_onset, pca

    % --- Default parameters ---
    if nargin < 3 || isempty(n_components)
        n_components = 10;
    end

    [n_channels, n_subjects] = size(ann1);

    % --- Replace 0 (no burst) with NaN ---
    ann1(ann1 == 0) = NaN;
    ann2(ann2 == 0) = NaN;

    % --- Sort onsets within each subject (ascending) ---
    % NaN values are placed at the end by sort
    sorted_pre = zeros(n_channels, n_subjects);
    sorted_post = zeros(n_channels, n_subjects);

    for i = 1:n_subjects
        sorted_pre(:, i) = sort(ann1(:, i), 'ascend', 'MissingPlacement', 'last');
        sorted_post(:, i) = sort(ann2(:, i), 'ascend', 'MissingPlacement', 'last');
    end

    % --- Run PCA for each half ---
    pca_results.pre = run_pca_with_nan(sorted_pre, n_components);
    pca_results.post = run_pca_with_nan(sorted_post, n_components);
    pca_results.n_components = n_components;

    fprintf('PCA complete.\n');
    fprintf('  Pre-event:  top %d PCs explain %.1f%% variance.\n', ...
        length(pca_results.pre.explained), sum(pca_results.pre.explained));
    fprintf('  Post-event: top %d PCs explain %.1f%% variance.\n', ...
        length(pca_results.post.explained), sum(pca_results.post.explained));
end


function result = run_pca_with_nan(sorted_data, n_components)
% RUN_PCA_WITH_NAN Handle NaN imputation and all-NaN column removal before PCA

    % PCA input: [subjects x channels]
    data = sorted_data';

    % Remove columns that are entirely NaN (sorted ranks with no bursts)
    all_nan_cols = all(isnan(data), 1);
    data_clean = data(:, ~all_nan_cols);

    if isempty(data_clean)
        warning('All onset values are NaN. PCA cannot be performed.');
        result.coeff = [];
        result.score = [];
        result.explained = [];
        result.mu = [];
        result.removed_cols = find(all_nan_cols);
        return;
    end

    % Impute remaining NaN with column mean
    for c = 1:size(data_clean, 2)
        col = data_clean(:, c);
        nan_mask = isnan(col);
        if any(nan_mask)
            col(nan_mask) = mean(col, 'omitnan');
            data_clean(:, c) = col;
        end
    end

    % Run PCA
    n_comp = min(n_components, min(size(data_clean)) - 1);

    if n_comp < 1
        warning('Not enough data for PCA. Need at least 2 subjects or 2 channels.');
        result.coeff = [];
        result.score = [];
        result.explained = [];
        result.mu = [];
        result.removed_cols = find(all_nan_cols);
        return;
    end

    [coeff, score, ~, ~, explained, mu] = pca(data_clean, 'NumComponents', n_comp);

    result.coeff = coeff;
    result.score = score;
    result.explained = explained(1:n_comp);
    result.mu = mu;
    result.removed_cols = find(all_nan_cols);
end
