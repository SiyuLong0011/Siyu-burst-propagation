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
%                     (output from batch_detect_burst_onset or detect_burst_onset).
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
%           .coeff          - [channels x n_components] principal component
%                             coefficients (loadings) for pre-event.
%           .score          - [subjects x n_components] PC scores for pre-event.
%           .explained      - [n_components x 1] percentage of variance
%                             explained by each component.
%           .mu             - [1 x channels] mean used for centering.
%       .post
%           .coeff          - Same as above for post-event.
%           .score          - Same as above for post-event.
%           .explained      - Same as above for post-event.
%           .mu             - Same as above for post-event.
%       .n_components   - Number of components retained.
%
% Notes:
%   - Onset values of 0 are treated as missing data (no burst detected)
%     and replaced with NaN before sorting.
%   - For PCA, NaN values are imputed with the column mean (mean onset at
%     that sorted rank across subjects) to allow matrix decomposition.
%     The original NaN positions are preserved in sorted_pre/sorted_post.
%   - PCA is performed on [subjects x channels], so each row is one
%     subject's sorted onset distribution and each component captures a
%     pattern of onset ordering shared across subjects.
%
% Example:
%   [ann1, ann2] = batch_detect_burst_onset('/data/burst_results/');
%   [sorted_pre, sorted_post, pca_results] = burst_onset_pca(ann1, ann2);
%
%   % Plot variance explained
%   figure;
%   subplot(1,2,1); bar(pca_results.pre.explained);
%   title('Pre-event variance explained'); xlabel('PC'); ylabel('%');
%   subplot(1,2,2); bar(pca_results.post.explained);
%   title('Post-event variance explained'); xlabel('PC'); ylabel('%');
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

    % --- Prepare data for PCA ---
    % PCA input: [subjects x channels] (each row = one subject's sorted onsets)
    data_pre = sorted_pre';   % [subjects x channels]
    data_post = sorted_post'; % [subjects x channels]

    % Impute NaN with column mean for PCA (mean onset at each sorted rank)
    data_pre_imputed = impute_nan_colmean(data_pre);
    data_post_imputed = impute_nan_colmean(data_post);

    % --- Run PCA ---
    n_comp_pre = min(n_components, min(size(data_pre_imputed)) - 1);
    n_comp_post = min(n_components, min(size(data_post_imputed)) - 1);

    [coeff_pre, score_pre, ~, ~, explained_pre, mu_pre] = pca(data_pre_imputed, ...
        'NumComponents', n_comp_pre);
    [coeff_post, score_post, ~, ~, explained_post, mu_post] = pca(data_post_imputed, ...
        'NumComponents', n_comp_post);

    % --- Assemble output ---
    pca_results.pre.coeff = coeff_pre;
    pca_results.pre.score = score_pre;
    pca_results.pre.explained = explained_pre(1:n_comp_pre);
    pca_results.pre.mu = mu_pre;

    pca_results.post.coeff = coeff_post;
    pca_results.post.score = score_post;
    pca_results.post.explained = explained_post(1:n_comp_post);
    pca_results.post.mu = mu_post;

    pca_results.n_components = n_components;

    fprintf('PCA complete.\n');
    fprintf('  Pre-event:  top %d PCs explain %.1f%% variance.\n', ...
        n_comp_pre, sum(explained_pre(1:n_comp_pre)));
    fprintf('  Post-event: top %d PCs explain %.1f%% variance.\n', ...
        n_comp_post, sum(explained_post(1:n_comp_post)));
end


function data_out = impute_nan_colmean(data_in)
% IMPUTE_NAN_COLMEAN Replace NaN values with column mean
%
% For each column, NaN values are replaced with the mean of the non-NaN
% values in that column. If an entire column is NaN, it remains NaN.

    data_out = data_in;
    for c = 1:size(data_in, 2)
        col = data_in(:, c);
        nan_mask = isnan(col);
        if any(nan_mask) && ~all(nan_mask)
            col(nan_mask) = nanmean(col);
            data_out(:, c) = col;
        end
    end
end
