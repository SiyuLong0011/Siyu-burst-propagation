function beta_burst = detect_beta_bursts(tf_avg, ap_avg, percentile_thresh, min_burst_duration)
% DETECT_BETA_BURSTS Detect beta bursts from TF map after removing aperiodic component
%
% Subtracts the aperiodic (1/f) spectral component from the trial-averaged
% time-frequency map, then identifies beta bursts by thresholding at the
% 75th percentile and retaining only sustained bursts that exceed a minimum
% duration. Pre-event and post-event halves are processed independently.
%
% Usage:
%   beta_burst = detect_beta_bursts(tf_avg, ap_avg)
%   beta_burst = detect_beta_bursts(tf_avg, ap_avg, percentile_thresh, min_burst_duration)
%
% Inputs:
%   tf_avg              - [channels x timepoints] trial-averaged time-frequency
%                         map (output from extract_beta_tf.m). Assumed to be
%                         evenly split into pre-event and post-event halves.
%   ap_avg              - [channels x 2] trial-averaged aperiodic model from
%                         SPRiNT. Column 1 corresponds to the pre-event half,
%                         column 2 to the post-event half.
%   percentile_thresh   - (Optional) Percentile threshold for burst detection.
%                         Applied independently to each channel and each half.
%                         Default: 75.
%   min_burst_duration  - (Optional) Minimum number of consecutive samples
%                         above threshold to be considered a burst.
%                         Default: 100.
%
% Outputs:
%   beta_burst          - [channels x timepoints] binary matrix. 1 indicates
%                         a timepoint belonging to a detected beta burst,
%                         0 otherwise.
%
% Algorithm:
%   For each half (pre-event, post-event) independently:
%     1. Subtract the corresponding aperiodic component from each channel.
%     2. Compute the 75th percentile of the corrected signal per channel.
%     3. Binarize: timepoints >= threshold are set to 1.
%     4. Label contiguous above-threshold regions (connected components).
%     5. Remove bursts shorter than min_burst_duration samples.
%
% Example:
%   load('subject01_tf.mat');  % loads tf_avg [channels x 2800] and ap_avg [channels x 2]
%   beta_burst = detect_beta_bursts(tf_avg, ap_avg);
%
% See also: extract_beta_tf, bwlabel, regionprops

    % --- Default parameters ---
    if nargin < 3 || isempty(percentile_thresh)
        percentile_thresh = 75;
    end
    if nargin < 4 || isempty(min_burst_duration)
        min_burst_duration = 100;  % samples
    end

    % --- Get dimensions ---
    [n_channels, n_timepoints] = size(tf_avg);
    half_point = floor(n_timepoints / 2);

    % --- Subtract aperiodic component ---
    tf_corrected = zeros(n_channels, n_timepoints);
    tf_corrected(:, 1:half_point) = tf_avg(:, 1:half_point) - ap_avg(:, 1);
    tf_corrected(:, (half_point + 1):end) = tf_avg(:, (half_point + 1):end) - ap_avg(:, 2);

    % --- Preallocate output ---
    beta_burst = zeros(n_channels, n_timepoints);

    % --- Detect bursts (parallelized across channels) ---
    if isempty(gcp('nocreate'))
        parpool('local');
    end

    parfor i = 1:n_channels
        temp_result = zeros(1, n_timepoints);

        % --- Pre-event half ---
        x1 = tf_corrected(i, 1:half_point);
        thresh1 = prctile(x1, percentile_thresh);
        bin_x1 = x1 >= thresh1;

        [L1, num1] = bwlabel(bin_x1);
        if num1 > 0
            stats1 = regionprops(L1, 'Area', 'PixelIdxList');
            for j = 1:num1
                if stats1(j).Area >= min_burst_duration
                    temp_result(stats1(j).PixelIdxList) = 1;
                end
            end
        end

        % --- Post-event half ---
        x2 = tf_corrected(i, (half_point + 1):end);
        thresh2 = prctile(x2, percentile_thresh);
        bin_x2 = x2 >= thresh2;

        [L2, num2] = bwlabel(bin_x2);
        if num2 > 0
            stats2 = regionprops(L2, 'Area', 'PixelIdxList');
            for j = 1:num2
                if stats2(j).Area >= min_burst_duration
                    temp_result(half_point + stats2(j).PixelIdxList) = 1;
                end
            end
        end

        beta_burst(i, :) = temp_result;
    end

    fprintf('Burst detection complete: %d channels x %d timepoints.\n', ...
        n_channels, n_timepoints);
end
