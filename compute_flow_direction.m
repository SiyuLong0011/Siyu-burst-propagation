function subject_region_means = compute_flow_direction(source_data, onsets, window_length, event_sample)
% COMPUTE_FLOW_DIRECTION Compute mean optical flow direction per channel
%
% For each channel, extracts the optical flow direction angles within a
% time window starting at that channel's burst onset, computes the circular
% mean, and projects the result onto the range [0, pi] (ignoring flow
% directionality, treating opposite directions as equivalent).
%
% Usage:
%   subject_region_means = compute_flow_direction(source_data, onsets)
%   subject_region_means = compute_flow_direction(source_data, onsets, window_length, event_sample)
%
% Inputs:
%   source_data     - [channels x timepoints] optical flow direction matrix.
%                     Values are angles in radians.
%   onsets          - [channels x 1] burst onset times in ms relative to
%                     event onset (output from detect_burst_onset). Channels
%                     with onset = 0 (no burst) are skipped and assigned NaN.
%   window_length   - (Optional) Duration of the analysis window after onset
%                     in ms (= samples at 1 kHz). Default: 100.
%   event_sample    - (Optional) Sample index corresponding to time 0 ms in
%                     source_data. Default: 1401 (for 2800-sample data with
%                     -1400 to +1400 ms range).
%
% Outputs:
%   subject_region_means - [channels x 1] mean flow direction per channel
%                          in radians, projected to [0, pi]. NaN for channels
%                          with no burst detected.
%
% Algorithm:
%   For each channel:
%     1. Convert onset time (ms) to sample index in source_data.
%     2. Extract angles from onset to onset + window_length.
%     3. Double the angles (theta -> 2*theta) to map axial data to circular.
%     4. Compute circular mean of doubled angles using atan2(sum(sin), sum(cos)).
%     5. Halve the result to project back to [0, pi].
%
% Notes:
%   - The doubling-and-halving trick converts axial (undirected) data to
%     circular data for proper averaging. This ensures that angles near 0
%     and pi are treated as similar directions.
%   - Channels with onset = 0 are treated as missing (no burst detected).
%   - If the window extends beyond the data, it is truncated.
%
% Example:
%   load('subject01_optical_flow.mat');  % loads Source_Data [200 x 2800]
%   load('subject01_onsets.mat');        % loads onset_pre [200 x 1]
%   means = compute_flow_direction(Source_Data, onset_pre, 100);
%
% See also: detect_burst_onset

    % --- Default parameters ---
    if nargin < 3 || isempty(window_length)
        window_length = 100;  % ms (= samples at 1 kHz)
    end
    if nargin < 4 || isempty(event_sample)
        event_sample = 1401;
    end

    % --- Get dimensions ---
    [n_channels, n_timepoints] = size(source_data);
    subject_region_means = NaN(n_channels, 1);

    % --- Compute mean direction per channel ---
    for ch = 1:n_channels
        % Skip channels with no burst
        if onsets(ch) == 0
            continue;
        end

        % Convert onset time (ms relative to event) to sample index
        onset_sample = onsets(ch) + event_sample;

        % Define window
        win_start = onset_sample;
        win_end = min(onset_sample + window_length - 1, n_timepoints);

        if win_start < 1 || win_start > n_timepoints
            continue;
        end

        % Extract angles in the window
        angles = source_data(ch, win_start:win_end);

        % --- Circular mean projected to [0, pi] ---
        % Double angles to convert axial data to circular
        doubled = 2 * angles;
        cum_x = sum(cos(doubled));
        cum_y = sum(sin(doubled));
        mean_doubled = atan2(cum_y, cum_x);

        % Halve to project back to [0, pi]
        mean_angle = mean_doubled / 2;

        % Ensure result is in [0, pi]
        if mean_angle < 0
            mean_angle = mean_angle + pi;
        end

        subject_region_means(ch) = mean_angle;
    end
end
