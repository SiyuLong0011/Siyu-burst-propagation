function [tf_avg, event_onsets] = extract_beta_tf(data, fs, freq_range, epoch_window, min_event_interval, morlet_cycles)
% EXTRACT_BETA_TF Extract trial-averaged beta-band time-frequency representation
%
% Computes a Morlet wavelet-based time-frequency decomposition in the beta
% band for task-related MEG data. Events (trial onsets) are automatically
% detected from a trigger channel, and the resulting single-trial TF maps
% are averaged across trials.
%
% Usage:
%   [tf_avg, event_onsets] = extract_beta_tf(data, fs)
%   [tf_avg, event_onsets] = extract_beta_tf(data, fs, freq_range, epoch_window, min_event_interval, morlet_cycles)
%
% Inputs:
%   data               - [channels x timepoints] continuous MEG time series.
%                         The last 15 rows are assumed to be non-data channels
%                         (e.g., trigger/status channels). Row (end-14) is
%                         used as the event trigger channel.
%   fs                  - Sampling rate in Hz (e.g., 1000).
%   freq_range          - (Optional) Frequency range for TF decomposition in Hz.
%                         Default: 13:30 (beta band).
%   epoch_window        - (Optional) [pre post] epoch window in seconds
%                         relative to event onset. Default: [1.479 1.48].
%   min_event_interval  - (Optional) Minimum interval between events in
%                         samples to avoid duplicates. Default: 3000.
%   morlet_cycles       - (Optional) Number of cycles for Morlet wavelet.
%                         Default: 3.
%
% Outputs:
%   tf_avg              - [channels x timepoints x frequencies] trial-averaged
%                         time-frequency map. Trimmed to remove 80-sample
%                         edges on each side.
%   event_onsets        - [1 x N] vector of detected event onset sample indices.
%
% Dependencies:
%   morlet_transform.m  - Built-in function from Brainstorm toolbox
%                         (Tadel et al., 2011; https://neuroimage.usc.edu/brainstorm/).
%
% Notes:
%   - Events are detected as rising edges of the most frequently occurring
%     amplitude level among the top 3 unique values in the trigger channel.
%   - Subjects with fewer than 100 detected events are flagged (events are
%     cleared), which may require manual inspection.
%   - The trigger channel is taken as abs(data(end-14, :)).
%
% Example:
%   load('subject01_meg_task.mat');  % loads variable F [channels x time]
%   [tf_avg, events] = extract_beta_tf(F, 1000);
%
% See also: morlet_transform, detect_task_events

    % --- Default parameters ---
    if nargin < 3 || isempty(freq_range)
        freq_range = 13:30;
    end
    if nargin < 4 || isempty(epoch_window)
        epoch_window = [1.479, 1.48];  % seconds before and after event
    end
    if nargin < 5 || isempty(min_event_interval)
        min_event_interval = 3000;  % samples
    end
    if nargin < 6 || isempty(morlet_cycles)
        morlet_cycles = 3;
    end

    % --- Separate data channels from trigger channels ---
    n_total_channels = size(data, 1);
    n_data_channels = n_total_channels - 15;
    data_channels = data(1:n_data_channels, :);

    % --- Detect events from trigger channel ---
    event_onsets = detect_task_events(data, min_event_interval);

    if isempty(event_onsets)
        warning('No valid events detected. Returning empty output.');
        tf_avg = [];
        return;
    end

    % --- Compute epoch time vector ---
    pre_samples = round(epoch_window(1) * fs);
    post_samples = round(epoch_window(2) * fs);
    t = -epoch_window(1):(1/fs):epoch_window(2);

    % --- Accumulate TF across trials ---
    epoch_length = pre_samples + post_samples + 1;
    tf_sum = zeros(n_data_channels, epoch_length);
    n_valid_trials = 0;

    for k = 1:(length(event_onsets) - 1)
        onset = event_onsets(k);

        % Check epoch boundaries
        if (onset - pre_samples) < 1 || (onset + post_samples) > size(data_channels, 2)
            continue;
        end

        % Extract epoch
        epoch = data_channels(:, (onset - pre_samples):(onset + post_samples));

        % Morlet wavelet time-frequency decomposition
        tf = morlet_transform(epoch, t, freq_range, 1, morlet_cycles, 'y');

        % Average across frequencies to get [channels x time]
        tf = mean(tf, 3);

        tf_sum = tf_sum + tf;
        n_valid_trials = n_valid_trials + 1;
    end

    % --- Average across trials ---
    if n_valid_trials > 0
        tf_avg = tf_sum / n_valid_trials;
    else
        warning('No valid trials after boundary check.');
        tf_avg = [];
        return;
    end

    % --- Trim edges (remove 80 samples from each side) ---
    edge_trim = 80;
    tf_avg = tf_avg(:, (edge_trim + 1):(end - edge_trim));

    fprintf('Processed %d valid trials out of %d detected events.\n', ...
        n_valid_trials, length(event_onsets));
end
