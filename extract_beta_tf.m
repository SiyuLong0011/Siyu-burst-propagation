function [tf_avg, ap_avg, event_onsets] = extract_beta_tf(data, fs, freq_range, epoch_window, min_event_interval, morlet_cycles, sprint_opt)
% EXTRACT_BETA_TF Extract trial-averaged beta-band TF and aperiodic components
%
% Computes a Morlet wavelet-based time-frequency decomposition in the beta
% band for task-related MEG data, and simultaneously extracts the aperiodic
% (1/f) spectral component using SPRiNT. Events (trial onsets) are
% automatically detected from a trigger channel, and both TF maps and
% aperiodic models are averaged across trials.
%
% Usage:
%   [tf_avg, ap_avg, event_onsets] = extract_beta_tf(data, fs)
%   [tf_avg, ap_avg, event_onsets] = extract_beta_tf(data, fs, freq_range, ...
%       epoch_window, min_event_interval, morlet_cycles, sprint_opt)
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
%   sprint_opt          - (Optional) Struct of SPRiNT options. Fields that
%                         differ from SPRiNT defaults are set automatically:
%                           .sfreq      = 1000  (sampling rate, Hz)
%                           .WinLength  = 1     (STFT window length, s)
%                           .WinOverlap = 0     (window overlap, %)
%                           .WinAverage = 1     (windows averaged per time point)
%                         Any field provided in this struct will override
%                         the above defaults.
%
% Outputs:
%   tf_avg              - [channels x timepoints] trial-averaged time-frequency
%                         map (averaged across frequencies). Trimmed to remove
%                         80-sample edges on each side.
%   ap_avg              - [channels x time_windows] trial-averaged aperiodic
%                         model from SPRiNT, averaged across the frequency
%                         dimension (3rd dim of SPRiNT output).
%   event_onsets        - [1 x N] vector of detected event onset sample indices.
%
% Dependencies:
%   morlet_transform.m  - Built-in function from Brainstorm toolbox
%                         (Tadel et al., 2011; https://neuroimage.usc.edu/brainstorm/).
%   SPRiNT              - Spectral Parameterization Resolved in Time toolbox
%                         (Wilson et al., 2022; https://github.com/lucwilson/SPRiNT).
%
% Notes:
%   - Events are detected as rising edges of the most frequently occurring
%     amplitude level among the top 3 unique values in the trigger channel.
%   - Subjects with fewer than 100 detected events are flagged (events are
%     cleared), which may require manual inspection.
%   - The trigger channel is taken as abs(data(end-14, :)).
%   - For SPRiNT, two time segments are extracted from each epoch:
%     samples 181:1280 (pre-event) and 1681:2780 (post-event), concatenated
%     as input. The aperiodic model is averaged across the 3rd dimension.
%   - SPRiNT parameters were modified from the toolbox defaults:
%       sfreq      = 1000 (default: 500)
%       WinLength  = 1    (default: 1)
%       WinOverlap = 0    (default: 50)
%       WinAverage = 1    (default: 5)
%     See https://github.com/lucwilson/SPRiNT for original defaults.
%
% Example:
%   load('subject01_meg_task.mat');  % loads variable F [channels x time]
%   [tf_avg, ap_avg, events] = extract_beta_tf(F, 1000);
%
% See also: morlet_transform, detect_task_events, SPRiNT

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
    if nargin < 7 || isempty(sprint_opt)
        sprint_opt = struct();
    end
    % SPRiNT parameters (modified from SPRiNT defaults, see Notes)
    if ~isfield(sprint_opt, 'sfreq'),       sprint_opt.sfreq = 1000;       end
    if ~isfield(sprint_opt, 'WinLength'),   sprint_opt.WinLength = 1;      end
    if ~isfield(sprint_opt, 'WinOverlap'),  sprint_opt.WinOverlap = 0;     end
    if ~isfield(sprint_opt, 'WinAverage'),  sprint_opt.WinAverage = 1;     end

    % --- Separate data channels from trigger channels ---
    n_total_channels = size(data, 1);
    n_data_channels = n_total_channels - 15;
    data_channels = data(1:n_data_channels, :);

    % --- Detect events from trigger channel ---
    event_onsets = detect_task_events(data, min_event_interval);

    if isempty(event_onsets)
        warning('No valid events detected. Returning empty output.');
        tf_avg = [];
        ap_avg = [];
        return;
    end

    % --- Compute epoch time vector ---
    pre_samples = round(epoch_window(1) * fs);
    post_samples = round(epoch_window(2) * fs);
    t = -epoch_window(1):(1/fs):epoch_window(2);

    % --- Accumulate TF and aperiodic components across trials ---
    epoch_length = pre_samples + post_samples + 1;
    tf_sum = zeros(n_data_channels, epoch_length);
    ap_sum = [];  % initialized on first valid trial (size depends on SPRiNT output)
    n_valid_trials = 0;

    for k = 1:(length(event_onsets) - 1)
        onset = event_onsets(k);

        % Check epoch boundaries
        if (onset - pre_samples) < 1 || (onset + post_samples) > size(data_channels, 2)
            continue;
        end

        % Extract epoch
        epoch = data_channels(:, (onset - pre_samples):(onset + post_samples));

        % --- Morlet wavelet time-frequency decomposition ---
        tf = morlet_transform(epoch, t, freq_range, 1, morlet_cycles, 'y');
        tf = mean(tf, 3);  % average across frequencies
        tf_sum = tf_sum + tf;

        % --- SPRiNT aperiodic component ---
        % Extract two time segments: pre-event (181:1280) and post-event (1681:2780)
        xs = [epoch(:, 181:1280), epoch(:, 1681:2780)];
        opt = sprint_opt;
        SPRiNT;  % runs SPRiNT using variables 'xs' and 'opt' in workspace
        ap = mean(s_data.SPRiNT.aperiodic_models, 3);

        if isempty(ap_sum)
            ap_sum = zeros(size(ap));
        end
        ap_sum = ap_sum + ap;

        n_valid_trials = n_valid_trials + 1;
        clear ap s_data;
    end

    % --- Average across trials ---
    if n_valid_trials > 0
        tf_avg = tf_sum / n_valid_trials;
        ap_avg = ap_sum / n_valid_trials;
    else
        warning('No valid trials after boundary check.');
        tf_avg = [];
        ap_avg = [];
        return;
    end

    % --- Trim edges (remove 80 samples from each side) ---
    edge_trim = 80;
    tf_avg = tf_avg(:, (edge_trim + 1):(end - edge_trim));

    fprintf('Processed %d valid trials out of %d detected events.\n', ...
        n_valid_trials, length(event_onsets));
end
