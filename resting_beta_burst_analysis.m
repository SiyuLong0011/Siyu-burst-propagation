function [relative_onsets, ref_region, burst_counts] = resting_beta_burst_analysis(data, fs, params)
% RESTING_BETA_BURST_ANALYSIS Full resting-state beta burst onset pipeline
%
% Performs the complete resting-state beta burst analysis in one call:
% segmentation, TF decomposition, aperiodic removal, burst detection,
% and relative onset computation.
%
% Usage:
%   [relative_onsets, ref_region, burst_counts] = resting_beta_burst_analysis(data, fs)
%   [relative_onsets, ref_region, burst_counts] = resting_beta_burst_analysis(data, fs, params)
%
% Inputs:
%   data            - [channels x timepoints] continuous resting-state MEG
%                     data (data channels only, no trigger channels).
%   fs              - Sampling rate in Hz (e.g., 1000).
%   params          - (Optional) Struct with analysis parameters. Any missing
%                     field uses its default value. Fields:
%                       .epoch_duration     - Epoch length in seconds. Default: 4.
%                       .freq_range         - Frequency range (Hz). Default: 13:30.
%                       .morlet_cycles      - Morlet wavelet cycles. Default: 3.
%                       .sprint_opt         - SPRiNT options struct. See extract_resting_tf.
%                       .percentile_thresh  - Burst threshold percentile. Default: 75.
%                       .min_burst_duration - Min burst duration (samples). Default: 100.
%                       .max_onset_diff     - Max onset difference for pairing (ms). Default: 500.
%
% Outputs:
%   relative_onsets - [channels x 1] cell array. Each cell contains a
%                     vector of relative onset times (in ms) for that
%                     channel, computed as (channel_onset - nearest_ref_onset).
%                     Positive = after reference, negative = before.
%                     The reference region cell contains all zeros.
%   ref_region      - Scalar index of the reference region (channel with
%                     highest total burst count across all epochs).
%   burst_counts    - [channels x 1] total number of bursts detected per
%                     channel across all epochs.
%
% Algorithm:
%   1. Segment continuous data into non-overlapping epochs.
%   2. For each epoch:
%      a. Morlet wavelet TF decomposition (beta band).
%      b. SPRiNT aperiodic parameterization.
%      c. Subtract aperiodic component, threshold at 75th percentile,
%         retain bursts >= min_burst_duration.
%   3. Across all epochs:
%      a. Identify reference region (most total bursts).
%      b. For each burst in non-reference channels, find nearest burst
%         in reference channel (same epoch) and compute time difference.
%      c. Exclude pairs with |difference| > max_onset_diff.
%
% Dependencies:
%   morlet_transform.m (Brainstorm), SPRiNT, Image Processing Toolbox
%   (bwlabel, regionprops)
%
% Example:
%   load('subject01_resting.mat');
%   data_channels = F(1:end-15, :);
%   [relative_onsets, ref_region, counts] = resting_beta_burst_analysis(data_channels, 1000);
%
%   % Plot mean relative onset per channel
%   mean_onset = cellfun(@mean, relative_onsets);
%   figure; bar(mean_onset); xlabel('Channel'); ylabel('Mean relative onset (ms)');
%
% See also: extract_beta_tf, detect_beta_bursts

    % --- Default parameters ---
    if nargin < 3 || isempty(params)
        params = struct();
    end
    if ~isfield(params, 'epoch_duration'),     params.epoch_duration = 4;         end
    if ~isfield(params, 'freq_range'),         params.freq_range = 13:30;         end
    if ~isfield(params, 'morlet_cycles'),      params.morlet_cycles = 3;          end
    if ~isfield(params, 'percentile_thresh'),  params.percentile_thresh = 75;     end
    if ~isfield(params, 'min_burst_duration'), params.min_burst_duration = 100;   end
    if ~isfield(params, 'max_onset_diff'),     params.max_onset_diff = 500;       end
    if ~isfield(params, 'sprint_opt'),         params.sprint_opt = struct();      end

    % SPRiNT defaults
    if ~isfield(params.sprint_opt, 'sfreq'),       params.sprint_opt.sfreq = 1000;       end
    if ~isfield(params.sprint_opt, 'WinLength'),   params.sprint_opt.WinLength = 1;      end
    if ~isfield(params.sprint_opt, 'WinOverlap'),  params.sprint_opt.WinOverlap = 0;     end
    if ~isfield(params.sprint_opt, 'WinAverage'),  params.sprint_opt.WinAverage = 1;     end

    % --- Step 1: Segment into epochs ---
    epoch_samples = round(params.epoch_duration * fs);
    [n_channels, total_samples] = size(data);
    n_epochs = floor(total_samples / epoch_samples);

    if n_epochs == 0
        error('Data is shorter than one epoch (%d s).', params.epoch_duration);
    end

    fprintf('Segmenting into %d epochs of %d s...\n', n_epochs, params.epoch_duration);

    % --- Step 2-3: Process each epoch (TF + aperiodic + burst detection + onset extraction) ---
    all_onsets = cell(n_channels, n_epochs);
    burst_counts = zeros(n_channels, 1);
    t = (0:(epoch_samples - 1)) / fs;

    for e = 1:n_epochs
        if mod(e, 50) == 0 || e == 1
            fprintf('  Processing epoch %d/%d...\n', e, n_epochs);
        end

        % Extract epoch
        start_idx = (e - 1) * epoch_samples + 1;
        end_idx = e * epoch_samples;
        epoch = data(:, start_idx:end_idx);

        % --- TF decomposition ---
        tf = morlet_transform(epoch, t, params.freq_range, 1, params.morlet_cycles, 'y');
        tf_epoch = mean(tf, 3);  % [channels x time]

        % --- SPRiNT aperiodic ---
        xs = epoch;
        opt = params.sprint_opt;
        SPRiNT;
        ap_epoch = mean(mean(s_data.SPRiNT.aperiodic_models, 3), 2);  % [channels x 1]
        clear s_data;

        % --- Subtract aperiodic and detect bursts ---
        tf_corrected = tf_epoch - repmat(ap_epoch, 1, epoch_samples);

        for ch = 1:n_channels
            x = tf_corrected(ch, :);
            thresh = prctile(x, params.percentile_thresh);
            bin_x = x >= thresh;

            % Label contiguous regions and filter by duration
            [L, num] = bwlabel(bin_x);
            burst_signal = zeros(1, epoch_samples);
            if num > 0
                stats = regionprops(L, 'Area', 'PixelIdxList');
                for j = 1:num
                    if stats(j).Area >= params.min_burst_duration
                        burst_signal(stats(j).PixelIdxList) = 1;
                    end
                end
            end

            % Detect rising edges as burst onsets
            onsets = find(diff([0, burst_signal]) == 1);
            all_onsets{ch, e} = onsets;
            burst_counts(ch) = burst_counts(ch) + length(onsets);
        end
    end

    % --- Step 4: Identify reference region ---
    [~, ref_region] = max(burst_counts);
    fprintf('Reference region: channel %d (%d total bursts).\n', ...
        ref_region, burst_counts(ref_region));

    % --- Step 5: Compute relative onset times ---
    relative_onsets = cell(n_channels, 1);

    for ch = 1:n_channels
        if ch == ref_region
            relative_onsets{ch} = zeros(1, burst_counts(ch));
            continue;
        end

        onset_diffs = [];
        for e = 1:n_epochs
            ch_onsets = all_onsets{ch, e};
            ref_onsets = all_onsets{ref_region, e};

            if isempty(ch_onsets) || isempty(ref_onsets)
                continue;
            end

            for b = 1:length(ch_onsets)
                diffs = abs(ch_onsets(b) - ref_onsets);
                [min_diff, nearest_idx] = min(diffs);
                time_diff = ch_onsets(b) - ref_onsets(nearest_idx);

                if abs(time_diff) <= params.max_onset_diff
                    onset_diffs = [onset_diffs, time_diff]; %#ok<AGROW>
                end
            end
        end
        relative_onsets{ch} = onset_diffs;
    end

    % --- Summary ---
    n_valid = sum(cellfun(@length, relative_onsets)) - burst_counts(ref_region);
    n_non_ref_bursts = sum(burst_counts) - burst_counts(ref_region);
    fprintf('\nResting-state analysis complete.\n');
    fprintf('  %d epochs, %d channels.\n', n_epochs, n_channels);
    fprintf('  Total bursts: %d\n', sum(burst_counts));
    fprintf('  Valid onset pairs: %d / %d non-reference bursts.\n', n_valid, n_non_ref_bursts);
    fprintf('  Excluded (>%d ms): %d\n', params.max_onset_diff, n_non_ref_bursts - n_valid);
end
