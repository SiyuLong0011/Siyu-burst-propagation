function [onset_pre, onset_post] = detect_burst_onset(beta_burst, pre_window, post_window, search_range)
% DETECT_BURST_ONSET Find burst onset closest to event onset in each channel
%
% For each channel, identifies the onset of the beta burst closest to the
% event onset (time = 0). Pre-event and post-event windows are processed
% independently:
%   - Pre-event:  finds the LAST burst onset (closest to time 0)
%   - Post-event: finds the FIRST burst onset (closest to time 0)
%
% If a burst in the pre-event window is already active at the start of the
% search range, the function extends the search backward to find its true
% onset.
%
% Usage:
%   [onset_pre, onset_post] = detect_burst_onset(beta_burst)
%   [onset_pre, onset_post] = detect_burst_onset(beta_burst, pre_window, post_window, search_range)
%
% Inputs:
%   beta_burst      - [channels x timepoints] binary burst matrix (output
%                     from detect_beta_bursts.m). 1 = burst active, 0 = no burst.
%                     Assumed to span -1400 to +1400 ms at 1 ms resolution
%                     (2800 timepoints total).
%   pre_window      - (Optional) [start end] sample indices for the pre-event
%                     search window. Default: [401 1400] (corresponding to
%                     -1000 to 0 ms).
%   post_window     - (Optional) [start end] sample indices for the post-event
%                     search window. Default: [1401 2400] (corresponding to
%                     0 to +1000 ms).
%   search_range    - (Optional) [start end] sample indices for backward
%                     extension when a burst is already active at the start
%                     of the pre-event window. Default: [1 400].
%
% Outputs:
%   onset_pre       - [channels x 1] onset time of the last burst in the
%                     pre-event window (closest to event), in ms relative
%                     to event onset (negative values). 0 if no burst found.
%   onset_post      - [channels x 1] onset time of the first burst in the
%                     post-event window (closest to event), in ms relative
%                     to event onset (positive values). 0 if no burst found.
%
% Notes:
%   - Time conversion assumes 1 kHz sampling rate (1 sample = 1 ms) and
%     that sample 1401 corresponds to time 0 ms.
%
% Example:
%   load('subject01_bursts.mat');  % loads beta_burst [channels x 2800]
%   [onset_pre, onset_post] = detect_burst_onset(beta_burst);
%
% See also: detect_beta_bursts

    % --- Default parameters ---
    if nargin < 2 || isempty(pre_window)
        pre_window = [401, 1400];
    end
    if nargin < 3 || isempty(post_window)
        post_window = [1401, 2400];
    end
    if nargin < 4 || isempty(search_range)
        search_range = [1, 400];
    end

    % --- Get dimensions ---
    n_channels = size(beta_burst, 1);
    onset_pre = zeros(n_channels, 1);
    onset_post = zeros(n_channels, 1);

    % Sample 1401 = time 0 ms
    total_samples = size(beta_burst, 2);
    event_sample = floor(total_samples / 2) + 1;  % 1401

    for j = 1:n_channels

        % ===== PRE-EVENT WINDOW =====
        % Find the LAST rising edge (closest to event onset)
        an_pre = beta_burst(j, pre_window(1):pre_window(2));
        pre_len = length(an_pre);

        % Search from end to beginning for the last rising edge
        last_onset_found = false;
        for l = pre_len:-1:2
            if an_pre(l) == 1 && an_pre(l - 1) == 0
                onset_sample = pre_window(1) + l - 1;  % global sample index
                onset_pre(j) = onset_sample - event_sample;  % convert to ms
                last_onset_found = true;
                break;
            end
        end

        % Check if burst is active at the very first sample of the window
        if ~last_onset_found && an_pre(1) == 1
            % This means the entire window is one continuous burst,
            % or the only burst starts before the window.
            % Extend backward to find the true onset.
            expanded = beta_burst(j, search_range(1):search_range(2));
            zero_pos = find(expanded == 0, 1, 'last');

            if isempty(zero_pos)
                % Burst extends all the way to the beginning of the data
                onset_sample = search_range(1);
            else
                onset_sample = search_range(1) + zero_pos;  % one sample after last zero
            end
            onset_pre(j) = onset_sample - event_sample;
            last_onset_found = true;
        end

        if ~last_onset_found
            onset_pre(j) = 0;  % no burst found
        end

        % ===== POST-EVENT WINDOW =====
        % Find the FIRST rising edge (closest to event onset)
        an_post = beta_burst(j, post_window(1):post_window(2));

        if an_post(1) == 1
            % Burst is already active at time 0
            onset_post(j) = 0;
        else
            burst_found = false;
            for l = 2:length(an_post)
                if an_post(l) == 1 && an_post(l - 1) == 0
                    onset_sample = post_window(1) + l - 1;  % global sample index
                    onset_post(j) = onset_sample - event_sample;  % convert to ms
                    burst_found = true;
                    break;
                end
            end
            if ~burst_found
                onset_post(j) = 0;
            end
        end
    end

    fprintf('Burst onset detection complete: %d channels.\n', n_channels);
end
