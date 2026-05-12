function event_onsets = detect_task_events(data, min_interval, min_event_count)
% DETECT_TASK_EVENTS Detect task event onsets from MEG trigger channel
%
% Identifies trial onset times by detecting rising edges in the trigger
% channel. The trigger channel is assumed to be row (end-14) of the data
% matrix. The function automatically determines the threshold by selecting
% the most frequently occurring value among the top 3 unique amplitude
% levels in the trigger channel.
%
% Usage:
%   event_onsets = detect_task_events(data)
%   event_onsets = detect_task_events(data, min_interval, min_event_count)
%
% Inputs:
%   data            - [channels x timepoints] raw MEG data matrix including
%                     trigger/status channels. Row (end-14) is used as the
%                     trigger channel.
%   min_interval    - (Optional) Minimum interval between consecutive events
%                     in samples, used to remove spurious detections.
%                     Default: 3000.
%   min_event_count - (Optional) Minimum number of events required. If fewer
%                     events are detected, the function returns an empty
%                     array and issues a warning. Default: 100.
%
% Outputs:
%   event_onsets    - [1 x N] vector of sample indices corresponding to
%                     detected event onsets. Returns empty if fewer than
%                     min_event_count events are found.
%
% Algorithm:
%   1. Extract absolute values from trigger channel (row end-14).
%   2. Identify the top 3 unique amplitude values.
%   3. Select the one with the highest occurrence count as threshold.
%   4. Detect rising edges: transitions from non-threshold to threshold.
%   5. Remove events closer than min_interval samples apart.
%
% Example:
%   load('subject01_meg_task.mat');  % loads variable F
%   events = detect_task_events(F, 3000, 100);
%
% See also: extract_beta_tf

    % --- Default parameters ---
    if nargin < 2 || isempty(min_interval)
        min_interval = 3000;
    end
    if nargin < 3 || isempty(min_event_count)
        min_event_count = 100;
    end

    % --- Extract trigger channel ---
    trigger = abs(data(end - 14, :));

    % --- Determine threshold automatically ---
    % Find the top unique values by magnitude (up to 3)
    unique_vals = unique(trigger);
    n_top = min(3, length(unique_vals));

    if n_top == 0
        warning('Trigger channel is empty. Returning empty.');
        event_onsets = [];
        return;
    end

    top_vals = maxk(unique_vals, n_top);

    % Count occurrences of each top value
    counts = zeros(1, n_top);
    for j = 1:n_top
        counts(j) = sum(trigger == top_vals(j));
    end

    % Use the most frequent among the top values as the event threshold
    [~, idx] = max(counts);
    threshold = top_vals(idx);

    % --- Detect rising edges ---
    % A rising edge is defined as a transition from a non-threshold value
    % to the threshold value
    raw_onsets = [];
    n = 1;
    for k = 2:length(trigger)
        if trigger(k - 1) ~= threshold && trigger(k) == threshold
            raw_onsets(n) = k; %#ok<AGROW>
            n = n + 1;
        end
    end

    % --- Check minimum event count ---
    if length(raw_onsets) < min_event_count
        warning('Only %d events detected (threshold: %d). Returning empty.', ...
            length(raw_onsets), min_event_count);
        event_onsets = [];
        return;
    end

    % --- Remove events too close together ---
    event_onsets = raw_onsets(1);
    for k = 2:length(raw_onsets)
        if raw_onsets(k) - event_onsets(end) >= min_interval
            event_onsets = [event_onsets, raw_onsets(k)]; %#ok<AGROW>
        end
    end

    fprintf('Detected %d events (threshold value: %.2f).\n', ...
        length(event_onsets), threshold);
end
