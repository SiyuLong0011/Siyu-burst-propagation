function batch_extract_beta_tf(input_dir, output_dir, fs, freq_range, epoch_window, min_event_interval, morlet_cycles)
% BATCH_EXTRACT_BETA_TF Batch process beta-band TF extraction across subjects
%
% Iterates over subject folders in the input directory, loads MEG task data,
% extracts trial-averaged beta-band time-frequency representations, and
% saves the results to the output directory.
%
% Usage:
%   batch_extract_beta_tf(input_dir, output_dir)
%   batch_extract_beta_tf(input_dir, output_dir, fs, freq_range, epoch_window, min_event_interval, morlet_cycles)
%
% Inputs:
%   input_dir           - Path to the parent directory containing subject
%                         folders. Each subject folder should contain
%                         'meg_task/data_block001.mat' with variable F.
%   output_dir          - Path to save output .mat files. Each file is named
%                         after the subject folder and contains the variable
%                         'res' (trial-averaged TF map).
%   fs                  - (Optional) Sampling rate in Hz. Default: 1000.
%   freq_range          - (Optional) Frequency range in Hz. Default: 13:30.
%   epoch_window        - (Optional) [pre post] in seconds. Default: [1.479 1.48].
%   min_event_interval  - (Optional) Min interval between events in samples.
%                         Default: 3000.
%   morlet_cycles       - (Optional) Morlet wavelet cycles. Default: 3.
%
% Outputs:
%   Saves one .mat file per subject in output_dir, containing:
%     res  - [channels x timepoints] trial-averaged TF representation.
%
% Dependencies:
%   extract_beta_tf.m, detect_task_events.m, morlet_transform.m
%
% Example:
%   batch_extract_beta_tf('/data/MEG_study/subjects/', '/data/MEG_study/results/');
%
% See also: extract_beta_tf, detect_task_events

    % --- Default parameters ---
    if nargin < 3 || isempty(fs),                  fs = 1000;              end
    if nargin < 4 || isempty(freq_range),           freq_range = 13:30;     end
    if nargin < 5 || isempty(epoch_window),         epoch_window = [1.479, 1.48]; end
    if nargin < 6 || isempty(min_event_interval),   min_event_interval = 3000;    end
    if nargin < 7 || isempty(morlet_cycles),         morlet_cycles = 3;     end

    % --- Create output directory if it does not exist ---
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % --- List subject folders ---
    file_list = dir(input_dir);
    % Remove '.', '..', and hidden files
    file_list = file_list(~ismember({file_list.name}, {'.', '..'}));
    file_list = file_list([file_list.isdir]);

    fprintf('Found %d subject folders in %s\n', length(file_list), input_dir);

    % --- Process each subject ---
    for i = 1:length(file_list)
        subject_name = file_list(i).name;
        data_path = fullfile(input_dir, subject_name, 'meg_task', 'data_block001.mat');

        fprintf('\n--- Processing subject %d/%d: %s ---\n', i, length(file_list), subject_name);

        % Check if data file exists
        if ~isfile(data_path)
            warning('Data file not found for %s. Skipping.', subject_name);
            continue;
        end

        % Load data
        loaded = load(data_path, 'F');
        if ~isfield(loaded, 'F')
            warning('Variable F not found in %s. Skipping.', data_path);
            continue;
        end

        % Extract beta TF
        [res, event_onsets] = extract_beta_tf(loaded.F, fs, freq_range, ...
            epoch_window, min_event_interval, morlet_cycles);

        if isempty(res)
            warning('No valid result for %s. Skipping save.', subject_name);
            continue;
        end

        % Save result
        output_path = fullfile(output_dir, [subject_name '.mat']);
        save(output_path, 'res');
        fprintf('Saved: %s\n', output_path);

        % Clear per-subject variables
        clear res event_onsets loaded;
    end

    fprintf('\n=== Batch processing complete. ===\n');
end
