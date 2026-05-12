function results = batch_beta_burst_pipeline(input_dir, output_dir, params)
% BATCH_BETA_BURST_PIPELINE Full pipeline for beta burst analysis
%
% Runs all analysis steps sequentially:
%   Step 1: Extract trial-averaged TF and aperiodic components
%   Step 2: Detect beta bursts (binarization after removing aperiodic)
%   Step 3: Detect burst onset times closest to event
%   Step 4: Sort onsets within subjects and perform group-level PCA
%
% Usage:
%   results = batch_beta_burst_pipeline(input_dir, output_dir)
%   results = batch_beta_burst_pipeline(input_dir, output_dir, params)
%
% Inputs:
%   input_dir   - Path to parent directory containing subject folders.
%                 Each subject folder should contain
%                 'meg_task/data_block001.mat' with variable F.
%   output_dir  - Path to save all output files. Subdirectories are created
%                 automatically:
%                   output_dir/tf/        - TF and aperiodic results
%                   output_dir/bursts/    - binary burst matrices
%                   output_dir/onsets/    - burst onset results
%                   output_dir/pca/       - PCA results
%   params      - (Optional) Struct with analysis parameters. Any missing
%                 field uses its default value. Fields:
%
%                 Step 1 (extract_beta_tf):
%                   .fs                 - Sampling rate (Hz). Default: 1000.
%                   .freq_range         - Frequency range (Hz). Default: 13:30.
%                   .epoch_window       - [start end] in seconds. Default: [-1.479 1.48].
%                   .min_event_interval - Min event interval (samples). Default: 3000.
%                   .morlet_cycles      - Morlet cycles. Default: 3.
%                   .sprint_opt         - SPRiNT options struct. Default: see extract_beta_tf.
%
%                 Step 2 (detect_beta_bursts):
%                   .percentile_thresh  - Burst threshold percentile. Default: 75.
%                   .min_burst_duration - Min burst duration (samples). Default: 100.
%
%                 Step 3 (detect_burst_onset):
%                   .pre_window         - Pre-event search window. Default: [401 1400].
%                   .post_window        - Post-event search window. Default: [1401 2400].
%                   .search_range       - Backward extension range. Default: [1 400].
%
%                 Step 4 (burst_onset_pca):
%                   .n_components       - Number of PCA components. Default: 10.
%
% Outputs:
%   results     - Struct with fields:
%       .ann1           - [channels x valid_subjects] pre-event burst onsets (ms).
%       .ann2           - [channels x valid_subjects] post-event burst onsets (ms).
%       .sorted_pre     - [channels x valid_subjects] sorted pre-event onsets.
%       .sorted_post    - [channels x valid_subjects] sorted post-event onsets.
%       .pca_results    - PCA output struct (see burst_onset_pca).
%       .subject_names  - {1 x valid_subjects} cell array of processed subject names.
%       .skipped        - {1 x N} cell array of skipped subject names.
%
% Dependencies:
%   extract_beta_tf.m, detect_task_events.m, detect_beta_bursts.m,
%   detect_burst_onset.m, burst_onset_pca.m,
%   morlet_transform.m (Brainstorm), SPRiNT
%
% Example:
%   results = batch_beta_burst_pipeline('/data/subjects/', '/data/output/');
%
%   % With custom parameters
%   params.percentile_thresh = 80;
%   params.n_components = 5;
%   results = batch_beta_burst_pipeline('/data/subjects/', '/data/output/', params);
%
% See also: extract_beta_tf, detect_beta_bursts, detect_burst_onset, burst_onset_pca

    % --- Default parameters ---
    if nargin < 3 || isempty(params)
        params = struct();
    end

    % Step 1 defaults
    if ~isfield(params, 'fs'),                  params.fs = 1000;                       end
    if ~isfield(params, 'freq_range'),          params.freq_range = 13:30;              end
    if ~isfield(params, 'epoch_window'),        params.epoch_window = [-1.479, 1.48];   end
    if ~isfield(params, 'min_event_interval'),  params.min_event_interval = 3000;       end
    if ~isfield(params, 'morlet_cycles'),       params.morlet_cycles = 3;               end
    if ~isfield(params, 'sprint_opt'),          params.sprint_opt = struct();            end

    % Step 2 defaults
    if ~isfield(params, 'percentile_thresh'),   params.percentile_thresh = 75;          end
    if ~isfield(params, 'min_burst_duration'),  params.min_burst_duration = 100;        end

    % Step 3 defaults
    if ~isfield(params, 'pre_window'),          params.pre_window = [401, 1400];        end
    if ~isfield(params, 'post_window'),         params.post_window = [1401, 2400];      end
    if ~isfield(params, 'search_range'),        params.search_range = [1, 400];         end

    % Step 4 defaults
    if ~isfield(params, 'n_components'),        params.n_components = 10;               end

    % --- Create output subdirectories ---
    tf_dir = fullfile(output_dir, 'tf');
    burst_dir = fullfile(output_dir, 'bursts');
    onset_dir = fullfile(output_dir, 'onsets');
    pca_dir = fullfile(output_dir, 'pca');

    dirs = {tf_dir, burst_dir, onset_dir, pca_dir};
    for d = 1:length(dirs)
        if ~exist(dirs{d}, 'dir'), mkdir(dirs{d}); end
    end

    % --- List subject folders ---
    file_list = dir(input_dir);
    file_list = file_list(~ismember({file_list.name}, {'.', '..'}));
    file_list = file_list([file_list.isdir]);
    n_subjects = length(file_list);

    fprintf('========================================\n');
    fprintf('  Beta Burst Analysis Pipeline\n');
    fprintf('  %d subjects found\n', n_subjects);
    fprintf('========================================\n');

    % --- Track valid and skipped subjects ---
    valid_names = {};
    skipped_names = {};
    ann1_list = {};
    ann2_list = {};

    % --- Process each subject ---
    for i = 1:n_subjects
        subject_name = file_list(i).name;
        data_path = fullfile(input_dir, subject_name, 'meg_task', 'data_block001.mat');

        fprintf('\n=== Subject %d/%d: %s ===\n', i, n_subjects, subject_name);

        % Check if data file exists
        if ~isfile(data_path)
            warning('Data file not found for %s. Skipping.', subject_name);
            skipped_names{end + 1} = subject_name; %#ok<AGROW>
            continue;
        end

        % --- Step 1: Extract TF and aperiodic ---
        fprintf('  Step 1: Extracting TF and aperiodic...\n');
        loaded = load(data_path, 'F');
        if ~isfield(loaded, 'F')
            warning('Variable F not found in %s. Skipping.', data_path);
            skipped_names{end + 1} = subject_name; %#ok<AGROW>
            continue;
        end

        [res, ap_avg, ~] = extract_beta_tf(loaded.F, params.fs, params.freq_range, ...
            params.epoch_window, params.min_event_interval, params.morlet_cycles, ...
            params.sprint_opt);

        if isempty(res)
            warning('No valid TF result for %s. Skipping.', subject_name);
            skipped_names{end + 1} = subject_name; %#ok<AGROW>
            continue;
        end

        save(fullfile(tf_dir, [subject_name '.mat']), 'res', 'ap_avg');
        clear loaded;

        % --- Step 2: Detect beta bursts ---
        fprintf('  Step 2: Detecting beta bursts...\n');
        beta_burst = detect_beta_bursts(res, ap_avg, ...
            params.percentile_thresh, params.min_burst_duration);

        save(fullfile(burst_dir, [subject_name '.mat']), 'beta_burst');
        clear res ap_avg;

        % --- Step 3: Detect burst onsets ---
        fprintf('  Step 3: Detecting burst onsets...\n');
        [onset_pre, onset_post] = detect_burst_onset(beta_burst, ...
            params.pre_window, params.post_window, params.search_range);

        save(fullfile(onset_dir, [subject_name '.mat']), 'onset_pre', 'onset_post');

        % Accumulate valid subjects
        valid_names{end + 1} = subject_name; %#ok<AGROW>
        ann1_list{end + 1} = onset_pre; %#ok<AGROW>
        ann2_list{end + 1} = onset_post; %#ok<AGROW>

        clear beta_burst onset_pre onset_post;
    end

    % --- Assemble onset matrices from valid subjects only ---
    n_valid = length(valid_names);
    if n_valid == 0
        error('No valid subjects processed. Cannot proceed to PCA.');
    end

    n_channels = length(ann1_list{1});
    ann1 = zeros(n_channels, n_valid);
    ann2 = zeros(n_channels, n_valid);
    for i = 1:n_valid
        ann1(:, i) = ann1_list{i};
        ann2(:, i) = ann2_list{i};
    end

    % --- Step 4: Group-level PCA ---
    fprintf('\n=== Step 4: Group-level PCA (%d valid subjects) ===\n', n_valid);
    [sorted_pre, sorted_post, pca_results] = burst_onset_pca(ann1, ann2, params.n_components);

    subject_names = valid_names;
    save(fullfile(pca_dir, 'group_pca_results.mat'), ...
        'ann1', 'ann2', 'sorted_pre', 'sorted_post', 'pca_results', 'subject_names');

    % --- Assemble output ---
    results.ann1 = ann1;
    results.ann2 = ann2;
    results.sorted_pre = sorted_pre;
    results.sorted_post = sorted_post;
    results.pca_results = pca_results;
    results.subject_names = valid_names;
    results.skipped = skipped_names;

    fprintf('\n========================================\n');
    fprintf('  Pipeline complete.\n');
    fprintf('  Processed: %d subjects\n', n_valid);
    fprintf('  Skipped:   %d subjects\n', length(skipped_names));
    fprintf('  Results saved to: %s\n', output_dir);
    fprintf('========================================\n');
end
