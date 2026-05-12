# Beta Burst Propagation Analysis

## Overview

MATLAB code for analyzing beta-band (13–30 Hz) burst propagation patterns in MEG data during task and resting states. The pipeline detects beta bursts from source-reconstructed time-frequency data, extracts burst onset timing, and characterizes propagation patterns using PCA and optical flow analysis.

This code accompanies the following publication:

> [Your Name et al. (Year). Paper title. *Journal Name*. DOI: xxx]

## Requirements

### Software

- **MATLAB** (tested on R2021b+)
- **Parallel Computing Toolbox** (for `parfor` in burst detection)
- **Image Processing Toolbox** (for `bwlabel`, `regionprops`)
- **Statistics and Machine Learning Toolbox** (for `pca`, `prctile`)

### External Toolboxes

- **[Brainstorm](https://neuroimage.usc.edu/brainstorm/)** (Tadel et al., 2011)
  - Used for MEG preprocessing, source reconstruction, and Morlet wavelet time-frequency decomposition (`morlet_transform`).
  - Not included in this repository. Install separately and add to MATLAB path.

- **[SPRiNT](https://github.com/lucwilson/SPRiNT)** (Wilson et al., 2022)
  - Used for time-resolved spectral parameterization of aperiodic (1/f) components.
  - Not included in this repository. Install separately and add to MATLAB path.
  - **Modified parameters** (set programmatically in our code, no modification of SPRiNT source required):

    | Parameter | Value used | SPRiNT default | Description |
    |-----------|-----------|----------------|-------------|
    | `sfreq` | 1000 | 500 | Input sampling rate (Hz) |
    | `WinLength` | 1 | 1 | STFT window length (s) |
    | `WinOverlap` | 0 | 50 | Window overlap (%) |
    | `WinAverage` | 1 | 5 | Windows averaged per time point |

## Data

This study used publicly available data from the [Cambridge Centre for Ageing and Neuroscience (Cam-CAN)](https://cam-can.mrc-cbu.cam.ac.uk/dataset/) repository:

> Taylor, J.R., Williams, N., Cusack, R. et al. (2017). The Cambridge Centre for Ageing and Neuroscience (Cam-CAN) data repository. *NeuroImage*, 144, 262–269.

## Functions

### Task Pipeline

| Function | Description |
|----------|-------------|
| `detect_task_events.m` | Detects trial onset times from MEG trigger channel via automatic threshold selection and rising-edge detection. |
| `extract_beta_tf.m` | Extracts trial-averaged beta-band time-frequency maps and aperiodic components (via SPRiNT) for task data. Calls `detect_task_events` internally. |
| `detect_beta_bursts.m` | Subtracts aperiodic component from TF map, applies 75th percentile threshold, and identifies sustained bursts (≥100 ms). Pre-event and post-event halves processed independently. |
| `detect_burst_onset.m` | Finds burst onset closest to event onset: last onset in the pre-event window, first onset in the post-event window. |
| `burst_onset_pca.m` | Sorts burst onsets within each subject, then performs group-level PCA on the sorted onset distributions. |
| `batch_beta_burst_pipeline.m` | Master batch script that runs all task pipeline steps (Steps 1–4) sequentially across all subjects. |

### Resting-State Pipeline

| Function | Description |
|----------|-------------|
| `resting_beta_burst_analysis.m` | All-in-one resting-state analysis: segments continuous data into 4 s non-overlapping epochs, computes TF and aperiodic per epoch, detects bursts, identifies the reference region (most bursts), and computes relative onset times. |

### Optical Flow Analysis

| Function | Description |
|----------|-------------|
| `compute_flow_direction.m` | Computes mean optical flow direction per channel within a 100 ms window following each channel's burst onset, projected to [0, π] using axial circular statistics. |

## Usage

### Task Pipeline (single call)

```matlab
% Run full task pipeline across all subjects
results = batch_beta_burst_pipeline('/path/to/subjects/', '/path/to/output/');

% With custom parameters
params.percentile_thresh = 80;
params.n_components = 5;
results = batch_beta_burst_pipeline('/path/to/subjects/', '/path/to/output/', params);
```

### Task Pipeline (step by step)

```matlab
% Step 1: Extract TF and aperiodic
load('subject01/meg_task/data_block001.mat');  % loads F
[tf_avg, ap_avg, events] = extract_beta_tf(F, 1000);

% Step 2: Detect bursts
beta_burst = detect_beta_bursts(tf_avg, ap_avg);

% Step 3: Detect burst onsets
[onset_pre, onset_post] = detect_burst_onset(beta_burst);

% Step 4: Group-level PCA (after collecting all subjects)
[sorted_pre, sorted_post, pca_results] = burst_onset_pca(ann1, ann2);
```

### Resting-State Pipeline

```matlab
% Load resting-state data and remove trigger channels
load('subject01_resting.mat');
data_channels = F(1:end-15, :);

% Run full resting-state analysis
[relative_onsets, ref_region, burst_counts] = resting_beta_burst_analysis(data_channels, 1000);
```

### Optical Flow Direction

```matlab
% Compute mean flow direction per channel
load('subject01_optical_flow.mat');  % loads Source_Data [200 x 2800]
means = compute_flow_direction(Source_Data, onset_pre, 100);
```

## Analysis Parameters

Key parameters used in this study (all configurable via function arguments):

| Parameter | Value | Description |
|-----------|-------|-------------|
| Sampling rate | 1000 Hz | MEG data sampling rate |
| Frequency range | 13–30 Hz | Beta band for TF decomposition |
| Morlet cycles | 3 | Wavelet cycles for TF decomposition |
| Epoch window (task) | [−1.479, 1.48] s | Time window around event onset |
| Burst threshold | 75th percentile | Applied per channel per half |
| Min burst duration | 100 ms | Minimum consecutive samples above threshold |
| Resting epoch duration | 4 s | Non-overlapping segments |
| Max onset difference (resting) | 500 ms | Threshold for excluding spurious burst pairs |
| Optical flow window | 100 ms | Post-onset window for direction analysis |

## Notes

- Event detection is fully automatic for most subjects. For a subset of participants with inconsistent trigger coding, event markers were manually inspected and corrected (see Methods section of the paper for details).
- MEG preprocessing and source reconstruction were performed using the Brainstorm GUI. Detailed processing steps and parameters are described in the Materials and Methods section of the paper.

## License

[MIT License](LICENSE)
