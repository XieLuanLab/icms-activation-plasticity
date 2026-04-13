%% Data paths
% Set REPO_ROOT to the root of the icms-plasticity-code repository.
REPO_ROOT = fileparts(fileparts(mfilename('fullpath')));
if isempty(REPO_ROOT)
    % Fallback if mfilename fails (e.g., running from editor temp dir)
    REPO_ROOT = pwd;
end

DATA_ROOT = fullfile(REPO_ROOT, 'data', 'matlab');

% Fig 2 / SFig2: volumetric imaging source data
FIG2_DATA = fullfile(DATA_ROOT, 'figure_source_data.mat');

% Fig 3: planar imaging results
FIG3_DATA = fullfile(DATA_ROOT, 'fig3_final_results.mat');
FIG3_RASTER = fullfile(DATA_ROOT, 'fig3b_raster.mat');
FIG3_PLOT = fullfile(DATA_ROOT, 'fig3def_plot_data.mat');

% SFig2: control vs behavioral source data
SFIG2D_CONTROL = fullfile(DATA_ROOT, 'control_density_source.mat');
SFIG2D_BEHAVIOR = fullfile(DATA_ROOT, 'behavior_density_source.mat');
SFIG2E_CONTROL = fullfile(DATA_ROOT, 'control_count_source.mat');
SFIG2E_BEHAVIOR = fullfile(DATA_ROOT, 'behavior_count_source.mat');

% Only needed to reprocess from scratch. Not required for figure generation.
RAW_DATA_ROOT = '';

% Add utility scripts to path
addpath(fullfile(REPO_ROOT, 'matlab'));
