% plot thresholds as function of week
if ~exist('DATA_ROOT','var'), run(fullfile(pwd, 'matlab', 'config.m')); end
load(FIG2_DATA);

weeks = [];
thresholds = [];
for animal_index = 1:6
    dayThresholds = daysThresholdsAll{animal_index}(:);
    daysTrained = ones(3, 1) * daysTrainedAll{animal_index};
    daysTrained = daysTrained(:);
    
    % Filter out zeros and keep only first 4 weeks
    weeks_of_training = ceil(daysTrained / 7);
    valid_idx = dayThresholds ~= 0 & weeks_of_training < 5;
    thresholds = [thresholds; dayThresholds(valid_idx)];
    weeks = [weeks; weeks_of_training(valid_idx)];
end

figure;
scatter(weeks, thresholds, 50, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('Week of Training');
ylabel('Threshold');
title('Thresholds vs Training Week');
xlim([-0.5, 4.5])
xticks(0:4)

%% Statistical analysis

early_thresh = thresholds(weeks < 2);
late_thresh = thresholds(weeks >= 2);
early_thresh_num = length(early_thresh); 
late_thresh_num = length(late_thresh);

[p_rank, h_rank] = ranksum(early_thresh, late_thresh);
if h_rank
    fprintf('Ranksum: p=%.4f (significant)\n', p_rank);
else
    fprintf('Ranksum: p=%.4f (not significant)\n', p_rank);
end

early_median = median(early_thresh);
late_median = median(late_thresh);
percent_change_median = ((late_median - early_median) / early_median) * 100;

fprintf('Percent change (median): %.1f%%\n', percent_change_median);

statsCSV = fullfile(fileparts(mfilename('fullpath')), 'stats.csv');
report_stat(statsCSV, 'Fig2A', 'threshold early(0-1) vs late(2-4)', early_thresh, late_thresh);