%% Fig 3I - dF/Fmax distribution early vs late
if ~exist('DATA_ROOT','var'), run(fullfile(pwd, 'matlab', 'config.m')); end
load(fullfile(DATA_ROOT, 'fig3def_plot_data.mat'));

% Pool all data points from weeks 0-1 and weeks 3-4
earlyVals = reshape(mergedDFF(:,1:2), 1, []);
earlyVals = earlyVals(isfinite(earlyVals));

lateVals = reshape(mergedDFF(:,4:5), 1, []);
lateVals = lateVals(isfinite(lateVals));

% Threshold on late data
thresh = mean(lateVals) + std(lateVals);

fprintf('Early (weeks 0-1): n=%d, mean=%.3f\n', numel(earlyVals), mean(earlyVals));
fprintf('Late (weeks 3-4): n=%d, mean=%.3f\n', numel(lateVals), mean(lateVals));
fprintf('Threshold (mean+STD): %.3f\n', thresh);
fprintf('Above threshold: %d\n', sum(lateVals > thresh));
fprintf('Below threshold: %d\n', sum(lateVals <= thresh));

% Build edges so threshold falls exactly on a bin boundary
binWidth = 0.08;
edgesBelow = fliplr(thresh:-binWidth:0.2);
edgesAbove = thresh:binWidth:3.7;
edges = unique([edgesBelow, edgesAbove]);

% Top: weeks 0-1
figure('Position', [100 200 400 400]);
subplot(2,1,1);
histogram(earlyVals, edges, 'FaceColor', [0.5 0.5 0.5]);
hold on
xline(mean(earlyVals), 'r', 'LineWidth', 1.5);
xlim([0.2 3.7]); ylim([0 30]);
xticks([0.5 1.5 2.5 3.5]); yticks([0 15 30]);
title('Weeks 0-1');
xlabel('dF/F_{max}'); ylabel('Count');
box off;

% Bottom: weeks 3-4, blue below threshold, red above
subplot(2,1,2);
histogram(lateVals(lateVals <= thresh), edges, 'FaceColor', 'b');
hold on
histogram(lateVals(lateVals > thresh), edges, 'FaceColor', 'r');
xline(mean(lateVals), 'k-', 'LineWidth', 1.5);
xline(thresh, 'k--', 'LineWidth', 1.5);
xlim([0.2 3.7]); ylim([0 30]);
xticks([0.5 1.5 2.5 3.5]); yticks([0 15 30]);
title('Weeks 3-4');
xlabel('dF/F_{max}'); ylabel('Count');
box off;
