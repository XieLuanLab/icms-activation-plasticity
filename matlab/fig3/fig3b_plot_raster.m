%% Fig 3B - Plot activation raster from saved data
SHUFFLE = true;  % true = shuffle within groups

if ~exist('DATA_ROOT','var'), run(fullfile(pwd, 'matlab', 'config.m')); end
load(fullfile(DATA_ROOT, 'fig3b_raster.mat'));

nROI = size(binROIRaster, 1);
fSums = sum(binROIRaster, 2);

if SHUFFLE
    % Sort by first appearance, shuffle within each group
    startInds = zeros(nROI, 1);
    for i = 1:nROI
        f = find(binROIRaster(i,:), 1, 'first');
        if isempty(f), startInds(i) = 6; else, startInds(i) = f; end
    end
    rng(42);
    order = [];
    for s = unique(startInds)'
        grpIdx = find(startInds == s);
        order = [order; grpIdx(randperm(numel(grpIdx)))];
    end
    sorted = binROIRaster(order, :);
    fSumsSorted = fSums(order);
    colorData = zeros(size(sorted));
    for i = 1:size(colorData,1)
        for w = 1:5
            if sorted(i,w) == 1
                colorData(i,w) = fSumsSorted(i);
            end
        end
    end
    dataFlipped = flipud(colorData);
else
    % Use Roy's pre-sorted heatRaster directly
    dataFlipped = flipud(heatRaster(1:nROI, :));
end

blueMap = [1 1 1;
           0.60 0.78 0.95;
           0.40 0.60 0.85;
           0.20 0.42 0.72;
           0.08 0.25 0.55;
           0.02 0.10 0.38];

figure('Position', [100 100 400 600]);
imagesc(dataFlipped);
colormap(blueMap);
caxis([-0.5 5.5]);
cb = colorbar;
cb.Ticks = [1 2 3 4 5];
cb.TickLabels = {'1','2','3','4','5'};
cb.Limits = [0.5 5.5];
cb.Label.String = 'Weeks activated';
xlabel('Week'); ylabel('Tracked neuron ID');
xticks(1:5); xticklabels({'0','1','2','3','4'});
nRows = size(dataFlipped, 1);
ypos = sort(unique([1, nRows:-200:1]));
yticks(ypos); yticklabels(string(nRows - ypos));
title(sprintf('n = %d', nROI));
set(gca, 'FontName', 'Arial', 'FontSize', 8);
box on;

% Print summary
wkActive = sum(binROIRaster, 2);
fprintf('Total: %d\n', nROI);
for w = 5:-1:1
    fprintf('  %d/5 weeks: %d\n', w, sum(wkActive == w));
end
