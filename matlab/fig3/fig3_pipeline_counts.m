%% Analyze Fig 3 neuron/response counts across the full pipeline
% Requires Plotter2Results and allRasterData
% Outputs a summary of counts at each processing stage
if ~exist('DATA_ROOT','var'), run(fullfile(pwd, 'matlab', 'config.m')); end

fprintf('=== Loading data ===\n');
load(FIG3_DATA);
load(FIG3_RASTER);

%% ========== STAGE 1: Fig 3B Raster (allRasterData.mat) ==========
fprintf('\n========================================\n');
fprintf('STAGE 1: Fig 3B Raster\n');
fprintf('========================================\n');
fprintf('Source: allRasterData.mat (fActBin)\n');
fprintf('Each row = unique neuron, columns = 5 weekly bins\n\n');

weekCounts_raster = sum(fActBin, 2);
fprintf('Total tracked neurons: %d\n', size(fActBin, 1));
for w = 5:-1:1
    fprintf('  Active %d/5 weeks: %d neurons\n', w, sum(weekCounts_raster == w));
end
fprintf('  Active 4+/5 weeks: %d neurons\n', sum(weekCounts_raster >= 4));

%% ========== STAGE 2: Spatial matching via PatchWarp (allROI2Comp) ==========
fprintf('\n========================================\n');
fprintf('STAGE 2: Spatially tracked ROIs (PatchWarp)\n');
fprintf('========================================\n');
fprintf('Source: Plotter2Results (allROI2Comp)\n');
fprintf('Neurons with pixel overlap across >minSetCnt sessions\n\n');

animalNames = {'ICMS92', 'ICMS98', 'ICMS93'};
minSetCnts = [4, 3, 3]; % from roiPLotterBackgroundGeneration.m
nSessions = [7, 7, 5];
total_stage2 = 0;

for a = 1:3
    r = allROI2Comp{a};
    nROI = size(r, 1);
    nSess = size(r, 2);
    sessPresent = sum(r > 0, 2);
    total_stage2 = total_stage2 + nROI;
    fprintf('  %s: %d tracked ROIs (%d sessions, minSetCnt=%d)\n', ...
        animalNames{a}, nROI, nSess, minSetCnts(a));
    for w = nSess:-1:1
        if sum(sessPresent == w) > 0
            fprintf('    %d/%d sessions: %d ROIs\n', w, nSess, sum(sessPresent == w));
        end
    end
end
fprintf('  Total spatially tracked: %d\n', total_stage2);

%% ========== STAGE 3: Response units after activation filtering (allDFF) ==========
fprintf('\n========================================\n');
fprintf('STAGE 3: Activation-validated response units\n');
fprintf('========================================\n');
fprintf('Source: Plotter2Results (allDFF)\n');
fprintf('Each entry = ROI x channel x current that passed:\n');
fprintf('  - 2x CV amplitude threshold\n');
fprintf('  - >5 consecutive suprathreshold timepoints\n');
fprintf('  - Onset < 700 ms\n');
fprintf('  - Active in >minSetCnt sessions for that stim condition\n\n');

total_stage3 = 0;
for a = 1:3
    nUnits = size(allDFF{a}, 1);
    nWeeks = size(allDFF{a}, 2);
    total_stage3 = total_stage3 + nUnits;
    fprintf('  %s: %d response units x %d weeks\n', animalNames{a}, nUnits, nWeeks);
end
fprintf('  Total response units: %d\n', total_stage3);

% Check unique ROIs for animal 3 (only one with allTrendsRoiCh saved)
if exist('allTrendsRoiCh', 'var')
    uniqueROIs_a3 = numel(unique(allTrendsRoiCh(:,1)));
    fprintf('\n  Animal 3 detail (allTrendsRoiCh available):\n');
    fprintf('    %d response units from %d unique ROIs\n', size(allDFF{3},1), uniqueROIs_a3);
    fprintf('    Duplication rate: %.1f%%\n', (1 - uniqueROIs_a3/size(allDFF{3},1)) * 100);

    % Extrapolate
    ratio = size(allDFF{3},1) / uniqueROIs_a3;
    est_total_unique = round(total_stage3 / ratio);
    fprintf('    Estimated total unique ROIs across all animals: ~%d (using ratio %.2f)\n', est_total_unique, ratio);
end

%% ========== STAGE 4: Merged across animals (mergedDFF) ==========
fprintf('\n========================================\n');
fprintf('STAGE 4: Merged into weekly bins (mergedDFF)\n');
fprintf('========================================\n');
fprintf('Source: Plotter2Results (mergedDFF)\n');
fprintf('Sessions merged into 5 weekly bins per response unit\n\n');

fprintf('  mergedDFF size: %d x %d\n', size(mergedDFF));
fprintf('  Total response units: %d\n', size(mergedDFF, 1));

% Check data presence per week
for w = 1:5
    valid_w = isfinite(mergedDFF(:,w)) & mergedDFF(:,w) ~= 0;
    fprintf('  Week %d: %d units with valid data\n', w-1, sum(valid_w));
end

%% ========== STAGE 5: Valid early + late for D-F stats ==========
fprintf('\n========================================\n');
fprintf('STAGE 5: Valid early + late (D-F stats)\n');
fprintf('========================================\n');
fprintf('Source: mean(weeks 0-1) vs mean(weeks 2-4), both finite\n\n');

metrics = {mergedDFF, mergedTrendsOnset, mergedTrendsDur};
metricNames = {'dF/Fmax', 'Onset time', 'Duration'};

for t = 1:3
    metricTrend = metrics{t};
    earlyData = mean(metricTrend(:, 1:2), 2, 'omitnan');
    lateData = mean(metricTrend(:, 3:5), 2, 'omitnan');
    valid = isfinite(earlyData) & isfinite(lateData);
    fprintf('  %s: %d of %d have valid early+late data\n', metricNames{t}, sum(valid), size(metricTrend,1));
end

%% ========== STAGE 6: Subset vs remaining (Fig 3I-J) ==========
fprintf('\n========================================\n');
fprintf('STAGE 6: Subset vs remaining (Fig 3I-N)\n');
fprintf('========================================\n');
fprintf('Source: subsetMerged\n\n');

subset_idx = subsetMerged == 1;
remaining_idx = subsetMerged == 0;
fprintf('  High-intensity subgroup: %d\n', sum(subset_idx));
fprintf('  Remaining: %d\n', sum(remaining_idx));
fprintf('  Total: %d\n', numel(subsetMerged));

% Per-week counts for subset vs remaining
subDFF = mergedDFF(subset_idx, :);
remDFF = mergedDFF(remaining_idx, :);
fprintf('\n  Per-week valid data (subset / remaining):\n');
for w = 1:5
    sub_valid = sum(isfinite(subDFF(:,w)) & subDFF(:,w) ~= 0);
    rem_valid = sum(isfinite(remDFF(:,w)) & remDFF(:,w) ~= 0);
    fprintf('    Week %d: %d / %d\n', w-1, sub_valid, rem_valid);
end

% Early vs late for subset
earSub = mean(subDFF(:, 1:2), 2, 'omitnan'); earSub_valid = sum(isfinite(earSub));
latSub = mean(subDFF(:, 3:5), 2, 'omitnan'); latSub_valid = sum(isfinite(latSub));
fprintf('\n  Subset early vs late: n_early=%d, n_late=%d\n', earSub_valid, latSub_valid);

earRem = mean(remDFF(:, 1:2), 2, 'omitnan'); earRem_valid = sum(isfinite(earRem));
latRem = mean(remDFF(:, 3:5), 2, 'omitnan'); latRem_valid = sum(isfinite(latRem));
fprintf('  Remaining early vs late: n_early=%d, n_late=%d\n', earRem_valid, latRem_valid);

%% ========== STAGE 7: End-recruited neurons (Fig 3N) ==========
fprintf('\n========================================\n');
fprintf('STAGE 7: End-recruited neurons (Fig 3N)\n');
fprintf('========================================\n');
fprintf('Source: mergedEndDFF\n\n');

endVals = mergedEndDFF(:);
endVals_valid = endVals(isfinite(endVals) & endVals ~= 0);
fprintf('  Total end-recruited entries: %d\n', numel(mergedEndDFF));
fprintf('  Valid (non-zero, finite): %d\n', numel(endVals_valid));

%% ========== SUMMARY ==========
fprintf('\n========================================\n');
fprintf('SUMMARY: Full pipeline\n');
fprintf('========================================\n');
fprintf('Fig 3B raster:           %d unique neurons\n', size(fActBin, 1));
fprintf('  Active 5/5 weeks:      %d\n', sum(weekCounts_raster == 5));
fprintf('  Active 4/5 weeks:      %d\n', sum(weekCounts_raster == 4));
fprintf('  Active 4+/5 weeks:     %d\n', sum(weekCounts_raster >= 4));
fprintf('  Active 1/5 weeks:      %d\n', sum(weekCounts_raster == 1));
fprintf('Spatially tracked (PW):  %d ROIs\n', total_stage2);
fprintf('Activation-filtered:     %d response units\n', total_stage3);
fprintf('D-F stats (early+late):  %d\n', sum(isfinite(mean(mergedDFF(:,1:2),2,'omitnan')) & isfinite(mean(mergedDFF(:,3:5),2,'omitnan'))));
fprintf('Subset (high-intensity): %d\n', sum(subset_idx));
fprintf('Remaining:               %d\n', sum(remaining_idx));
fprintf('End-recruited:           %d\n', numel(endVals_valid));

fprintf('\n=== NOTE ===\n');
fprintf('The Fig 3B raster (3528) and D-F analysis (289/238) come from\n');
fprintf('PARALLEL pipelines. The 238 is NOT a subset of the 692 (4+/5 weeks).\n');
fprintf('Fig 3B uses weekly binary activation across any stim condition.\n');
fprintf('Fig 3D-F uses per-condition activation validation on spatially tracked ROIs.\n');
