%% Export recruitment-figure data so the Python plotter can render in matching style.
% Inputs: figures/.../final/recruited_traces.mat (already exists)
% Outputs (same dir): recruitment_data.mat, recruitment_traces.csv, recruitment_counts.csv

clear; clc;
WIN_FRAMES = 5;   % 3 = first ~101 ms (strict pre-movement); 5 = first ~169 ms
suffix     = sprintf('first%df', WIN_FRAMES);

scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fullfile(scriptDir, '..', '..');
outDir    = fullfile(repoRoot, 'figures', 'wheel_movement', 'first5_recruitment_control', 'final');
load(fullfile(outDir, 'recruited_traces.mat'), ...
     'traces','t_ms','thr','M_aidx','M_week','M_day','M_chan','M_row','M_cur');

% --- example pair (same indices the MATLAB fig uses) ---
earlyCross   = 38;     % thr ~0.73, clean crosser (peak 2.2)
earlyNoCross = 4000;   % thr ~0.47, late riser, misses first-5

trace_recruited     = traces(earlyCross,:);
trace_notrecruited  = traces(earlyNoCross,:);
thr_recruited       = thr(earlyCross);
thr_notrecruited    = thr(earlyNoCross);

% --- first-N window (N = WIN_FRAMES) ---
frameMs = 33.73;
[~, onsetCol] = min(abs(t_ms - 0));
first5Cols = onsetCol + (0:WIN_FRAMES-1);   % kept variable name "first5Cols" for back-compat in Python
first5_ms  = t_ms(first5Cols(end));

% --- count panel (POOL_CURRENT = true: # neurons recruited per ch x session) ---
fN_peak        = max(traces(:, onsetCol:onsetCol+WIN_FRAMES-1), [], 2, 'omitnan');
first5_exceeds = fN_peak > thr(:);

T = table; T.animal = M_aidx(:); T.day = M_day(:); T.week = M_week(:);
T.chan = M_chan(:); T.row = M_row(:); T.first5 = first5_exceeds(:);
[Gu, uA, uDay, uWk, uCh, ~] = findgroups(T.animal, T.day, T.week, T.chan, T.row);
unit_rec = splitapply(@any, T.first5, Gu);
recTab   = table(uA, uDay, uWk, uCh, double(unit_rec), ...
                 'VariableNames', {'animal_idx','day','week','chan','rec'});

[Gc, cA, cDay, cWk, cCh] = findgroups(recTab.animal_idx, recTab.day, recTab.week, recTab.chan);
n_rec = splitapply(@sum, recTab.rec, Gc);
countTable = table(cA, cDay, cWk, cCh, n_rec, ...
                  'VariableNames', {'animal_idx','day','week','chan','n_rec'});

earlyVals = countTable.n_rec(countTable.week <= 1);
lateVals  = countTable.n_rec(countTable.week >= 2);
p_value   = ranksum(earlyVals, lateVals);

% --- amplitude panel: weekly mean+/-error of max dF/F (first 5 frames) for the
%     increasing subset (Roy's allSubset, already filtered in subset_traces.mat).
S = load(fullfile(outDir, 'subset_traces.mat'), ...
         'traces','t_ms','M_aidx','M_row','M_chan','M_cur','M_week','M_day');
[~, ocS] = min(abs(S.t_ms - 0));
amp = max(S.traces(:, ocS:ocS+WIN_FRAMES-1), [], 2, 'omitnan');
keyCols = [S.M_aidx(:), S.M_row(:), S.M_chan(:), S.M_cur(:)];
[~,~,Gcond] = unique(keyCols, 'rows');
nC = max(Gcond);
amp_weeks = 0:4; nW = numel(amp_weeks);
[Gcw, cC, cW] = findgroups(Gcond, S.M_week(:));
mAmp = splitapply(@(x) mean(x,'omitnan'), amp, Gcw);
Mat = nan(nC, nW);
for i = 1:numel(mAmp)
  w = find(amp_weeks == cW(i));
  if ~isempty(w); Mat(cC(i), w) = mAmp(i); end
end

amp_mean = nan(1,nW); amp_std = nan(1,nW); amp_sem = nan(1,nW); amp_n = zeros(1,nW);
for w = 1:nW
  col = Mat(:,w); col = col(~isnan(col));
  amp_n(w)    = numel(col);
  amp_mean(w) = mean(col);
  amp_std(w)  = std(col);
  amp_sem(w)  = amp_std(w) / sqrt(max(amp_n(w),1));
end

% --- per channel-session amplitude for grouped early-vs-late panel ---
% each dot = one (animal, session, channel); value = mean amp across subset
% conditions at that channel-session.  Same unit as panel 2.
[Gcs, ~, ~, ~] = findgroups(S.M_aidx(:), S.M_day(:), S.M_chan(:));
amp_per_cs    = splitapply(@(x) mean(x,'omitnan'), amp,           Gcs);
cs_week       = splitapply(@(x) x(1),             S.M_week(:),    Gcs);
amp_chsess_early = amp_per_cs(cs_week <= 1);
amp_chsess_late  = amp_per_cs(cs_week >= 2);
amp_p_value      = ranksum(amp_chsess_early, amp_chsess_late);

% --- save .mat (-v7 -> scipy.io.loadmat compatible) ---
save(fullfile(outDir, sprintf('recruitment_data_%s.mat', suffix)), ...
     'WIN_FRAMES', ...
     't_ms','trace_recruited','trace_notrecruited', ...
     'thr_recruited','thr_notrecruited', ...
     'first5_ms','frameMs','onsetCol','first5Cols', ...
     'earlyVals','lateVals','p_value', ...
     'amp_weeks','amp_mean','amp_std','amp_sem','amp_n', ...
     'amp_chsess_early','amp_chsess_late','amp_p_value','-v7');

% --- CSVs (handy for inspection / external use) ---
T_ex = table(t_ms(:), trace_recruited(:), trace_notrecruited(:), ...
             'VariableNames', {'t_ms','recruited','not_recruited'});
writetable(T_ex,       fullfile(outDir, 'recruitment_traces.csv'));
writetable(countTable, fullfile(outDir, 'recruitment_counts.csv'));

fprintf('Saved recruitment_data.mat, recruitment_traces.csv, recruitment_counts.csv\n');
fprintf('thr_rec=%.3f, thr_nonrec=%.3f, p=%.3g  (early n=%d, late n=%d)\n', ...
        thr_recruited, thr_notrecruited, p_value, numel(earlyVals), numel(lateVals));
fprintf('amp_first5 weekly mean: '); fprintf('%.3f ', amp_mean); fprintf('\n');
fprintf('amp_first5 weekly n   : '); fprintf('%d ', amp_n); fprintf('\n');
