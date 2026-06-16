%% Export trial-averaged dF/F traces for the INCREASING-SUBSET conditions
% Parallel to export_recruited_traces.m, but the population is Roy's predefined
% increasing subset (allSubset / allChanData), NOT the full-window-recruited pool.
%   allChanData{a} cols = [roiRow (tracked roiMap row), channel, current]
%   allSubset{a}        = logical, true for the increasing subset
% For each subset condition x session:
%   - trial-averaged, span-5 smoothed dF/F over PRE+MID+POST, aligned so t=0 is
%     stim onset (pre = negative time)
%   - amp_full   = max of the MID (stim-window) trace
%   - amp_first5 = max of the first 5 MID frames (~169 ms; pre-movement window)
% Identity for dedup is the tracked roiMap row (M_row), stable across sessions.
%
% Output: figures/.../final/subset_traces.mat (+ _meta.csv)

clear; clc;
sourceMat = 'S:\Roy\NoLimitsVariablesV2_ROIPLotter2Data.mat';
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fullfile(scriptDir, '..', '..');
outDir    = fullfile(repoRoot, 'figures', 'wheel_movement', 'first5_recruitment_control', 'final');
if ~exist(outDir,'dir'); mkdir(outDir); end
S = load(sourceMat,'allROI2Comp','allDays','allSubset','allChanData');
allROI2Comp=S.allROI2Comp; allDays=S.allDays; allSubset=S.allSubset; allChanData=S.allChanData;

animals = struct([]);
animals(1).name='ICMS92'; animals(1).files={'S:\ICMS92\9-6-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS92\9-8-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS92\9-12-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS92\9-14-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS92\9-19-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS92\9-21-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS92\9-25-23\PlanarTrialsG8\SessionMetrics.mat'};
animals(2).name='ICMS98'; animals(2).files={'S:\ICMS98\10-20-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS98\10-24-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS98\10-26-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS98\10-31-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS98\11-2-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS98\11-7-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS98\11-17-23\PlanarTrialsG8\SessionMetrics.mat'};
animals(3).name='ICMS93'; animals(3).files={'S:\ICMS93\8-30-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS93\9-6-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS93\9-14-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS93\9-22-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS93\9-26-23\PlanarTrialsG8\SessionMetrics.mat'};

frameMs=33.73; KPRE=20;   % frames of pre-stim baseline to keep (~675 ms)

trCell={}; amp_full=[]; amp_first5=[]; M_aidx=[]; M_week=[]; M_day=[]; M_chan=[]; M_cur=[]; M_row=[]; M_roi=[];

for a=1:numel(animals)
    fprintf('=== %s ===\n', animals(a).name);
    roiMap=allROI2Comp{a}; days=allDays{a}; weeks=ceil(days./7); nSess=numel(animals(a).files);
    cd_=allChanData{a}; sub=logical(allSubset{a}(:));
    valid = cd_(:,1)>0 & ~isnan(cd_(:,1)) & sub;   % subset conditions only
    conds = cd_(valid,:);                           % [roiRow, channel, current]
    fprintf('  %d subset conditions\n', size(conds,1));

    for s=1:nSess
        sess=load(animals(a).files{s},'allTraces','TrialInfo','regions');
        traces=sess.allTraces; trialInfo=sess.TrialInfo; nRoi=numel(sess.regions);
        nTr=min(100,size(trialInfo,1));
        completeTrace=[];
        for t=1:nTr; if trialInfo(t,6)~=0; completeTrace=[completeTrace,[traces{t,1},traces{t,2},traces{t,3}]]; end; end %#ok<AGROW>
        if isempty(completeTrace); continue; end
        regionMeans=mean(completeTrace,2);

        for k=1:size(conds,1)
            roiRow=conds(k,1); ch=conds(k,2); cu=conds(k,3);
            if roiRow<1 || roiRow>size(roiMap,1); continue; end
            roiId=roiMap(roiRow,s);
            if roiId<1 || roiId>nRoi; continue; end
            inds=find(trialInfo(:,1)==ch & trialInfo(:,3)==cu & trialInfo(:,6)~=0);
            if isempty(inds); continue; end
            Full=[]; Mid=[];
            for i=1:numel(inds)
                t=inds(i); pre=traces{t,1}; mid=traces{t,2}; post=traces{t,3};
                if isempty(mid); continue; end
                F=[pre(roiId,:), mid(roiId,:), post(roiId,:)];
                dff=smooth((F-regionMeans(roiId))./regionMeans(roiId),5)';
                onset=size(pre,2)+1; nmid=size(mid,2);
                startI=onset-KPRE;
                if startI>=1; v=dff(startI:end); else; v=[nan(1,1-startI), dff]; end
                Full=stack_pad(Full, v);
                Mid =stack_pad(Mid,  dff(onset:onset+nmid-1));
            end
            if isempty(Full); continue; end
            curMid=mean(Mid,1,'omitnan'); trace=mean(Full,1,'omitnan');
            trCell{end+1}=trace; %#ok<SAGROW>
            amp_full(end+1)  = max(curMid,[],'omitnan');
            amp_first5(end+1)= max(curMid(1:min(5,end)),[],'omitnan');
            M_aidx(end+1)=a; M_week(end+1)=weeks(s); M_day(end+1)=days(s);
            M_chan(end+1)=ch; M_cur(end+1)=cu; M_row(end+1)=roiRow; M_roi(end+1)=roiId;
        end
    end
end

T=max(cellfun(@numel,trCell));
traces=nan(numel(trCell),T);
for i=1:numel(trCell); v=trCell{i}; traces(i,1:numel(v))=v; end
t_ms=((1:T)-(KPRE+1))*frameMs;   % t=0 at stim onset (column KPRE+1); pre = negative
amp_full=amp_full(:); amp_first5=amp_first5(:); animal=string({animals(M_aidx).name})';

save(fullfile(outDir,'subset_traces.mat'),'traces','t_ms','amp_full','amp_first5','animal', ...
     'M_aidx','M_week','M_day','M_chan','M_cur','M_row','M_roi','KPRE','-v7.3');
meta=table(animal, M_aidx', M_week', M_day', M_chan', M_cur', M_row', M_roi', amp_full, amp_first5, ...
     'VariableNames',{'animal','animal_idx','week','day','chan','current','roiRow','roiId','amp_full','amp_first5'});
writetable(meta, fullfile(outDir,'subset_traces_meta.csv'));

early=M_week(:)<=1; late=M_week(:)>=2;
fprintf('\nN subset traces: %d (early %d, late %d); T=%d frames, t in [%.0f, %.0f] ms\n',...
    numel(trCell), sum(early), sum(late), T, t_ms(1), t_ms(end));
fprintf('amp_full   mean: early %.3f, late %.3f\n', mean(amp_full(early)), mean(amp_full(late)));
fprintf('amp_first5 mean: early %.3f, late %.3f\n', mean(amp_first5(early)), mean(amp_first5(late)));
fprintf('Saved subset_traces.mat + subset_traces_meta.csv\n');

function M=stack_pad(M,v)
    v=v(:)'; if isempty(M); M=v; return; end
    w=max(size(M,2),numel(v));
    if size(M,2)<w; M(:,end+1:w)=NaN; end
    if numel(v)<w; v(end+1:w)=NaN; end
    M=[M;v];
end
