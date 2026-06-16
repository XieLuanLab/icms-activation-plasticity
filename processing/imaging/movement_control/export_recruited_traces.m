%% Export full trial-averaged dF/F traces (pre+mid+post) for recruited ROIs
% For each tracked ROI x stim channel (4-6 uA pooled), per session:
%   - trial-averaged, span-5 smoothed dF/F over PRE + MID + POST, aligned so
%     t=0 is stim onset (pre = negative time)
%   - per-ROI threshold thr = 2*CV (CV = regionSTD/regionMean, full-trial)
%   - recruitment flags: full-window gate (>5 consec mid frames >2*CV, onset<700)
%     and first-5 gate (any of mid frames 1-5 >2*CV)
% Stores ROIs recruited by the FULL-window gate. Plot the mean trace and overlay
% mean(thr) as the average activation threshold line.
%
% Output: figures/.../final/recruited_traces.mat (+ _meta.csv)
%   traces [N x T]  trial-avg dF/F (pre+mid+post), NaN-padded, aligned at onset
%   full_trial_traces {N x 1} per-trial aligned full dF/F traces for each row
%   t_ms   [1 x T]  time from stim onset (ms); negative = pre-stim
%   thr    [N x 1]  2*CV activation threshold (dF/F units) for that ROI
%   week_* arrays: ready-to-plot weekly mean/SEM traces from traces
%   meta: animal, M_aidx, M_week, M_day, M_chan, M_roi, M_rf5

clear; clc;
sourceMat = 'S:\Roy\NoLimitsVariablesV2_ROIPLotter2Data.mat';
scriptDir = fileparts(mfilename('fullpath'));
repoRoot  = fullfile(scriptDir, '..', '..');
outDir    = fullfile(repoRoot, 'figures', 'wheel_movement', 'first5_recruitment_control', 'final');
if ~exist(outDir,'dir'); mkdir(outDir); end
S = load(sourceMat,'allROI2Comp','allDays'); allROI2Comp=S.allROI2Comp; allDays=S.allDays;

animals = struct([]);
animals(1).name='ICMS92'; animals(1).chans=[9 11 12];
animals(1).files={'S:\ICMS92\9-6-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS92\9-8-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS92\9-12-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS92\9-14-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS92\9-19-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS92\9-21-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS92\9-25-23\PlanarTrialsG8\SessionMetrics.mat'};
animals(2).name='ICMS98'; animals(2).chans=[4 5 6];
animals(2).files={'S:\ICMS98\10-20-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS98\10-24-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS98\10-26-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS98\10-31-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS98\11-2-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS98\11-7-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS98\11-17-23\PlanarTrialsG8\SessionMetrics.mat'};
animals(3).name='ICMS93'; animals(3).chans=[9 11 12];
animals(3).files={'S:\ICMS93\8-30-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS93\9-6-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS93\9-14-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS93\9-22-23\PlanarTrialsG8\SessionMetrics.mat';'S:\ICMS93\9-26-23\PlanarTrialsG8\SessionMetrics.mat'};

analysisCurrs=[4 5 6]; frameMs=33.73; KPRE=20;   % frames of pre-stim baseline to keep (~675 ms)

trCell={}; fullTrialCell={}; thr=[]; M_aidx=[]; M_week=[]; M_day=[]; M_chan=[]; M_cur=[]; M_roi=[]; M_row=[]; M_rf5=[];

for a=1:numel(animals)
    fprintf('=== %s ===\n', animals(a).name);
    roiMap=allROI2Comp{a}; days=allDays{a}; weeks=ceil(days./7); nSess=numel(animals(a).files);
    for s=1:nSess
        sess=load(animals(a).files{s},'allTraces','TrialInfo','regions');
        traces=sess.allTraces; trialInfo=sess.TrialInfo; nRoi=numel(sess.regions);
        nTr=min(100,size(trialInfo,1));
        completeTrace=[];
        for t=1:nTr; if trialInfo(t,6)~=0; completeTrace=[completeTrace,[traces{t,1},traces{t,2},traces{t,3}]]; end; end %#ok<AGROW>
        if isempty(completeTrace); continue; end
        regionMeans=mean(completeTrace,2); regionSTDs=std(completeTrace,[],2);
        validRows=find(roiMap(:,s)>0 & roiMap(:,s)<=nRoi);

        for ch=animals(a).chans
          for cu=analysisCurrs           % per-current condition (no pooling across 4-6 uA)
            inds=find(trialInfo(:,1)==ch & trialInfo(:,3)==cu & trialInfo(:,6)~=0);
            if isempty(inds); continue; end
            for ii=1:numel(validRows)
                roiRow=validRows(ii);            % stable tracked-neuron identity (roiMap row)
                roiId=roiMap(validRows(ii),s);   % session-local region label
                curStd=regionSTDs(roiId)/regionMeans(roiId);
                if ~isfinite(curStd)||curStd<=0; continue; end
                Full=[]; Mid=[];
                for i=1:numel(inds)
                    t=inds(i); pre=traces{t,1}; mid=traces{t,2}; post=traces{t,3};
                    if isempty(mid); continue; end
                    F=[pre(roiId,:), mid(roiId,:), post(roiId,:)];
                    dff=smooth((F-regionMeans(roiId))./regionMeans(roiId),5)';
                    onset=size(pre,2)+1; nmid=size(mid,2);
                    % aligned full trace: KPRE frames before onset .. end
                    startI=onset-KPRE;
                    if startI>=1; v=dff(startI:end);
                    else; v=[nan(1,1-startI), dff]; end       % pad front if pre<KPRE
                    Full=stack_pad(Full, v);
                    Mid=stack_pad(Mid, dff(onset:onset+nmid-1));   % mid only for the gate
                end
                if isempty(Full); continue; end
                curMid=mean(Mid,1,'omitnan'); trace=mean(Full,1,'omitnan');
                aTP=curMid>2*curStd;
                rec_full = max_run(aTP)>5 && ~isempty(find(aTP,1)) && find(aTP,1)*frameMs<700;
                rec_f5   = any(curMid(1:min(5,end))>2*curStd);
                if rec_full
                    trCell{end+1}=trace; fullTrialCell{end+1}=Full; thr(end+1)=2*curStd; %#ok<SAGROW>
                    M_aidx(end+1)=a; M_week(end+1)=weeks(s); M_day(end+1)=days(s);
                    M_chan(end+1)=ch; M_cur(end+1)=cu; M_roi(end+1)=roiId; M_row(end+1)=roiRow; M_rf5(end+1)=rec_f5;
                end
            end
          end
        end
    end
end

T=max(cellfun(@numel,trCell));
traces=nan(numel(trCell),T);
for i=1:numel(trCell); v=trCell{i}; traces(i,1:numel(v))=v; end
full_trial_traces=cell(numel(fullTrialCell),1);
for i=1:numel(fullTrialCell); full_trial_traces{i}=pad_cols(fullTrialCell{i},T); end
t_ms=((1:T)-(KPRE+1))*frameMs;   % t=0 at stim onset (column KPRE+1); pre = negative
thr=thr(:); animal=string({animals(M_aidx).name})';

week_values=unique(M_week(:))';
week_mean_traces=nan(numel(week_values),T);
week_sem_traces=nan(numel(week_values),T);
week_n_traces=zeros(numel(week_values),1);
week_thr_mean=nan(numel(week_values),1);
for wi=1:numel(week_values)
    idx=M_week(:)==week_values(wi);
    week_n_traces(wi)=sum(idx);
    week_mean_traces(wi,:)=mean(traces(idx,:),1,'omitnan');
    week_sem_traces(wi,:)=std(traces(idx,:),0,1,'omitnan') ./ sqrt(max(week_n_traces(wi),1));
    week_thr_mean(wi)=mean(thr(idx),'omitnan');
end

save(fullfile(outDir,'recruited_traces.mat'),'traces','full_trial_traces','t_ms','thr','animal', ...
     'M_aidx','M_week','M_day','M_chan','M_cur','M_roi','M_row','M_rf5','KPRE', ...
     'week_values','week_mean_traces','week_sem_traces','week_n_traces','week_thr_mean','-v7.3');
meta=table(animal, M_aidx', M_week', M_day', M_chan', M_cur', M_roi', M_row', logical(M_rf5'), thr, ...
     'VariableNames',{'animal','animal_idx','week','day','chan','current','roiId','roiRow','rec_first5','thr_2CV'});
writetable(meta, fullfile(outDir,'recruited_traces_meta.csv'));

early=M_week(:)<=1; late=M_week(:)>=2;
me=mean(traces(early,:),1,'omitnan'); ml=mean(traces(late,:),1,'omitnan');
[~,pe]=max(me); [~,pl]=max(ml);
fprintf('\nN recruited traces: %d (early %d, late %d); T=%d frames, t in [%.0f, %.0f] ms\n',...
    numel(trCell), sum(early), sum(late), T, t_ms(1), t_ms(end));
fprintf('Mean 2*CV threshold: %.3f dF/F (early %.3f, late %.3f)\n', mean(thr), mean(thr(early)), mean(thr(late)));
fprintf('Population-mean dF/F peak: early %.0f ms, late %.0f ms\n', t_ms(pe), t_ms(pl));
fprintf('Saved recruited_traces.mat + recruited_traces_meta.csv\n');
fprintf('Weekly means are in week_values, week_mean_traces, week_sem_traces, and week_n_traces\n');

function M=stack_pad(M,v)
    v=v(:)'; if isempty(M); M=v; return; end
    w=max(size(M,2),numel(v));
    if size(M,2)<w; M(:,end+1:w)=NaN; end
    if numel(v)<w; v(end+1:w)=NaN; end
    M=[M;v];
end
function M=pad_cols(M,T)
    if size(M,2)<T; M(:,end+1:T)=NaN; end
    if size(M,2)>T; M=M(:,1:T); end
end
function r=max_run(tp)
    tp=tp(:)'; tp(~isfinite(tp))=false; tp=logical(tp);
    if ~any(tp); r=0; return; end
    d=[true,diff(tp)~=0]; vals=tp(d); runs=diff(find([d,true])); runs=runs(vals==1); r=max(runs);
end
