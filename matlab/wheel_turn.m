addpath('\\10.129.151.108\xieluanlabs\xl_stimulation\Robin\lib\neuroshare\')
nevFile = 'block1_H8M41_ephys1035.nev';
expTypeDir = 'C:/data/ICMS92/Behavior/30-Aug-2023/ChannelVolumetric';
[~, hFile] = ns_OpenFile(fullfile(expTypeDir,nevFile)); 
% [~, nsFileInfo] = ns_GetFileInfo(fullfile(source_dir,nevFile))

time_str = extractBetween(string(nevFile),'_','_ephys');
fileString = 'block1_H8M41';
matfileName = [fileString,'_behavior.mat'];
matFilePath = fullfile(expTypeDir,matfileName);
load(matFilePath,'expStruct');
depth_idx = expStruct.depthIndices(expStruct.channelIndex);
intan_ch = depth2intan(double(depth_idx));
trialData = expStruct.trialData;

if isempty(hFile.FileInfo(1).MemoryMap)
    return
end

all_stim_times = double(hFile.FileInfo(1).MemoryMap.Data.TimeStamp);
packetIDs = double(hFile.FileInfo(1).MemoryMap.Data.PacketID);
classIDs = double(hFile.FileInfo(1).MemoryMap.Data.Class);
if length(unique(packetIDs)) > 2
    fprintf('Ignoring spike data...\n');
    goodidx = find(packetIDs == 0 | packetIDs > 5000);
    packetIDs = packetIDs(goodidx);
    classIDs = classIDs(goodidx);
    all_stim_times = all_stim_times(goodidx);
end
classes = unique(classIDs);
img_id = find(strcmp({hFile.Entity.Label}, 'analog 1'));
rev_id = find(strcmp({hFile.Entity.Label}, 'analog 4'));
[~, ~, imgAnalogRaw] = ns_GetAnalogData(hFile, img_id, 1, hFile.TimeSpan);
[~, ~, revAnalogRaw] = ns_GetAnalogData(hFile, rev_id, 1, hFile.TimeSpan);

init_ts = find(diff(imgAnalogRaw) > 1000); % find where derivative is greater than 1000

dt = diff(init_ts); % samples between derivative threshold crossings 
img_ts = [];
for i = 1:numel(dt)
    if dt(i) > 1000 % assume 33 ms frame period 
        img_ts = [img_ts; init_ts(i)];
    end
end

frame_freq = 1/median(diff(img_ts)/3e4); % Hz
sprintf('Frame rate: %.2f\n',frame_freq);
% Extract imaging timestamps
tsVec = [];
group = [];
for i = 1:numel(classes)
    ts_i = all_stim_times(classIDs == classes(i));
    tsVec = [tsVec; ts_i];
    group = [group; i*ones(numel(ts_i),1)];
end
% 
tsVec = [tsVec; img_ts];
group = [group; 4*ones(numel(img_ts),1)];
assert(size(group,1)==size(tsVec,1));

[~,b] = sort(tsVec); % sort by time 
tsVec = tsVec(b)/3e4;
group = group(b);

trialStarts = tsVec(group == 3); % img start stop 
trialStarts = trialStarts(1:4:end);

% calculate trials by finding tsVec between trialStart and next trialStart
trials = [];
% remove timestamps before first trial 
goodidx = find(tsVec >= trialStarts(1));
tsVec = tsVec(goodidx);
group = group(goodidx);
assert(numel(tsVec) == numel(group));

for i = 1:numel(trialStarts)
    startTs = trialStarts(i);
    if i == numel(trialStarts)
        endTs = max(tsVec) + 0.001;
    else
        endTs = trialStarts(i+1);
    end

    idx = find((tsVec >= startTs) & (tsVec < endTs));
    n = numel(idx);
    x = i * ones(n,1);
    trials = [trials; x];
end

try 
    assert(numel(trials)==numel(group));
catch
    numel(trials)
    numel(group) 
    error('Different number of trials and groups for %s: %d trials vs %d groups',nevFile, numel(trials), numel(group))
end

nTRIALS = max(trialData(:,1));
dataCellArr = cell(nTRIALS,1);
trialDataStruct = struct();
for i = 1:nTRIALS

stim_ts = tsVec(trials == i & group == 1);
img_start_stop = tsVec(trials == i & group == 3);
stim_start_stop = tsVec(trials == i & group == 2);
img_ts = tsVec(trials == i & group == 4);

img_start_stop = img_start_stop(1:2:end);
img_ts = img_ts(1:end);
stim_start_stop = stim_start_stop(1:2:end);
stim_ts = stim_ts(1:end);

current = trialData(i,2);
hit = trialData(i,3);
rt = diff(stim_start_stop);
if rt >= 0.695
    rt = nan;
end

if ~isempty(stim_start_stop)
    if stim_start_stop(1) > max(img_start_stop) 
        good_trial = 0;
    else
        good_trial = 1;
    end
else
    good_trial = 0;
end

dataCellArr{i} = {current; hit; rt; ...
    good_trial; img_start_stop;img_ts;stim_start_stop;stim_ts};
trial_i_cell = dataCellArr{i};

% save as struct with names
names = {'Current', 'Response', 'ResponseTime', 'isGoodTrial', 
    'ImageStartStop', 'FrameTimestamps', 
    'StimStartStop', 'StimTimestamps'};

for j = 1:11
    trialDataStruct(i).(names{j}) = trial_i_cell{j};
end
end


%%
function sig_out = lowpass_filter(sig_in, cutoff, order)
    Wn = cutoff / (30000 / 2);
    [b, a] = butter(order, Wn, 'low');
    sig_out = filtfilt(b, a, sig_in);
end

plot_flag = false;
angle_thresholds = [];
good_segments = {};
good_response_times = [];
for i = 1:100
    trial = trialDataStruct(i); 
    isGoodTrial = trial.isGoodTrial;

    if isGoodTrial && trial.Current > 0
        stimStartStops = trial.StimStartStop;      
        stim = trial.StimTimestamps;    
        stimOnset = stimStartStops(1);
        responseTime = trial.ResponseTime * 30000;

        % Segment around stimulus
        segStart = round((stimStartStops(1) - 0.7) * 30000);
        segEnd = round((stim(end) + 0.7) * 30000);
        ref = segStart;

        % Extract and filter
        rev_seg = revAnalogRaw(segStart : segEnd);
        filt_rev_seg = lowpass_filter(rev_seg, 50, 4);

        % Compute sample index
        stimOnset_samples = round(stimOnset * 30000);
        index = round(stimOnset_samples + responseTime - ref);
        stimOnset_ms = (stimOnset_samples - ref) / fs * 1000;
        response_time = (stimOnset_samples + responseTime - ref) / fs * 1000 - stimOnset_ms;

        if index > 0 && index <= length(filt_rev_seg)
            threshold = filt_rev_seg(index);
            angle_thresholds(end+1) = threshold;
            if max(filt_rev_seg) < 1500
                good_segments{end+1} = filt_rev_seg;
                good_response_times(end+1) = response_time;
            end

        else
            warning('Index out of bounds for trial %d', i);
        end

        
        % Optional plot
        if plot_flag
            figure;
            fs = 30000;
        
            % Create time vector in ms, centered at stimulus onset
            stimOnset_ms = (stimOnset_samples - ref) / fs * 1000;
            x = (0:length(filt_rev_seg)-1) / fs * 1000 - stimOnset_ms;
        
            plot(x, filt_rev_seg)
        
            if ~isempty(stim)
                stim_lines = (stim * fs - ref) / fs * 1000 - stimOnset_ms;
                xline(stim_lines);
            end
        
            if ~isnan(responseTime)
                response_line = (stimOnset_samples + responseTime - ref) / fs * 1000 - stimOnset_ms;
                xline(response_line, 'r');
            end
        
            xlabel('Time (ms)');
            ylabel('Voltage (mV)');
            xlim([-700, 1400]);
            title(['Trial ' num2str(i)]);
        end
    end
end

%%
colors = parula(numel(good_segments)+2);

for i = 1:numel(good_segments)
    segment = good_segments{i}; 
    segment = segment - min(segment);
    response_time = good_response_times(i);
    fs = 30000;
    x = (0:length(segment)-1) / fs * 1000 - stimOnset_ms;
    
    plot(x, segment, 'Color', colors(i, :));
    xline(response_time);
    hold on


    

end

xlabel('Time (ms)');
ylabel('Voltage (mV)');
xlim([-700, 1400]);
y = 90;
% plot([0, 700], [y, y], 'k-', 'LineWidth', 3);
ylim([0, 100]);

title('ICMS92 hit trials')