% Requires raw session data (not in repo).
% Internal paths must be configured manually for your local environment.

addpath(genpath('S:/Roy/ryhattori-PatchWarp-1.3.3.0'));

% calculate trace means and stds
% TODO - report baseline trace means over each session, i.e. see if there is a drifting baseline in each session or for each roi
figGen=0;
    

fprintf('\n=== SWEEPING ALL ANIMALS IN MASTER PIPELINE ===\n');
for animal =  1:1
target_current = 5;  % <--- EDIT THIS TO CHANGE STIM CURRENT (e.g. 4, 5, 6)
% target_chan explicitly picks the 2nd structural channel automatically per animal below! 

     if(animal==1)
    icms = 'ICMS92';

    files2Load = {'S:\ICMS92\9-6-23\PlanarTrialsG8\SessionMetrics.mat';...
                'S:\ICMS92\9-8-23\PlanarTrialsG8\SessionMetrics.mat';...
                'S:\ICMS92\9-12-23\PlanarTrialsG8\SessionMetrics.mat';...
                'S:\ICMS92\9-14-23\PlanarTrialsG8\SessionMetrics.mat';...
                'S:\ICMS92\9-19-23\PlanarTrialsG8\SessionMetrics.mat';...
                'S:\ICMS92\9-21-23\PlanarTrialsG8\SessionMetrics.mat';...
                'S:\ICMS92\9-25-23\PlanarTrialsG8\SessionMetrics.mat'};
            
    minSetCnt = 4;
    refInd = 3;
    chans = [9,11,12]; % hardcoded for this animal - update for each animal
    daysTrained = [7,9,13,15,20,22,26];
    
    elseif(animal==2)
        icms = 'ICMS98';

        files2Load = {'S:\ICMS98\10-20-23\PlanarTrialsG8\SessionMetrics.mat';...
                      'S:\ICMS98\10-24-23\PlanarTrialsG8\SessionMetrics.mat';...
                      'S:\ICMS98\10-26-23\PlanarTrialsG8\SessionMetrics.mat';...
                      'S:\ICMS98\10-31-23\PlanarTrialsG8\SessionMetrics.mat';...
                      'S:\ICMS98\11-2-23\PlanarTrialsG8\SessionMetrics.mat';...
                      'S:\ICMS98\11-7-23\PlanarTrialsG8\SessionMetrics.mat';...
                      'S:\ICMS98\11-17-23\PlanarTrialsG8\SessionMetrics.mat'};
        refInd = 5;
        minSetCnt = 3;
        chans = [4,5,6]; % hardcoded for this animal - update for each animal
        daysTrained = [0,4,6,11,13,18,28];


    elseif(animal==3)
        icms = 'ICMS93';

        files2Load = {'S:\ICMS93\8-30-23\PlanarTrialsG8\SessionMetrics.mat';...
                      'S:\ICMS93\9-6-23\PlanarTrialsG8\SessionMetrics.mat';...
                      'S:\ICMS93\9-14-23\PlanarTrialsG8\SessionMetrics.mat';...
                      'S:\ICMS93\9-22-23\PlanarTrialsG8\SessionMetrics.mat';...
                      'S:\ICMS93\9-26-23\PlanarTrialsG8\SessionMetrics.mat'};
        refInd = 2;
        minSetCnt = 3;
        chans = [9,11,12]; % hardcoded for this animal - update for each animal
        daysTrained = [0,7,15,23,27];     
    end
     % Automatically lock onto the primary middle stimulation channel for this specific animal!
     target_chan = chans(2);

    % Instantly target pre-existing files in root `planar_analysis`!
    curFld = fullfile('..', '..', 'planar_analysis', icms, filesep); 
    if ~exist(curFld, 'dir'), mkdir(curFld); end
    
    srcFld = fullfile(curFld, 'Source', filesep); 
    if ~exist(srcFld, 'dir'), mkdir(srcFld); end


    %% perform initial pass of alignment of session data (if not already performed)
    alignFile = strcat(curFld,icms,'_alignment.mat');
    refData = load(files2Load{refInd});
    refMean = refData.meanPre./max(refData.meanPre,[],'all'); 
    if(~exist(alignFile, 'file'))
        % perform alignment
        ptForms = cell(numel(files2Load),1);
        for f = 1:numel(files2Load) 
            curData = load(files2Load{f});
            image2_mean = curData.meanPre./max(curData.meanPre,[],'all'); 
            image1_mean = refMean;

            % Normalize image intensity
            image1_mean = 1000*reshape(normalize(image1_mean(:), 'range'), size(image1_mean));
            image2_mean = 1000*reshape(normalize(image2_mean(:), 'range'), size(image2_mean));

%             if(f==5)
%                 image2_mean(1:200,:)=0.0001;
%             end
            
            transform1 = 'euclidean';
            transform2 = 'affine';
            warp_blocksize = 6;
            warp_overlap_pix_frac = 0.25;
            norm_radius = 0;     % Set to 0 when the signals are sparse or the result does not look good.
            alignVal = inf;
            % run 3 iterations (instead of 30) for a massive 10x speedup!
            for t = 1:3

                patchwarp_results = patchwarp_across_sessions(image1_mean, image2_mean,transform1, transform2, warp_blocksize, warp_overlap_pix_frac, norm_radius);
                shiftedplane = spatial_interp_patchwarp(image2_mean, patchwarp_results.warp1_cell{1}, 'euclidean', 1:512, 1:512);

                FOV = ones(512,512);
                shiftedFOV = spatial_interp_patchwarp(FOV, patchwarp_results.warp1_cell{1}, 'euclidean', 1:512, 1:512);

                refComp = image1_mean./max(image1_mean,[],'all');
                curShift = shiftedplane./max(shiftedplane,[],'all');
                
                % Trim comparision volumes
                compVol = refComp - curShift;
                compVol(shiftedFOV==0)=NaN;
                compCnt = sum(shiftedFOV,'all','omitnan');
                curAV = abs(mean(compVol,'all','omitnan'));
                curAV = curAV/compCnt; % normalize by volume compared


                if(curAV < alignVal)
                    alignVal = curAV;
                    ptForms{f} = patchwarp_results.warp1_cell{1};
                end
            end        
            if(figGen==1)
                curFig = figure('NumberTitle', 'off', 'Name', strcat('Sess',num2str(refInd),'--vs--Sess',num2str(f)));
                shiftedplane = spatial_interp_patchwarp(image2_mean, ptForms{f}, 'euclidean', 1:512, 1:512);
                imshow(imfuse(shiftedplane./max(shiftedplane,[],'all'),image1_mean./max(image1_mean,[],'all'),'falsecolor','Scaling','joint','ColorChannels','red-cyan'));
                txt = strcat('Sess',num2str(refInd),'vsSess',num2str(f));
                % combineFileName = strcat(curFld,'Alignment-',txt,'.fig');
                % saveas(curFig,combineFileName);  
                % combineFileName = strcat(curFld,'Alignment-',txt,'.tif');
                % saveas(curFig,combineFileName);
                % close(curFig)
                % 
                curFig = figure('NumberTitle', 'off', 'Name', strcat('Sess',num2str(refInd),'--vs--Sess',num2str(f)));
                shiftedplane = spatial_interp_patchwarp(image2_mean, ptForms{f}, 'euclidean', 1:512, 1:512);
                imshowpair(shiftedplane./max(shiftedplane,[],'all'),image1_mean./max(image1_mean,[],'all'),'montage');
                txt = strcat('Sess',num2str(refInd),'vsSess',num2str(f));
                % combineFileName = strcat(curFld,'AlignmentM-',txt,'.fig');
                % saveas(curFig,combineFileName);  
                % combineFileName = strcat(curFld,'AlignmentM-',txt,'.tif');
                % saveas(curFig,combineFileName);
                % close(curFig)
            end
        end

        %save alignment
        % save(alignFile,'ptForms');

    else

        % Alignment already performed, load alignment for current datasets
        temp = load(alignFile);
        ptForms = temp.ptForms;
    end


    
    
    %% Load electrode position for each session and channel
    chanPositions = NaN(3,2);
%     for f = 1:numel(files2Load) 
        curData = load(files2Load{refInd});
        selpath = files2Load{refInd};
        idcs = strfind(selpath,'\');
        pullPath = selpath(1:idcs(numel(idcs)-1)-1);
        ePath = strcat(pullPath,'\ElecPos.mat');
        if exist(ePath, 'file')
            elecDat = load(ePath);
            electDat = elecDat.contactPos;
            for ch = 1:numel(chans)
                chanPositions(ch,:) = electDat(chans(ch),1:2);
            end
        end

    AllData = cell(numel(files2Load),1);
    for s = 1:numel(files2Load)
        AllData{s} = load(files2Load{s});
    end

    %% create updated region positions for each dataset
    updatedRegions = cell(numel(files2Load),1);
    processedRegions =  cell(numel(files2Load),1);
    shiftedSum = zeros(512,512);
    shiftedSumV2 = zeros(512,512);
    for f = 1:numel(files2Load) 
        curData = AllData{f};
        regions = curData.regions;
        regionsV2 = cell(numel(regions),1);

        for r = 1:numel(regions)
            zProjectedBW = false(512,512);
            curPixels = regions(r).PixelList;
            for m = 1:size(curPixels,1)
                zProjectedBW(curPixels(m,2),curPixels(m,1))=1;
            end
            shiftedplane = spatial_interp_patchwarp(zProjectedBW, ptForms{f}, 'euclidean', 1:512, 1:512);
            shiftedplane = shiftedplane>0;

            if(sum(shiftedplane,'all')>1)
                regionsV2{r} = find(shiftedplane);
                shiftedSum = shiftedSum + shiftedplane;
            end
        end
        updatedRegions{f} = regionsV2;

        pR = false(numel(regions,1)); % mark invalid regions as processed
        for r = 1:numel(regions)
            pR(r) = isempty(regionsV2);
        end
        processedRegions{f} = pR; 
    end
    processedRegionsOrig = processedRegions;
    processedRegionsEnd = processedRegions;
    % Plot overlapping regions
    if(figGen==1)
        curFig = figure();
        imshow(shiftedSum./max(shiftedSum,[],'all'))
        % combineFileName = strcat(curFld,'OverlapMap.fig');
        % saveas(curFig,combineFileName);  
        % combineFileName = strcat(curFld,'OverlapMap.tif');
        % saveas(curFig,combineFileName);
        % close(curFig)
    end

  

    %% find overlapping regions across sessions
    roiSets = cell(1,1);
    roiSetsPx = cell(1,1);
    rsC = 0;
    for f = 1:numel(files2Load)
        curRegions = updatedRegions{f};
        curProc = processedRegions{f};
        for r1 = 1:numel(curRegions) % iterate through all regions and find subsequent regions that overlap
            if(curProc(r1)==0) % only process valid regions that have not been processed yet
                % initialze current set
                curSet = zeros(numel(files2Load),1);
                curSet(f) = r1;

                for f2 = f+1:numel(files2Load)
                    searchRegions = updatedRegions{f2};
                    searchProc = processedRegions{f2};
                    overlap = 1;
                    optInd = 0;
                    for r2 = 1:numel(searchRegions)
                        if(searchProc(r2)==0)
                            oCnt = ismember(curRegions{r1},searchRegions{r2});
                            if(sum(oCnt)>overlap) % find largest region of overlap
                                overlap = oCnt;
                                optInd = r2;
                            end
                        end
                    end

                    if(optInd>0) % if these was an overlap add to current set and make region as processed
                        curSet(f2) = optInd;
                        searchProc(optInd) = 1;
                        processedRegions{f2} = searchProc;
                    end
                end

                % add current set of ROIs to master list if there are sufficent
                % number of multiple sessions for the ROI
                setMask = curSet>0;
                if(sum(setMask) >= minSetCnt)
                    rsC = rsC+1;
                    roiSets{rsC} = curSet;
                    roiSetsPx{rsC} = curRegions{r1};
                end
            end
        end
    end

    
        
    %% Plot aligned datasets and selected ROIs 
    if(figGen==1)
        for f = 1:numel(files2Load) 
            curData = AllData{f};
            curRegions = curData.regions;
            image2_mean = curData.meanPre./max(curData.meanPre,[],'all'); 
    %         shiftedplane = spatial_interp_patchwarp(image2_mean, ptForms{f}, 'euclidean', 1:512, 1:512);

            curFig = figure('NumberTitle', 'off', 'Name', strcat('Sess',num2str(f)),'WindowState','maximized');
            imshow(image2_mean./max(image2_mean,[],'all'))
            hold on
            for r = 1:numel(roiSets)
                curROIInds = roiSets{r};
                curSessROIind = curROIInds(f);
                if(curSessROIind>0)
                    zProjectedBW = false(512,512);

                    curPixels = curRegions(curSessROIind).PixelList;
                    for m = 1:size(curPixels,1)
                        zProjectedBW(curPixels(m,2),curPixels(m,1))=1;
                    end
                    ry = mean(curPixels(:,1));
                    rx = mean(curPixels(:,2));
                    boundaries = bwboundaries(zProjectedBW);
                    for k=1:numel(boundaries)
                        b = boundaries{k};
                        plot(b(:,2),b(:,1),'r','LineWidth',1);
                    end
    %                 text(ry,rx,num2str(r),'Color','y')
                end
            end

            txt = strcat('Sess',num2str(f),'ROI');
    %         combineFileName = strcat(curFld,'AlignmentROIs-',txt,'.fig');
    %         saveas(curFig,combineFileName);  
            % combineFileName = strcat(curFld,'AlignmentROIs-',txt,'.tif');
            % saveas(curFig,combineFileName);
            % close(curFig)
        end
    end
    
    % Plot 
    if(figGen==1)
        for f = 1:numel(files2Load) 
            curData = AllData{s};
            curRegions = curData.regions;
            image2_mean = curData.meanPre./max(curData.meanPre,[],'all'); 
            shiftedplane = spatial_interp_patchwarp(image2_mean, ptForms{f}, 'euclidean', 1:512, 1:512);

            curFig = figure('NumberTitle', 'off', 'Name', strcat('Sess',num2str(f)),'WindowState','maximized');
            imshow(shiftedplane./max(shiftedplane,[],'all'))
            hold on
            for r = 1:numel(roiSets)
                
                % Plot consistent ROIs
                curROIInds = roiSets{r};
                if(curROIInds(f)~=0)
                    fInd = find(curROIInds,1);
                    curSessROIind = curROIInds(fInd);
                    zProjectedBW = false(512,512);
                    curData = AllData{fInd};
                    
                    curRegions2 = curData.regions;
            
                    curPixels = curRegions2(curSessROIind).PixelList;
                    for m = 1:size(curPixels,1)
                        zProjectedBW(curPixels(m,2),curPixels(m,1))=1;
                    end
                    zProjectedBW = spatial_interp_patchwarp(zProjectedBW, ptForms{fInd}, 'euclidean', 1:512, 1:512);
                    ry = mean(curPixels(:,1));
                    rx = mean(curPixels(:,2));
                    boundaries = bwboundaries(zProjectedBW);
                    for k=1:numel(boundaries)
                        b = boundaries{k};
                        plot(b(:,2),b(:,1),'r','LineWidth',1);
                    end
                end
                % Plot non-consistent ROIs
            end

            txt = strcat('Sess',num2str(f));
            % combineFileName = strcat(curFld,'AlignmentTransformedROIs-',txt,'.fig');
            % saveas(curFig,combineFileName);  
            % combineFileName = strcat(curFld,'AlignmentTransformedROIs-',txt,'.tif');
            % saveas(curFig,combineFileName);
            % close(curFig)
        end
    end
    
    
    
    
    
    % %% Generate averaged ROIs for showing alignment
    % aveROIImg = false(512,512,numel(roiSets));
    % for r = 1:numel(roiSets)
    %     curAve = zeros(512,512);
    %     % Plot consistent ROIs
    %     curROIInds = roiSets{r};
    %     n = 0;
    %     for f = 1:numel(files2Load)
    %         if(curROIInds(f)~=0)
    %             zProjectedBW = false(512,512);
    %             curData = load(files2Load{f});
    %             curRegions2 = curData.regions;
    % 
    %             curPixels = curRegions2(curROIInds(f)).PixelList;
    %             for m = 1:size(curPixels,1)
    %                 zProjectedBW(curPixels(m,2),curPixels(m,1))=1;
    %             end
    %             zProjectedBW = spatial_interp_patchwarp(zProjectedBW, ptForms{f}, 'euclidean', 1:512, 1:512);
    %             curAve = curAve+zProjectedBW;
    %             boundaries = bwboundaries(zProjectedBW);
    %             n=n+1;
    %         end
    %     end
    %     aveROIImg(:,:,r) = (curAve/n) >= 0.33;
    % end
    
    















    

    %% calculate trace means and stds
    allCurrs = [];
    regionMeans = cell(numel(files2Load),1);
    regionSTDs = cell(numel(files2Load),1);
    for s = 1:numel(files2Load)
        curData = AllData{s};
        trialInfo = curData.TrialInfo;
        currs = unique(trialInfo(:,3));
        regions = curData.regions;
        allCurrs = [allCurrs; currs];
        cN = 0;
        cPN = 0;
        
        curAllTraces = curData.allTraces;
        cleanAllTraces = cell(size(curAllTraces));
        completeTrace = zeros(numel(regions),1);
        completePreTrace = zeros(numel(regions),1);
        for t = 1:100
            if(trialInfo(t,6)~=0)
                temp = [curAllTraces{t,1},curAllTraces{t,2},curAllTraces{t,3}];
                completeTrace(:,cN+1:cN+size(temp,2)) = temp;
                completePreTrace(:,cPN+1:cPN+size(curAllTraces{t,1},2)) = curAllTraces{t,1};
                cN = cN + size(temp,2);
                cPN = cPN + size(curAllTraces{t,1},2);
            end
        end
        if(cN>0)
            regionMeans{s} = mean(completeTrace,2);
            regionSTDs{s} = std(completeTrace,[],2);
        end
    end
    allCurrs = unique(allCurrs);
    allCurrs(allCurrs==0)=[];



    %% validate all rois across all sessions - check they're active to final standards
    % this might take a hot second

    validActive = cell(numel(files2Load),1);
    for s = 1:numel(files2Load)

        % load session processed data results - get da raster
        selpath = files2Load{s};
        idcs = strfind(selpath,'\');
        pullPath = selpath(1:idcs(numel(idcs)-1)-1);
        dirData = dir(pullPath);
        tn = 0;
        for d = 1:numel(dirData)
            curName = dirData(d).name;
            if(contains(curName,'G8_'))
                tn=d;
            end
        end
        sessPath = strcat(pullPath,'\',dirData(tn).name,'\OverallResults\SessionMetrics.mat');
        sessData = load(sessPath);
        curUCurrs = sessData.uCurrs;
        curUChans = sessData.uChans;
        allRasterH = sessData.allRasterH;
        newRaster = false(size(allRasterH));

        curData = AllData{s};
        curAllTraces = curData.allTraces;
        regions = curData.regions;
        trialInfo = curData.TrialInfo;
        curRegSTDs = regionSTDs{s};
        curRegMeans = regionMeans{s};
        


        for ch = 1:numel(chans)
            for cu = 1:numel(allCurrs)
                m1 = trialInfo(:,1) == chans(ch);
                m2 = trialInfo(:,3) == allCurrs(cu);

                inds = find(m1 & m2);
                if(numel(inds)>0)
                    % load raster to find ROIs to check
                    cuI = find(curUCurrs == allCurrs(cu));
                    chI = find(curUChans == chans(ch));
                    curRaster = allRasterH(chI,cuI,:);
                    currRois = find(curRaster);
                    for r = 1:numel(currRois)
                    
                        midOnlyAve = NaN(numel(inds),10);
                        for i = 1:numel(inds)
                           midOnly = (curAllTraces{inds(i),2}-regionMeans{s})./regionMeans{s};
                           midOnlyAve(i,size(midOnly,2))=NaN;
                           midOnlyAve(i,1:size(midOnly,2))=smooth(midOnly(currRois(r),:),5);
                        end

                        % check if current trace has sufficently amplitude
                        % above threshold to be considered a spike.  This
                        % is to reduce the amount of high noise being
                        % processed.
                        curMidTrace = mean(midOnlyAve,1,'omitnan');
                        curSTD = curRegSTDs(currRois(r))./curRegMeans(currRois(r));
                        activeTP = curMidTrace > (2*curSTD);

                        % find largest continous region in mid trace
                        d = [true, diff(activeTP) ~= 0];   % TRUE if values change
                        b = activeTP(d);                   % Elements without repetitions
                        n = diff(find([d, true]));  % Number of repetitions
                        n = n(b == 1);   % Care for runs of 1's only (active)
                        m = max(n);
                        if(n > 5) % Valid if there is a peak not a single datapoint
                            sTime = find(activeTP,1,'first')*33;
                            if(sTime < 700) % start time must be before end of stimulation
                                newRaster(chI,cuI,currRois(r)) = 1;
                            end
                        end
                    end
                end
            end
        end
        validActive{s} = newRaster;
    end












































    
    
    %% Plot averaged volume of each stimulation scenario across all sessions longitudinally tracked
    % This will show the trial averaged fourescence intensity aligned and
    % with ROI's active during that session Localized.  This is as raw of
    % data as I can generate aside from the traces
    bmiFld = strcat(curFld,'BaseMeanImages\'); % creat folder in current directory
    mkdir(bmiFld)
    
    for f = 1:numel(files2Load) 
        curData = AllData{f};
        curRegions = curData.regions;
        sessionMean = (curData.meanMid + curData.meanPre + curData.meanPost)./3;

        selpath = files2Load{f};
        idcs = strfind(selpath,'\');
        pullPath = selpath(1:idcs(numel(idcs)-1)-1);
        dirData = dir(pullPath);
        
        % load up trial key and localize which trials belong to each
        % stimulation scenario.  Average those together
        fprintf('\nProcessing Dataset: ');
        disp(pullPath)

        keypathO = strcat(pullPath,'\KEY2D');

        % Check if Key folder exists, if it does, continue processing, if
        % not go to next folder for processing
        if ~exist(keypathO, 'dir')
            % Jump to next folder in source folder list
            continue
        end
        selpath = strcat(pullPath,'\RAW2D');
        segListing = dir(selpath);
        segListing(1:2) = []; % remove directory entries
         
        % load data 
        trialPath = strcat(pullPath,'\PlanarTrialsG8'); %construct folder for trial specific data
        destPath = strcat(trialPath,'\SessionMetrics.mat'); % trial key filename
        sessMetrics = load(destPath);
        
        
        % load processed data results
        tn = 0;
        for d = 1:numel(dirData)
            curName = dirData(d).name;
            if(contains(curName,'G8_'))
                tn=d;
            end
        end
        sessPath = strcat(pullPath,'\',dirData(tn).name,'\OverallResults\SessionMetrics.mat');
        sessData = load(sessPath);
        curUCurrs = sessData.uCurrs;
        curUChans = sessData.uChans;
        % allRasterH = sessData.allRasterH; % replaced with updated raster
        allRasterH = validActive{f};

        % Load metrics from session        
        dataCellArr = sessMetrics.dataCellArr;
        TrialInfo = sessMetrics.TrialInfo;
        validTrials = TrialInfo(:,6)~=0;
        
        % midStimFrames = NaN(512,512,100);
        % % Load trials into memory
        % for i = 1:100
        %     if(validTrials(i))
        %         curTrialData = load(strcat(trialPath,'\t',num2str(i),'.mat'),'midStim');
        %         midStimFrames(:,:,i) = curTrialData.midStim;
        %     end
        % end
        % 

        for ch = 1:numel(curUChans)
            for cu = 1:numel(curUCurrs)
            
                curChan = curUChans(ch);
                curCurr = curUCurrs(cu);

                if(curChan==target_chan)
                    if(curCurr==target_current)
    
                        a = TrialInfo(:,1) == curChan;
                        b = TrialInfo(:,3) == curCurr;
                        targetInds = find(a&b);
        
                        if(~isempty(targetInds))
                            aveVolOrig = NaN(512,512,numel(targetInds));
                            aveVol = NaN(512,512,numel(targetInds));

                            for t = 1:numel(targetInds)
                                curTrialData = load(strcat(trialPath,'\t',num2str(targetInds(t)),'.mat'),'midStim');
                                temp = curTrialData.midStim;

                                aveVol(:,:,t) =  (temp-sessionMean)./sessionMean;
                                aveVolOrig(:,:,t) =  temp;
                            end
                            aveFrame = mean(aveVol,3);
                            aveFrameOrig = mean(aveVolOrig,3);

                        end % if ~isempty(targetInds)
                        
                        % --- IN-NATIVE-LOOP AUTO EXPORT SCRIPT ---
                        % Generate ROIs safely knowing we have exactly the right native channel/curr indexes
                        curRaster = allRasterH(ch, cu, :);
                        curRInds = find(curRaster);

                        % Extract exact raw regions directly
                        curDataRaw = AllData{f};
                        curRegionsRaw = curDataRaw.regions;

                        outFld = 'S:\Roy\robin - Copy\figures\Fig3A\Weekly_AllActive_ROIs';
                        if ~exist(outFld, 'dir'), mkdir(outFld); end

                        baseName = sprintf('%s_AllActiveROI_Session%d_Chan%d_Curr%duA', icms, f, curChan, curCurr);

                        light_blue_vec = [76, 150, 220] / 255;

                        SE_close = strel('disk', 2);
                        SE_smooth = strel('disk', 1);

                        % Pre-compute all ROI boundaries
                        allBounds = {};
                        for r = 1:numel(curRInds)
                            curSesRoi = curRInds(r);
                            localBW = false(512, 512);
                            curPixels = curRegionsRaw(curSesRoi).PixelList;
                            for m = 1:size(curPixels, 1)
                                localBW(curPixels(m,2), curPixels(m,1)) = 1;
                            end
                            localBW = spatial_interp_patchwarp(localBW, ptForms{f}, 'euclidean', 1:512, 1:512);
                            localBW = localBW > 0.3;
                            localBW = imclose(localBW, SE_close);
                            localBW = imopen(localBW, SE_smooth);

                            boundaries = bwboundaries(localBW);
                            for k = 1:numel(boundaries)
                                b = boundaries{k};
                                if size(b,1) < 4, continue; end
                                bSmooth = [smoothdata(b(:,1), 'gaussian', 5), ...
                                           smoothdata(b(:,2), 'gaussian', 5)];
                                allBounds{end+1} = bSmooth;
                            end
                        end

                        % === 1) ROI-only SVG (no background) ===
                        curFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 512, 512]);
                        hold on; axis ij; axis equal;
                        xlim([1 512]); ylim([1 512]); axis off;
                        for k = 1:numel(allBounds)
                            bSmooth = allBounds{k};
                            patch(bSmooth(:,2), bSmooth(:,1), light_blue_vec, ...
                                'EdgeColor', 'none', 'FaceAlpha', 0.7);
                        end
                        saveas(curFig, fullfile(outFld, [baseName '.svg']));
                        close(curFig);

                        % === 2) Outlines overlaid on trial-averaged dF/F ===
                        if exist('aveFrame', 'var') && ~isempty(aveFrame)
                            % Warp dF/F to reference frame
                            aveFrameWarped = spatial_interp_patchwarp(aveFrame, ptForms{f}, 'euclidean', 1:512, 1:512);

                            % Normalize for display
                            dispImg = aveFrameWarped;
                            pLo = prctile(dispImg(:), 1);
                            pHi = prctile(dispImg(:), 99);
                            dispImg = (dispImg - pLo) / (pHi - pLo);
                            dispImg = max(0, min(1, dispImg));

                            curFig2 = figure('Visible', 'off', 'Color', 'k', 'Position', [100, 100, 512, 512]);
                            imagesc(dispImg); axis ij; axis equal; axis off;
                            colormap('gray'); hold on;

                            for k = 1:numel(allBounds)
                                bSmooth = allBounds{k};
                                plot(bSmooth(:,2), bSmooth(:,1), '-', ...
                                    'Color', light_blue_vec, 'LineWidth', 1.2);
                            end

                            overlayName = sprintf('%s_Overlay_Session%d_Chan%d_Curr%duA', icms, f, curChan, curCurr);
                            saveas(curFig2, fullfile(outFld, [overlayName '.png']));
                            saveas(curFig2, fullfile(outFld, [overlayName '.svg']));
                            close(curFig2);
                        end
                        
                    end % if curCurr==target_current
                end % if curChan==target_chan
            end % for cu
        end % for ch
        
    end % for f = 1:numel...

fprintf('\nSuccess %s! Exported beautifully to: %s\n', icms, outFld);

end % <-- CLOSES THE MASTER 'for animal = 1:3' LOOP
