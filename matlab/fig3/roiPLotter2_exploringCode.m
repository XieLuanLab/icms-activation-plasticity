% Requires raw session data (not in repo).
% Internal paths must be configured manually for your local environment.

% calculate trace means and stds
% TODO - report baseline trace means over each session, i.e. see if there is a drifting baseline in each session or for each roi
figGen=0;
    fMaxAll = [];
bMaxAll = [];
% allDFF = cell(3,1);
% allOnset = cell(3,1);
% allDur = cell(3,1);
% allSpikeNum = cell(3,1);
% allSpikeOnset = cell(3,1);
% allVol = cell(3,1);
% allSubset = cell(3,1);
% allDays = cell(3,1);
% distanceSub = cell(3,1);
% distanceNon = cell(3,1);
% distanceSub_ver2 = cell(3,1);
% distanceNon_ver2 = cell(3,1);
% 
% endDFF = cell(3,1);
% endOnset = cell(3,1);
% endDur = cell(3,1);
% endSpikeNum = cell(3,1);
% endSpikeOnset = cell(3,1);
% endVol = cell(3,1);





for animal =  1:3
    if(animal==1)
    icms = 'ICMS92';
%     files2Load = {'Z:\ICMS92\9-6-23\PlanarTrials\SessionMetrics.mat';...
%                 'Z:\ICMS92\9-8-23\PlanarTrials\SessionMetrics.mat';...
%                 'Z:\ICMS92\9-12-23\PlanarTrials\SessionMetrics.mat';...
%                 'Z:\ICMS92\9-14-23\PlanarTrials\SessionMetrics.mat';...
%                 'Z:\ICMS92\9-19-23\PlanarTrials\SessionMetrics.mat';...
%                 'Z:\ICMS92\9-21-23\PlanarTrials\SessionMetrics.mat';...
%                 'Z:\ICMS92\9-25-23\PlanarTrials\SessionMetrics.mat'};
    files2Load = {'Z:\xl_stimulation\ICMS92\9-6-23\PlanarTrialsG8\SessionMetrics.mat';...
                'Z:\xl_stimulation\ICMS92\9-8-23\PlanarTrialsG8\SessionMetrics.mat';...
                'Z:\xl_stimulation\ICMS92\9-12-23\PlanarTrialsG8\SessionMetrics.mat';...
                'Z:\xl_stimulation\ICMS92\9-14-23\PlanarTrialsG8\SessionMetrics.mat';...
                'Z:\xl_stimulation\ICMS92\9-19-23\PlanarTrialsG8\SessionMetrics.mat';...
                'Z:\xl_stimulation\ICMS92\9-21-23\PlanarTrialsG8\SessionMetrics.mat';...
                'Z:\xl_stimulation\ICMS92\9-25-23\PlanarTrialsG8\SessionMetrics.mat'};
            
    minSetCnt = 4;
   % minSetCnt = 0;

    refInd = 3;
    chans = [9,11,12]; % hardcoded for this animal - update for each animal
    daysTrained = [7,9,13,15,20,22,26];
    
    elseif(animal==2)
        icms = 'ICMS98';
%         files2Load = {'Z:\ICMS98\10-20-23\PlanarTrials\SessionMetrics.mat';...
%                       'Z:\ICMS98\10-24-23\PlanarTrials\SessionMetrics.mat';...
%                       'Z:\ICMS98\10-26-23\PlanarTrials\SessionMetrics.mat';...
%                       'Z:\ICMS98\10-31-23\PlanarTrials\SessionMetrics.mat';...
%                       'Z:\ICMS98\11-2-23\PlanarTrials\SessionMetrics.mat';...
%                       'Z:\ICMS98\11-7-23\PlanarTrials\SessionMetrics.mat';...
%                       'Z:\ICMS98\11-20-23\PlanarTrials\SessionMetrics.mat';...
%                       'Z:\ICMS98\11-22-23\PlanarTrials\SessionMetrics.mat'};
        files2Load = {'Z:\xl_stimulation\ICMS98\10-20-23\PlanarTrialsG8\SessionMetrics.mat';...
                      'Z:\xl_stimulation\ICMS98\10-24-23\PlanarTrialsG8\SessionMetrics.mat';...
                      'Z:\xl_stimulation\ICMS98\10-26-23\PlanarTrialsG8\SessionMetrics.mat';...
                      'Z:\xl_stimulation\ICMS98\10-31-23\PlanarTrialsG8\SessionMetrics.mat';...
                      'Z:\xl_stimulation\ICMS98\11-2-23\PlanarTrialsG8\SessionMetrics.mat';...
                      'Z:\xl_stimulation\ICMS98\11-7-23\PlanarTrialsG8\SessionMetrics.mat';...
                      'Z:\xl_stimulation\ICMS98\11-17-23\PlanarTrialsG8\SessionMetrics.mat'};
        refInd = 5;
        minSetCnt = 3;
                % minSetCnt = 0;

        chans = [4,5,6]; % hardcoded for this animal - update for each animal
        daysTrained = [0,4,6,11,13,18,28];


    elseif(animal==3)
        icms = 'ICMS93';
%         files2Load = {'Z:\ICMS93\8-30-23\PlanarTrials\SessionMetrics.mat';...
%                       'Z:\ICMS93\9-6-23\PlanarTrials\SessionMetrics.mat';...
%                       'Z:\ICMS93\9-14-23\PlanarTrials\SessionMetrics.mat';...
%                       'Z:\ICMS93\9-22-23\PlanarTrials\SessionMetrics.mat';...
%                       'Z:\ICMS93\9-26-23\PlanarTrials\SessionMetrics.mat'};
        files2Load = {'Z:\xl_stimulation\ICMS93\8-30-23\PlanarTrialsG8\SessionMetrics.mat';...
                      'Z:\xl_stimulation\ICMS93\9-6-23\PlanarTrialsG8\SessionMetrics.mat';...
                      'Z:\xl_stimulation\ICMS93\9-14-23\PlanarTrialsG8\SessionMetrics.mat';...
                      'Z:\xl_stimulation\ICMS93\9-22-23\PlanarTrialsG8\SessionMetrics.mat';...
                      'Z:\xl_stimulation\ICMS93\9-26-23\PlanarTrialsG8\SessionMetrics.mat'};
        refInd = 2;
        minSetCnt = 3;
                % minSetCnt = 0;

        chans = [9,11,12]; % hardcoded for this animal - update for each animal
        daysTrained = [0,7,15,23,27];     
        
%     else
%         icms = 'ICMS101';
%         files2Load = {'Z:\\ICMS101Data\11-1-23\PlanarTrials\SessionMetrics.mat';...
%                         'Z:\\ICMS101Data\11-16-23\PlanarTrials\SessionMetrics.mat';...
%                         'Z:\\ICMS101Data\11-21-23\PlanarTrials\SessionMetrics.mat';...
%                         'Z:\\ICMS101Data\11-22-23\PlanarTrials\SessionMetrics.mat';...
%                         'Z:\\ICMS101Data\11-29-23\PlanarTrials\SessionMetrics.mat'};
%         refInd = 3;
%         minSetCnt = 3;
%         chans = [10,11,12]; % hardcoded for this animal - update for each animal
%         daysTrained = [6,21,26,27,34];   
    end

    % Create folder
    curFld = strcat(icms,'\'); % creat folder in current directory
    mkdir(curFld)
    srcFld = strcat(curFld,'Source\'); % creat folder in current directory
    mkdir(srcFld)


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
            % run multiple iteration of the fitting and accept the one with the
            % smallest mis-alignment score
            for t = 1:30

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

                if(t==1)
                    ptForms{f} = patchwarp_results.warp1_cell{1};
                    alignVal = curAV;
                else
                    % perform optimization check
                    if(alignVal>curAV)
                        alignVal = curAV;
                        ptForms{f} = patchwarp_results.warp1_cell{1};
                    end
                end
            end        
            if(figGen==1)
                curFig = figure('NumberTitle', 'off', 'Name', strcat('Sess',num2str(refInd),'--vs--Sess',num2str(f)));
                shiftedplane = spatial_interp_patchwarp(image2_mean, ptForms{f}, 'euclidean', 1:512, 1:512);
                imshow(imfuse(shiftedplane./max(shiftedplane,[],'all'),image1_mean./max(image1_mean,[],'all'),'falsecolor','Scaling','joint','ColorChannels','red-cyan'));
                txt = strcat('Sess',num2str(refInd),'vsSess',num2str(f));
                combineFileName = strcat(curFld,'Alignment-',txt,'.fig');
                % saveas(curFig,combineFileName);  
                combineFileName = strcat(curFld,'Alignment-',txt,'.tif');
                % saveas(curFig,combineFileName);
                 close(curFig)
                
                curFig = figure('NumberTitle', 'off', 'Name', strcat('Sess',num2str(refInd),'--vs--Sess',num2str(f)));
                shiftedplane = spatial_interp_patchwarp(image2_mean, ptForms{f}, 'euclidean', 1:512, 1:512);
                imshowpair(shiftedplane./max(shiftedplane,[],'all'),image1_mean./max(image1_mean,[],'all'),'montage');
                txt = strcat('Sess',num2str(refInd),'vsSess',num2str(f));
                combineFileName = strcat(curFld,'AlignmentM-',txt,'.fig');
                % saveas(curFig,combineFileName);  
                combineFileName = strcat(curFld,'AlignmentM-',txt,'.tif');
                % saveas(curFig,combineFileName);
                 close(curFig)
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
%     end
%     refInd
    
%     %% Align electrode positions across sessions
%     chanPositionsAligned = chanPositions;
%     figure
%     hold on
%     for f = 1:numel(files2Load)             
%         for ch = 1:numel(chans)
%             if(~isnan(chanPositions(ch,1,f)))
%                 zProjectedBW = zeros(512,512);
%                 zProjectedBW(round(chanPositions(ch,1,f)),round(chanPositions(ch,2,f)))=1;
%                 shiftedplane = spatial_interp_patchwarp(zProjectedBW, ptForms{f}, 'euclidean', 1:512, 1:512);
%                 [~,chanPositionsAligned(ch,1,f)] = max(max(shiftedplane,[],2));
%                 [~,chanPositionsAligned(ch,2,f)] = max(max(shiftedplane,[],1));
%             end
%         end
%         scatter(chanPositionsAligned(:,1,f),chanPositionsAligned(:,2,f))
%         ylim([1 512])
%         xlim([1 512])
%     end
%     legend()
    
    % Average aligned electrode positions across sessions for a agreed
    % reference
    

    

    %% create updated region positions for each dataset
    updatedRegions = cell(numel(files2Load),1);
    processedRegions =  cell(numel(files2Load),1);
    shiftedSum = zeros(512,512);
    shiftedSumV2 = zeros(512,512);
    for f = 1:numel(files2Load) 
        curData = load(files2Load{f});
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
        combineFileName = strcat(curFld,'OverlapMap.fig');
        % saveas(curFig,combineFileName);  
        combineFileName = strcat(curFld,'OverlapMap.tif');
        % saveas(curFig,combineFileName);
         close(curFig)
    end

    
    updatedRegionsEnd = updatedRegions;
    curprocessedRegionsProcEnd = processedRegions;
    
    
    
    
    


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
                if(sum(setMask) > minSetCnt)
                    rsC = rsC+1;
                    roiSets{rsC} = curSet;
                    roiSetsPx{rsC} = curRegions{r1};
                end
            end
        end
    end

    
        
    %% Plot aligned datasets and selected ROIs 
    if(figGen==3)
        for f = 1:numel(files2Load) 
            curData = load(files2Load{f});
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
            combineFileName = strcat(curFld,'AlignmentROIs-',txt,'.tif');
            saveas(curFig,combineFileName);
            close(curFig)
        end
    end
    
    % Plot 
    if(figGen==3)
        for f = 1:numel(files2Load) 
            curData = load(files2Load{f});
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
                    curData = load(files2Load{fInd});
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
            combineFileName = strcat(curFld,'AlignmentTransformedROIs-',txt,'.fig');
            % saveas(curFig,combineFileName);  
            combineFileName = strcat(curFld,'AlignmentTransformedROIs-',txt,'.tif');
            % saveas(curFig,combineFileName);
            close(curFig)
        end
    end
    
    
    
    
    
%     %% Plot aligned datasets and activated ROIs for different Stimulation scenarios
%     for f = 1:numel(files2Load) 
%         curData = load(files2Load{f});
%         curRegions = curData.regions;
%         image2_mean = curData.meanPre./max(curData.meanPre,[],'all'); 
%         
%         image2_mean = spatial_interp_patchwarp(image2_mean, ptForms{f}, 'euclidean', 1:512, 1:512);
%         
%         selpath = files2Load{f};
%         idcs = strfind(selpath,'\');
%         pullPath = selpath(1:idcs(numel(idcs)-1)-1);
%         dirData = dir(pullPath);
%         tn = 0;
%         for d = 1:numel(dirData)
%             curName = dirData(d).name;
%             if(contains(curName,'G7'))
%                 tn=d;
%             end
%         end
%         sessPath = strcat(pullPath,'\',dirData(tn).name,'\OverallResults\SessionMetrics.mat');
%         sessData = load(sessPath);
%         
%         curUCurrs = sessData.uCurrs;
%         curUChans = sessData.uChans;
%         allRasterH = sessData.allRasterH;
%         
%         for c1 = 1:numel(curUChans)
%             for c2 = 1:numel(curUCurrs)
%             curRaster = allRasterH(c1,c2,:);
%             
%             curFig = figure('NumberTitle', 'off', 'Name', strcat('Sess',num2str(f)),'WindowState','maximized');
%             imshow(image2_mean./max(image2_mean,[],'all'))
%             hold on
%             
%             curRInds = find(curRaster);
%             if(~isempty(curSessROIind>0))
%                 for r = 1:numel(curRInds)
%                     zProjectedBW = false(512,512);
% 
%                     curPixels = curRegions(curRInds(r)).PixelList;
%                     for m = 1:size(curPixels,1)
%                         zProjectedBW(curPixels(m,2),curPixels(m,1))=1;
%                     end
%                     zProjectedBW = spatial_interp_patchwarp(zProjectedBW, ptForms{f}, 'euclidean', 1:512, 1:512);
%                     zProjectedBW = zProjectedBW>0;
%                     
%                     
%                     boundaries = bwboundaries(zProjectedBW);
%                     for k=1:numel(boundaries)
%                         b = boundaries{k};
%                         plot(b(:,2),b(:,1),'r','LineWidth',1);
%                     end
%                 end
%             end
% 
%             txt = strcat('Sess',num2str(f),'Chan',num2str(curUChans(c1)),'Curr',num2str(curUCurrs(c2)));
%             combineFileName = strcat(curFld,'AlignmentROIs-',txt,'.tif');
%             saveas(curFig,combineFileName);
%             close(curFig)
%             end
%         end
%     end
    
    
    
    
    
     %% Plot overlapping aligned datasets and activated ROIs for different Stimulation scenarios
    
% 
%     for c1 = 1:3
%         for c2 = 3:6
%             curFig = figure('WindowState','maximized');
%             shiftedSum = zeros(512,512);
%             shiftedVol = false(512,512,numel(files2Load));
%             for f = 1:numel(files2Load) 
%                 curData = load(files2Load{f});
%                 curRegions = curData.regions;
%                 
%                 selpath = files2Load{f};
%                 idcs = strfind(selpath,'\');
%                 pullPath = selpath(1:idcs(numel(idcs)-1)-1);
%                 dirData = dir(pullPath);
%                 tn = 0;
%                 for d = 1:numel(dirData)
%                     curName = dirData(d).name;
%                     if(contains(curName,'G7'))
%                         tn=d;
%                     end
%                 end
%                 sessPath = strcat(pullPath,'\',dirData(tn).name,'\OverallResults\SessionMetrics.mat');
%                 sessData = load(sessPath);
%                 curUCurrs = sessData.uCurrs;
%                 curUChans = sessData.uChans;
%                 allRasterH = sessData.allRasterH;
%                 
%                 cInd = find(curUCurrs==c2);
%                 if(~isempty(cInd))
%                     curRaster = allRasterH(c1,cInd,:);
%                     curRInds = find(curRaster);
%                     if(~isempty(curSessROIind>0))
%                         curSliceROIs = false(512,512);
%                         for r = 1:numel(curRInds)
%                             zProjectedBW = false(512,512);
% 
%                             curPixels = curRegions(curRInds(r)).PixelList;
%                             for m = 1:size(curPixels,1)
%                                 zProjectedBW(curPixels(m,2),curPixels(m,1))=1;
%                             end
%                             zProjectedBW = spatial_interp_patchwarp(zProjectedBW, ptForms{f}, 'euclidean', 1:512, 1:512);
%                             zProjectedBW = zProjectedBW>0;
% 
%                             if(sum(zProjectedBW,'all')>1)
%                                 regionsV2{r} = find(shiftedplane);
%                                 shiftedSum = shiftedSum + zProjectedBW;
%                                 curSliceROIs = curSliceROIs + zProjectedBW;
%                             end
% 
%     %                         boundaries = bwboundaries(zProjectedBW);
%     %                         for k=1:numel(boundaries)
%     %                             b = boundaries{k};
%     %                             plot(b(:,2),b(:,1),'r','LineWidth',1);
%     %                         end
%                         end
%                         shiftedVol(:,:,f) = curSliceROIs;
%                     end
%                 end
%             end
%             imshow(shiftedSum./numel(files2Load))
%             txt = strcat('OverlappingROIsV2  Chan',num2str(curUChans(c1)),'Curr',num2str(c2));
%             combineFileName = strcat(curFld,'AlignmentROIs-',txt,'.tif');
%             saveas(curFig,combineFileName);
%             close(curFig)
% % %             
% % %             txt = strcat('OverlappingROIs  Chan',num2str(curUChans(c1)),'Curr',num2str(c2));
% % %             figure('NumberTitle', 'off', 'Name', txt,'WindowState','maximized');
% % %             imshow3D(shiftedVol)
%         end
%     end
%     
%     
%     
%     
%     
    
%     
%     
%     %%
%     for c1 = 1:3
%         for c2 = 3:6
%             
%             for f = 1:numel(files2Load) 
%                 curData = load(files2Load{f});
%                 curRegions = curData.regions;
%                 
%                 selpath = files2Load{f};
%                 idcs = strfind(selpath,'\');
%                 pullPath = selpath(1:idcs(numel(idcs)-1)-1);
%                 dirData = dir(pullPath);
%                 tn = 0;
%                 for d = 1:numel(dirData)
%                     curName = dirData(d).name;
%                     if(contains(curName,'G7'))
%                         tn=d;
%                     end
%                 end
%                 sessPath = strcat(pullPath,'\',dirData(tn).name,'\OverallResults\SessionMetrics.mat');
%                 sessData = load(sessPath);
%                 curUCurrs = sessData.uCurrs;
%                 curUChans = sessData.uChans;
%                 allRasterH = sessData.allRasterH;
%                 
%                 cInd = find(curUCurrs==c2);
%                 if(~isempty(cInd))
%                     curRaster = allRasterH(c1,cInd,:);
%                     curRInds = find(curRaster);
%                     if(~isempty(curSessROIind>0))
%                         curSliceROIs = false(512,512);
%                         for r = 1:numel(curRInds)
%                             zProjectedBW = false(512,512);
% 
%                             curPixels = curRegions(curRInds(r)).PixelList;
%                             for m = 1:size(curPixels,1)
%                                 zProjectedBW(curPixels(m,2),curPixels(m,1))=1;
%                             end
%                             zProjectedBW = spatial_interp_patchwarp(zProjectedBW, ptForms{f}, 'euclidean', 1:512, 1:512);
%                             zProjectedBW = zProjectedBW>0;
% 
%                             if(sum(zProjectedBW,'all')>1)
%                                 curSliceROIs = curSliceROIs + zProjectedBW;
%                             end
%                         end
%                         curFig = figure('WindowState','maximized');
%                         imshow(curSliceROIs)
%                         txt = strcat('OverlappingROIsV2  Chan',num2str(curUChans(c1)),'Curr',num2str(c2),'Sess',num2str(f));
%                         combineFileName = strcat(curFld,'AlignmentROIs-',txt,'.tif');
%                         saveas(curFig,combineFileName);
%                         close(curFig)
%                     end
%                 end
%             end
%         end
%     end
%     
%     
    
    %% Find ROIs uniquely activated in the last two sessions (i.e. never before activated)
    roiEndSets = cell(1,1);
    roiEndSetsPx = cell(1,1);
    rsC = 0;
    for f = 1:numel(files2Load)
        curRegions = updatedRegionsEnd{f};
        curProc = processedRegionsEnd{f};
        for r1 = 1:numel(curRegions) % iterate through all regions and find subsequent regions that overlap
            if(curProc(r1)==0) % only process valid regions that have not been processed yet
                % initialze current set
                curSet = zeros(numel(files2Load),1);
                curSet(f) = r1;

                for f2 = f+1:numel(files2Load)
                    searchRegions = updatedRegionsEnd{f2};
                    searchProc = processedRegionsEnd{f2};
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
                        processedRegionsEnd{f2} = searchProc;
                    end
                end

                % add current set of ROIs to master list if active only in
                % last two sessions for the ROI
                
                if(ceil(daysTrained(f)./7)==4)
                % if(f >= (numel(files2Load)-1))
                    setMask = curSet>0;
                    if(sum(setMask) == 1) % must be uniquely active
                        rsC = rsC+1;
                        roiEndSets{rsC} = curSet;
                        roiEndSetsPx{rsC} = curRegions{r1};
                    end
                end


            end
        end
    end

    
    %% Display selected ROIs
    if(figGen==1)
        curFig = figure('WindowState','maximized');
        imshow(refMean./max(refMean,[],'all'))
        hold on
        for r = 1:numel(roiSets)
            zProjectedBW = false(512,512);
            curPx = roiSetsPx{r};
            for m = 1:size(curPx,1)
                zProjectedBW(curPx(m))=1;
            end
            boundaries = bwboundaries(zProjectedBW);
            for k=1:numel(boundaries)
                b = boundaries{k};
                plot(b(:,2),b(:,1),'r','LineWidth',1);
            end
        end
        combineFileName = strcat(curFld,'SelectedROIs.fig');
        % saveas(curFig,combineFileName);  
        combineFileName = strcat(curFld,'SelectedROIs.tif');
        % saveas(curFig,combineFileName);
        close(curFig)

        curFig = figure('WindowState','maximized');
        imshow(refMean./max(refMean,[],'all'))
        hold on
        for r = 1:numel(roiSets)
            zProjectedBW = false(512,512);
            curPx = roiSetsPx{r};
            for m = 1:size(curPx,1)
                zProjectedBW(curPx(m))=1;
            end
            boundaries = bwboundaries(zProjectedBW);
            for k=1:numel(boundaries)
                b = boundaries{k};
                plot(b(:,2),b(:,1),'r','LineWidth',1);
            end
            [rx,ry] = ind2sub([512,512],curPx(1));
    %         text(ry,rx,num2str(r),'Color','y')
        end
        combineFileName = strcat(curFld,'SelectedROIsNumbered.fig');
        % saveas(curFig,combineFileName);  
        combineFileName = strcat(curFld,'SelectedROIsNumbered.tif');
        % saveas(curFig,combineFileName);
        close(curFig)
    
    end
    
    
    

    
    
    
    
%     %% Plot overlapping selected ROIs 
%     for r = 1:numel(roiSets)
%         curFig = figure();
% %         imshow(zeros(512,512))
%         hold on
%         curROIInds = roiSets{r};
%         for f = 1:numel(files2Load) 
%             curData = load(files2Load{f});
%             curRegions = curData.regions;
%             curSessROIind = curROIInds(f);
%             if(curSessROIind>0)
%                 zProjectedBW = false(512,512);
%                 
%                 curPixels = curRegions(curSessROIind).PixelList;
%                 for m = 1:size(curPixels,1)
%                     zProjectedBW(curPixels(m,2),curPixels(m,1))=1;
%                 end
%                 shiftedplane = spatial_interp_patchwarp(zProjectedBW, ptForms{f}, 'euclidean', 1:512, 1:512);
%                 shiftedplane = shiftedplane>0;
%                 
%                 boundaries = bwboundaries(shiftedplane);
%                 for k=1:numel(boundaries)
%                     b = boundaries{k};
%                     plot(b(:,2),b(:,1));
%                 end
%             end
%         end
%         curROIInds(curROIInds==0)=[];
%         legend()
% 
%         title(strcat('ROI Overlap',num2str(r)));
%        
% %         combineFileName = strcat(curFld,'OverlaptROIs-',num2str(r),'.fig');
% %         saveas(curFig,combineFileName);  
%         combineFileName = strcat(curFld,'OverlaptROIs-',num2str(r),'.tif');
%         saveas(curFig,combineFileName);
%         close(curFig)
%     end
    
    
    AllData = cell(numel(files2Load),1);
    for s = 1:numel(files2Load)
        AllData{s} = load(files2Load{s});
    end


    %% Generate traces for all ROIs selected
    ROIs2Compare = zeros(numel(roiSets),numel(files2Load));
    for r = 1:numel(roiSets) % convert keys to single matrix
        ROIs2Compare(r,:) = roiSets{r};
    end

    % calculate trace means and stds
    allCurrs = [];
    regionMeans = cell(numel(files2Load),1);
    regionSTDs = cell(numel(files2Load),1);
    regionPreStimMeans = cell(numel(files2Load),1);
    regionPreStimSTDs = cell(numel(files2Load),1);
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
        allSpikes = cell(size(curAllTraces));
        completeTrace = zeros(numel(regions),1);
        completePreTrace = zeros(numel(regions),1);
        for t = 1:100
            if(trialInfo(t,6)~=0)
%                 if(trialInfo(t,4)==1)
                    temp = [curAllTraces{t,1},curAllTraces{t,2},curAllTraces{t,3}];
                    completeTrace(:,cN+1:cN+size(temp,2)) = temp;
                    completePreTrace(:,cPN+1:cPN+size(curAllTraces{t,1},2)) = curAllTraces{t,1};
                    cN = cN + size(temp,2);
                    cPN = cPN + size(curAllTraces{t,1},2);
%                 end
            end
        end
        if(cN>0)
            % for r = 1:size(completeTrace,1)
            %     completeTrace(r,:) = smooth(completeTrace(r,:),5);
            % end

            regionMeans{s} = mean(completeTrace,2);
            regionSTDs{s} = std(completeTrace,[],2);
            
            regionPreStimMeans{s} = mean(completePreTrace,2);
            regionPreStimSTDs{s} = std(completePreTrace,[],2);
            
            % Perform deconvolution
            cleanTrace = zeros(size(completeTrace));
            spikes = zeros(size(completeTrace));
            for r = 1:size(completeTrace,1)
                curTrace = smooth(completeTrace(r,:),5);
                % curTrace = completeTrace(r,:);
                curTrace = (curTrace-mean(curTrace))./mean(curTrace);
                g = 0.95;
                [cleanTrace(r,:),temp , options] = deconvolveCa(curTrace, 'ar1', g,'thresholded', 'optimize_smin', 'optimize_pars', 'thresh_factor', 0.95);
                spikes(r,:) = temp > 0.1;
            end
            
            % break deconvoluted data down to individual trials
            cNR = 0;
            for t = 1:100
                if(trialInfo(t,6)~=0)
%                     if(trialInfo(t,4)==1)
                    % measure original segment size
                    s1 = size(curAllTraces{t,1},2);
                    s2 = size(curAllTraces{t,2},2);
                    s3 = size(curAllTraces{t,3},2);
                    
                    % pull those segments from deconvoluted data and store
                    % in structure matching original trace data
                    cleanAllTraces{t,1} = cleanTrace(:,cNR+1:cNR+s1);
                    cleanAllTraces{t,2} = cleanTrace(:,cNR+s1:cNR+s1+s2);
                    cleanAllTraces{t,3} = cleanTrace(:,cNR+s1+s2:cNR+s1+s2+s3);
                    
                    allSpikes{t,1} = spikes(:,cNR+1:cNR+s1);
                    allSpikes{t,2} = spikes(:,cNR+s1:cNR+s1+s2);
                    allSpikes{t,3} = spikes(:,cNR+s1+s2:cNR+s1+s2+s3);
                    
                    cNR = cNR + s1+s2+s3;
%                     end
                end
            end
        end
        
        curData.cleanTraces = cleanAllTraces;
        curData.deconvolute = allSpikes;
        AllData{s} = curData;
    end
    allCurrs = unique(allCurrs);
    allCurrs(allCurrs==0)=[];


    midTraces = cell(size(ROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    preTraces = cell(size(ROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeTraces = cell(size(ROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    preSpikeTraces = cell(size(ROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    midSTD = cell(size(ROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeOnset = zeros(size(ROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeVol = zeros(size(ROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeDur = zeros(size(ROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeMax = zeros(size(ROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeN20 = zeros(size(ROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeCount = zeros(size(ROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeMaxTime = zeros(size(ROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeInitialTime = zeros(size(ROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeDeconDur = zeros(size(ROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikePeaks = zeros(size(ROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    for curR = 1:size(ROIs2Compare,1)
        currRois = ROIs2Compare(curR,:);
        for s = 1:numel(files2Load)
            if(currRois(s)~=0)
                curData = AllData{s};
                curAllTraces = curData.allTraces;
                curAllTracesClean = curData.cleanTraces;
                regions = curData.regions;
                trialInfo = curData.TrialInfo;
                spikes = curData.deconvolute;
                curRegSTDs = regionSTDs{s};
                curRegMeans = regionMeans{s};
                

                for ch = 1:numel(chans)
                    for cu = 1:numel(allCurrs)
                        m1 = trialInfo(:,1) == chans(ch);
                        m2 = trialInfo(:,3) == allCurrs(cu);
                        inds = find(m1 & m2);
                        if(numel(inds)>0)
                            preAve = NaN(numel(inds),10);
                            midPostAve = NaN(numel(inds),10);
                            midOnlyAve = NaN(numel(inds),10);
                            preAveSpikes = NaN(numel(inds),10);
                            midPostAveSpikes = NaN(numel(inds),10);
                            midOnlyAveSpikes = NaN(numel(inds),10);
                            for i = 1:numel(inds)
                               pre = (curAllTraces{inds(i),1}-regionMeans{s})./regionMeans{s};
                               midPost = ([curAllTraces{inds(i),2},curAllTraces{inds(i),3}]-regionMeans{s})./regionMeans{s};

                               preAve(i,size(pre,2))=NaN;
                               preAve(i,1:size(pre,2))=smooth(flip(pre(currRois(s),:)),5);
                               midPostAve(i,size(midPost,2))=NaN;
                               midPostAve(i,1:size(midPost,2))=smooth(midPost(currRois(s),:),5);
                               midOnly = (curAllTraces{inds(i),2}-regionMeans{s})./regionMeans{s};
                               midOnlyAve(i,size(midOnly,2))=NaN;
                               midOnlyAve(i,1:size(midOnly,2))=smooth(midOnly(currRois(s),:),5);
                               
                               
                               preSpikes = spikes{inds(i),1};
                               midOnlySpikes = spikes{inds(i),2};
                               midPostSpikes = [spikes{inds(i),2},spikes{inds(i),3}];
                               
                               preAveSpikes(i,size(preSpikes,2))=NaN;
                               preAveSpikes(i,1:size(preSpikes,2))=flip(preSpikes(currRois(s),:));
                               midPostAveSpikes(i,size(midPostSpikes,2))=NaN;
                               midPostAveSpikes(i,1:size(midPostSpikes,2))=midPostSpikes(currRois(s),:);
                               midOnlyAveSpikes(i,size(midOnlySpikes,2))=NaN;
                               midOnlyAveSpikes(i,1:size(midOnlySpikes,2))=midOnlySpikes(currRois(s),:);
                            end

                            % check if current trace has sufficently amplitude
                            % above threshold to be considered a spike.  This
                            % is to reduce the amount of high noise being
                            % processed.
                            curMidTrace = mean(midOnlyAve,1,'omitnan');
                            curSTD = curRegSTDs(currRois(s))./curRegMeans(currRois(s));
                            activeTP = curMidTrace > (2*curSTD);
                            activeTP2 = curMidTrace > (curSTD);

                            % find largest continous region in mid trace
                            d = [true, diff(activeTP) ~= 0];   % TRUE if values change
                            b = activeTP(d);                   % Elements without repetitions
                            n = diff(find([d, true]));  % Number of repetitions
                            n = n(b == 1);   % Care for runs of 1's only (active)
                            m = max(n);
                            if(n > 5) % Valid if there is a peak not a single datapoint
                                sTime = find(activeTP,1,'first')*33;
                                if(sTime < 700) % start time must be before end of stimulation


                                    midPostTrace = mean(midPostAve,1,'omitnan');
                                    preTrace = mean(preAve,1,'omitnan');




                                    preTraces{curR,s,ch,cu} = flip(preTrace);
                                    midTraces{curR,s,ch,cu} = midPostTrace; % mid and post trace
                                    spikeTraces{curR,s,ch,cu} = mean(midPostAveSpikes,1,'omitnan');

                                    % find initial spike time from Deconvolved data
                                    iSpikes = zeros(size(midOnlyAveSpikes,1),1);
                                    eSpikes = zeros(size(midOnlyAveSpikes,1),1);
                                    for g = 1:size(midOnlyAveSpikes,1)
                                        if(isempty(find(midOnlyAveSpikes(g,:),1)))
                                            iSpikes(g) = NaN;
                                            eSpikes(g) = NaN;
                                        else
                                            iSpikes(g) = find(midOnlyAveSpikes(g,:),1);
                                            eSpikes(g) = find(midOnlyAveSpikes(g,:),1,'last');
                                        end
                                    end
                                    spikeInitialTime(curR,s,ch,cu) = mean(iSpikes,'omitnan')*33;


                                    spikeDeconDur(curR,s,ch,cu) =  mean(eSpikes - iSpikes,'omitnan')*33;
                                    temp = smooth( mean(midPostAveSpikes,1,'omitnan'),5);
                                    spikePeaks(curR,s,ch,cu) = numel(findpeaks(temp,'MinPeakProminence',0.05));

                                    midSpikeTraces = mean(midOnlyAveSpikes,1,'omitnan');
                                    preSpikeTraces{curR,s,ch,cu} = flip(mean(preAveSpikes,1,'omitnan'));
                                    % Determine onset and end times for current
                                    % trace
                                    % activeTP2 = midPostTrace > curSTD;
                                    % activeTP2 = activeTP;


                                    % save metrics of trace for analysis
                                    spikeOnset(curR,s,ch,cu) = find(activeTP2,1,'first')*33; 
                                    spikeDur(curR,s,ch,cu) = (find(activeTP2,1,'last') - find(activeTP2,1,'first'))*33;
                                    spikeVol(curR,s,ch,cu) = sum(midPostTrace(find(activeTP2,1,'first') : find(activeTP2,1,'last')));
                                    spikeMax(curR,s,ch,cu) = max(midPostTrace);
                                    spikeCount(curR,s,ch,cu) = sum(midSpikeTraces); % number of spikes during stimulation period


                                    % find max spike time
                                    [~,mInd] = max(midPostTrace);
                                    spikeMaxTime(curR,s,ch,cu) = mInd*33;
                                    midPostTraceN = midPostTrace./max(midPostTrace);
                                    iData = interp1(1:numel(midPostTraceN),midPostTraceN,1:0.1:numel(midPostTraceN)/2);
                                    iX = 1:0.1:numel(midPostTraceN)/2;
                                    iDataR = round(iData*10)./10;
                                    tInd = find(iDataR==0.2,1);
                                    if(~isempty(tInd))
                                        spikeN20(curR,s,ch,cu) = iX(tInd)*33;
                                    end


    %                                 % plot individual trace data and subsequent averages
    %                                 curFig = figure();
    %                                 dff = [preAve,midPostAve];
    % %                                 x = 33*((1:numel(dff))-numel(pre));
    %                                 h = heatmap(dff);
    %                                 h.GridVisible = 'off';
    %                                 colormap(parula)
    %                                 xlim([-1000 4000])
    %                                 xline(0,'--');
    %                                 set(gca,'xtick',[])
    %                                 xlabel('ms')
    %                                 ylabel('Trials')
    %                                 % Save figure
    %                                 txt = strcat('SourceHeatmap_ROI',num2str(curR),'Chan',num2str(chans(ch)),'Curr',num2str(allCurrs(cu)),'Session',num2str(s)); 
    %                                 combineFileName = strcat(srcFld,txt,'.tif');
    %                                 saveas(curFig,combineFileName);





                                    % % plot individual trace data and subsequent averages
                                    % curFig = figure();
                                    % hold on
                                    % legSess = [];
                                    % for z = 1:size(midPostAve,1)
                                    %    pre = flip(preAve(z,:));
                                    %    midPost = midPostAve(z,:);
                                    %    dff = [pre,midPost];
                                    %    x = 33*((1:numel(dff))-numel(pre));
                                    %    plot(x,dff,'color',[0 0 0 0.2],'linewidth',2)
                                    % end
                                    % 
                                    % dff = [flip(preTrace),midPostTrace];
                                    % x = 33*((1:numel(dff))-numel(preTrace));
                                    % plot(x,dff,'linewidth',3,'color',[0 0 0])
                                    % xlim([-1000 4000])
                                    % xline(0,'--');
                                    % yline(curSTD)
                                    % xlabel('ms')
                                    % ylabel('\DeltaF/F')
                                    % % txt = strcat('Source_ROI',num2str(curR),'Chan',num2str(chans(ch)),'Curr',num2str(allCurrs(cu)),'Session',num2str(s)); 
                                    % % combineFileName = strcat(srcFld,txt,'.tif');
                                    % % saveas(curFig,combineFileName);
                                    % 
                                    % % yyaxis right
                                    % spikesAves = mean(midPostAveSpikes,1)>=0.5;
                                    % spikeInds = find(spikesAves);
                                    % for si = 1:numel(spikeInds)
                                    %     xline(spikeInds(si)*33,'r');
                                    % end
                                    % % h = bar((1:numel(spikesAves))*33, spikesAves, 0.3,'FaceColor',[1 0 0]);
                                    % % % h.EdgeColor = 'none';
                                    % % ylim([0 1])
                                    % % % ylabel('Spiking Probability')
                                    % % ax = gca;
                                    % % ax.YAxis(2).Color = 'r';
                                    % % Save figure
                                    % txt = strcat('SourceDecon_ROI',num2str(curR),'Chan',num2str(chans(ch)),'Curr',num2str(allCurrs(cu)),'Session',num2str(s)); 
                                    % combineFileName = strcat(srcFld,txt,'.tif');
                                    % saveas(curFig,combineFileName);
                                    % combineFileName = strcat(srcFld,txt,'.fig');
                                    % saveas(curFig,combineFileName);
                                    % close(curFig)
                                end
                            end
                        end
                    end
                end
            end
        end
        
        
        
        
        
        
        
        
        if(figGen==1)
            for ch = 1:numel(chans)
                for cu = 1:numel(allCurrs)
                    runningPre = zeros(1,1);
                    runningMid = zeros(1,1);
                    curV=0;

                    for s = 1:numel(files2Load)
                        if(currRois(s)~=0)
                            curData = AllData{s};
                            curAllTraces = curData.allTraces;
                            curAllTracesClean = curData.cleanTraces;
                            regions = curData.regions;
                            trialInfo = curData.TrialInfo;
                            spikes = curData.deconvolute;
                            curRegSTDs = regionSTDs{s};
                            curRegMeans = regionMeans{s};

                            m1 = trialInfo(:,1) == chans(ch);
                            m2 = trialInfo(:,3) == allCurrs(cu);
                            inds = find(m1 & m2);
                            if(numel(inds)>0)
                                preAve = NaN(numel(inds),10);
                                midPostAve = NaN(numel(inds),10);
                                midOnlyAve = NaN(numel(inds),10);
                                preAveSpikes = NaN(numel(inds),10);
                                midPostAveSpikes = NaN(numel(inds),10);
                                midOnlyAveSpikes = NaN(numel(inds),10);
                                for i = 1:numel(inds)
                                   pre = (curAllTraces{inds(i),1}-regionMeans{s})./regionMeans{s};
                                   midPost = ([curAllTraces{inds(i),2},curAllTraces{inds(i),3}]-regionMeans{s})./regionMeans{s};
                                   midOnly = (curAllTraces{inds(i),2}-regionMeans{s})./regionMeans{s};

                                   midOnlyAve(i,size(midOnly,2))=NaN;
                                   midOnlyAve(i,1:size(midOnly,2))=smooth(midOnly(currRois(s),:),5);
                                   preAve(i,size(pre,2))=NaN;
                                   preAve(i,1:size(pre,2))=smooth(flip(pre(currRois(s),:)),5);
                                   midPostAve(i,size(midPost,2))=NaN;
                                   midPostAve(i,1:size(midPost,2))=smooth(midPost(currRois(s),:),5);

                                end

                                % check if current trace has sufficently amplitude
                                % above threshold to be considered a spike.  This
                                % is to reduce the amount of high noise being
                                % processed.
                                curMidTrace = mean(midOnlyAve,1,'omitnan');
                                curSTD = curRegSTDs(currRois(s))./curRegMeans(currRois(s));
                                % activeTP = curMidTrace > (1*curSTD);
                                activeTP = curMidTrace > (2*curSTD);%-curRegMeans(currRois(s));
                                activeTP2 = curMidTrace > (curSTD);%-curRegMeans(currRois(s));

                                % find largest continous region in mid trace
                                d = [true, diff(activeTP) ~= 0];   % TRUE if values change
                                b = activeTP(d);                   % Elements without repetitions
                                n = diff(find([d, true]));  % Number of repetitions
                                n = n(b == 1);   % Care for runs of 1's only
                                m = max(n);
                                if(n > 5)
                                    sTime = find(activeTP2,1,'first')*33;
                                    if(sTime < 700) % start time must be before end of stimulation
                                        if(size(runningPre,2)>size(preAve,2))
                                            preAve(size(preAve,1),size(runningPre,2))=0;
                                        else
                                            runningPre(size(runningPre,1),size(preAve,2))=0;
                                        end
                                        if(size(runningMid,2)>size(midPostAve,2))
                                            midPostAve(size(midPostAve,1),size(runningMid,2))=0;
                                        else
                                            runningMid(size(runningMid,1),size(midPostAve,2))=0;
                                        end

                                        runningPre = [runningPre; preAve];
                                        runningMid = [runningMid; midPostAve];
                                        runningPre(size(runningPre,1)+1,:)=min([min(runningMid,[],'all'),-0.5]);
                                        runningMid(size(runningMid,1)+1,:)= min([min(runningMid,[],'all'),-0.5]);
                                        curV = curV+1;
                                    end
                                end
                            end
                        end
                    end
                    if(figGen==1)
                        if(curV > minSetCnt)
                            % plot individual trace data and subsequent averages
                            curFig = figure();
                            dff = [flip(runningPre(:,1:min([30,size(runningPre,1)])),2),runningMid(:,1:min([121,size(runningMid,1)]))];
                            h = heatmap(dff);
                            h.GridVisible = 'off';
                            colormap(parula)
                            xlabel('ms')
                            ylabel('Trials')
                            % Save figure
                            txt = strcat('SourceHeatmap-ROI',num2str(curR),'Chan',num2str(chans(ch)),'Curr',num2str(allCurrs(cu))); 
                            title(txt)
                            combineFileName = strcat(srcFld,txt,'.tif');
                            % saveas(curFig,combineFileName);
                            combineFileName = strcat(srcFld,txt,'.fig');
                            % saveas(curFig,combineFileName);
                            close(curFig)
                        end
                    end
                end
            end
        end
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %%
        for ch = 1:numel(chans)
            for cu = 1:numel(allCurrs)
                traceplotData = cell(numel(files2Load),1);
                if(figGen==1)
                    curFig = figure();
                end
                for s = 1:numel(files2Load)
                    subplot(1,numel(files2Load),s)
                    runningPre = zeros(1,1);
                    runningMid = zeros(1,1);
                    
                    for r = 1:numel(roiSets)
                        currRois = roiSets{r};
                        if(currRois(s)~=0)
                            curData = AllData{s};
                            curAllTraces = curData.allTraces;
                            curAllTracesClean = curData.cleanTraces;
                            regions = curData.regions;
                            trialInfo = curData.TrialInfo;
                            spikes = curData.deconvolute;
                            curRegSTDs = regionSTDs{s};
                            curRegMeans = regionMeans{s};

                            m1 = trialInfo(:,1) == chans(ch);
                            m2 = trialInfo(:,3) == allCurrs(cu);
                            inds = find(m1 & m2);
                            if(numel(inds)>0)
                                preAve = NaN(numel(inds),10);
                                midPostAve = NaN(numel(inds),10);
                                midOnlyAve = NaN(numel(inds),10);
                                preAveSpikes = NaN(numel(inds),10);
                                midPostAveSpikes = NaN(numel(inds),10);
                                midOnlyAveSpikes = NaN(numel(inds),10);
                                for i = 1:numel(inds)
                                   pre = (curAllTraces{inds(i),1}-regionMeans{s})./regionMeans{s};
                                   midPost = ([curAllTraces{inds(i),2},curAllTraces{inds(i),3}]-regionMeans{s})./regionMeans{s};
                                   midOnly = (curAllTraces{inds(i),2}-regionMeans{s})./regionMeans{s};

                                   midOnlyAve(i,size(midOnly,2))=NaN;
                                   midOnlyAve(i,1:size(midOnly,2))=smooth(midOnly(currRois(s),:),5);
                                   preAve(i,size(pre,2))=NaN;
                                   preAve(i,1:size(pre,2))=smooth(flip(pre(currRois(s),:)),5);
                                   midPostAve(i,size(midPost,2))=NaN;
                                   midPostAve(i,1:size(midPost,2))=smooth(midPost(currRois(s),:),5);

                                end

                                % check if current trace has sufficently amplitude
                                % above threshold to be considered a spike.  This
                                % is to reduce the amount of high noise being
                                % processed.
                                curMidTrace = mean(midOnlyAve,1,'omitnan');
                                curSTD = curRegSTDs(currRois(s))./curRegMeans(currRois(s));
                                % activeTP = curMidTrace > (1*curSTD);
                                activeTP = curMidTrace > (2*curSTD);
                                % activeTP2 = curMidTrace > (curSTD)-curRegMeans(currRois(s));

                                % find largest continous region in mid trace
                                d = [true, diff(activeTP) ~= 0];   % TRUE if values change
                                b = activeTP(d);                   % Elements without repetitions
                                n = diff(find([d, true]));  % Number of repetitions
                                n = n(b == 1);   % Care for runs of 1's only
                                m = max(n);

                                preAveM = mean(preAve,1);
                                midPostAveM = mean(midPostAve,1);
                                
                                if(size(runningPre,2)>size(preAveM,2))
                                    preAveM(size(preAveM,1),size(runningPre,2))=0;
                                else
                                    runningPre(size(runningPre,1),size(preAveM,2))=0;
                                end
                                if(size(runningMid,2)>size(midPostAveM,2))
                                    midPostAveM(size(midPostAveM,1),size(runningMid,2))=0;
                                else
                                    runningMid(size(runningMid,1),size(midPostAveM,2))=0;
                                end

                                runningPre = [runningPre; preAveM];
                                runningMid = [runningMid; midPostAveM];
                            end
                        else
                            preAveM=0;
                            midPostAveM=0;
                            if(size(runningPre,2)>size(preAveM,2))
                                preAveM(size(preAveM,1),size(runningPre,2))=0;
                            else
                                runningPre(size(runningPre,1),size(preAveM,2))=0;
                            end
                            if(size(runningMid,2)>size(midPostAveM,2))
                                midPostAveM(size(midPostAveM,1),size(runningMid,2))=0;
                            else
                                runningMid(size(runningMid,1),size(midPostAveM,2))=0;
                            end

                            runningPre = [runningPre; preAveM];
                            runningMid = [runningMid; midPostAveM];
                        end
                    end
                    dff = [flip(runningPre(:,1:min([30,size(runningPre,2)])),2),runningMid(:,1:min([121,size(runningMid,2)]))];
                    if(figGen==1)
    
                        h = heatmap(dff,'ColorLimits',[0 3]);
                        h.GridVisible = 'off';
                        colormap(parula)
                        h.XDisplayLabels = nan(size(h.XDisplayData));
                        h.YDisplayLabels = nan(size(h.YDisplayData));
                        if(s ~= numel(files2Load))
                            h.ColorbarVisible = 'off';
                        end
                    end
                    traceplotData{s} = mean(dff,1);
                    
                end

                if(figGen==1)
                    % Save figure
                    txt = strcat('LongitudinalHeatmaps-Chan',num2str(chans(ch)),'Curr',num2str(allCurrs(cu))); 
                    combineFileName = strcat(srcFld,txt,'.tif');
                    % saveas(curFig,combineFileName);
                    combineFileName = strcat(srcFld,txt,'.fig');
                    % saveas(curFig,combineFileName);
                    close(curFig)
                    
                    
                    
                    
                    
                    % Save figure
                    curFig = figure();
                    for s = 1:numel(files2Load)
                        subplot(1,numel(files2Load),s)
                        if(~isempty(traceplotData{s}))
                            plot(traceplotData{s})
                        end
                    end
                    
                    txt = strcat('LongitudinalHeatmapsPlOT-Chan',num2str(chans(ch)),'Curr',num2str(allCurrs(cu))); 
                    combineFileName = strcat(srcFld,txt,'.tif');
                    % saveas(curFig,combineFileName);
                    combineFileName = strcat(srcFld,txt,'.fig');
                    % saveas(curFig,combineFileName);
                    close(curFig)
                end
            end
        end
        
        %%
        

        if(figGen==1)
        for ch = 1:numel(chans)
            for cu = 1:numel(allCurrs)
                curFig = figure();
                traceplotData = cell(numel(files2Load),1);
                for s = 1:numel(files2Load)
                    subplot(1,numel(files2Load),s)
                    runningPre = zeros(1,1);
                    runningMid = zeros(1,1);
                    
                    for r = 1:numel(roiEndSets)
                        currRois = roiEndSets{r};%  construct path
                        currRoiPxs = roiEndSetsPx{r};
                        if(currRois(s)~=0)
                            curData = AllData{s};
                            curAllTraces = curData.allTraces;
                            curAllTracesClean = curData.cleanTraces;
                            regions = curData.regions;
                            trialInfo = curData.TrialInfo;
                            spikes = curData.deconvolute;
                            curRegSTDs = regionSTDs{s};
                            curRegMeans = regionMeans{s};

                            m1 = trialInfo(:,1) == chans(ch);
                            m2 = trialInfo(:,3) == allCurrs(cu);
                            inds = find(m1 & m2);
                            if(numel(inds)>0)
                                preAve = NaN(numel(inds),10);
                                midPostAve = NaN(numel(inds),10);
                                midOnlyAve = NaN(numel(inds),10);
                                preAveSpikes = NaN(numel(inds),10);
                                midPostAveSpikes = NaN(numel(inds),10);
                                midOnlyAveSpikes = NaN(numel(inds),10);
                                for i = 1:numel(inds)
                                   pre = (curAllTraces{inds(i),1}-regionMeans{s})./regionMeans{s};
                                   midPost = ([curAllTraces{inds(i),2},curAllTraces{inds(i),3}]-regionMeans{s})./regionMeans{s};
                                   midOnly = (curAllTraces{inds(i),2}-regionMeans{s})./regionMeans{s};

                                   midOnlyAve(i,size(midOnly,2))=NaN;
                                   midOnlyAve(i,1:size(midOnly,2))=smooth(midOnly(currRois(s),:),5);
                                   preAve(i,size(pre,2))=NaN;
                                   preAve(i,1:size(pre,2))=smooth(flip(pre(currRois(s),:)),5);
                                   midPostAve(i,size(midPost,2))=NaN;
                                   midPostAve(i,1:size(midPost,2))=smooth(midPost(currRois(s),:),5);

                                end

                                % % check if current trace has sufficently amplitude
                                % % above threshold to be considered a spike.  This
                                % % is to reduce the amount of high noise being
                                % % processed.
                                % curMidTrace = mean(midOnlyAve,1,'omitnan');
                                % curSTD = curRegSTDs(currRois(s))./curRegMeans(currRois(s));
                                % % activeTP = curMidTrace > (1*curSTD);
                                % activeTP = curMidTrace > (2*curSTD);%-curRegMeans(currRois(s));
                                % % activeTP2 = curMidTrace > (curSTD)-curRegMeans(currRois(s));
                                % 
                                % % find largest continous region in mid trace
                                % d = [true, diff(activeTP) ~= 0];   % TRUE if values change
                                % b = activeTP(d);                   % Elements without repetitions
                                % n = diff(find([d, true]));  % Number of repetitions
                                % n = n(b == 1);   % Care for runs of 1's only
                                % m = max(n);

                                preAveM = mean(preAve,1);
                                midPostAveM = mean(midPostAve,1);
                                
                                if(size(runningPre,2)>size(preAveM,2))
                                    preAveM(size(preAveM,1),size(runningPre,2))=0;
                                else
                                    runningPre(size(runningPre,1),size(preAveM,2))=0;
                                end
                                if(size(runningMid,2)>size(midPostAveM,2))
                                    midPostAveM(size(midPostAveM,1),size(runningMid,2))=0;
                                else
                                    runningMid(size(runningMid,1),size(midPostAveM,2))=0;
                                end

                                runningPre = [runningPre; preAveM];
                                runningMid = [runningMid; midPostAveM];
                            end
                        else
                            
                            
                            
%                             
%                             % Load timing key data
%                             selpath = files2Load{s};
%                             idcs = strfind(selpath,'\');
%                             pullPath = selpath(1:idcs(numel(idcs)-1)-1);
%                             keypathO = strcat(pullPath,'\KEY2D');
%                             keyListing = dir(keypathO);
%                             keyFile = keyListing(3).name; % There should only be one key per session, grab it
%                             dataCellArr = load(strcat(keypathO,'\',keyFile));
%                             dataCellArr = dataCellArr.dataCellArr;
%                             eCell = zeros(numel(dataCellArr),1);
%                             for k = 1:numel(dataCellArr)
%                                 if(isempty(dataCellArr{k}))
%                                     eCell(k)=1;
%                                 end
%                             end
%                             dataCellArr(eCell==1)=[];
%                             
%                             trialMask = false(512,512);                            
%                             for m = 1:size(currRoiPxs,1)
%                                 curPixels = ind2sub([512, 512],currRoiPxs(m));
%                                 trialMask(curPixels(2),curPixels(1))=1;
%                             end
%                             
%                             for trialNum = 1:numel(dataCellArr) % process each relevent trial
% 
%                                 % Valid trial, process data. Trials with error flags will be ignored/discarded
%                                 trial_i = dataCellArr{trialNum};
%                                 validTrial = trial_i{7};
%                                 
%                                 if(validTrial)
%                                     trialCount = trialCount+1; % increment trial count
% 
%                                     % load imaged volume set for trial
%                                     destFile = strcat(pullPath,'PROCCESSED-PLANAR\Data_PlanarProcessed',num2str(trialNum),'.mat');
%                                     fileData = load(destFile,'imgBlock');
%                                     curSeq = fileData.imgBlock;
% 
%                                     % Quantify flouresence traces for given trial
%                                     roiTracesAll = zeros(1,size(curSeq,3));
% 
%                                     % process all frame data
%                                     for z = 1:size(curSeq,3)
%                                         cS = curSeq(:,:,z);
%                                         vals = cS(trialMask);
%                                         roiTracesAll(z) = mean(vals);
%                                     end
%                                 end
%                             end  
%                             
%                             
%                             
%                             trialCount = 0;
%                             for trialNum = 1:numel(dataCellArr) % process each relevent trial
% 
%                                 % Valid trial, process data. Trials with error flags will be ignored/discarded
%                                 trial_i = dataCellArr{trialNum};
%                                 validTrial = trial_i{7};
%                                 
%                                 if(validTrial)
%                                     if(trial_i{2} == chans(ch))
%                                         if(trial_i{1} == allCurrs(cu))
%                                             trialCount = trialCount+1; % increment trial count
% 
%                                             % load imaged volume set for trial
%                                             destFile = strcat(pullPath,'PROCCESSED-PLANAR\Data_PlanarProcessed',num2str(trialNum),'.mat');
%                                             fileData = load(destFile,'imgBlock');
%                                             curSeq = fileData.imgBlock;
% 
%                                             % Collect trial specific timing data
%                                             if(numel(trial_i{10})<2)
%                                                 stim_ts = [seconds(trial_i{10}), seconds(trial_i{10}+0.6)];
%                                             else
%                                                 stim_ts = seconds(trial_i{10});
%                                             end
%                                             imgStarts = seconds(trial_i{9});
%                                             imgStarts = imgStarts(1:size(curSeq,3)); % Trim Slices to limit of collected sequence
%                                             stimStart = stim_ts(1);
%                                             stimEnd = stim_ts(2) + seconds(1); % Add an additional 1000 ms to end of stim period to collect full flourescence activity
%                                             preSlices = imgStarts<stimStart;
%                                             postSlices = imgStarts>stimEnd;
%                                             midSlices = imgStarts>=stimStart;
%                                             b = imgStarts<=stimEnd;
%                                             midSlices(b==0)=0;
% 
%                                             % Assign frames to pre, mid, and post stim regions
%                                             preFrames = curSeq(:,:,preSlices);
%                                             midFrames = curSeq(:,:,midSlices);
%                                             postFrames = curSeq(:,:,postSlices);
% 
%                                             % Quantify flouresence traces for given trial
%                                             roiTracesPre = zeros(1,size(preSlices,3));
%                                             roiTracesMid = zeros(1,size(midSlices,3));
%                                             roiTracesPost = zeros(1,size(postSlices,3));
%                                             
% 
%                                             % process pre-stim data
%                                             for z = 1:size(preFrames,3)
%                                                 cS = preFrames(:,:,z);
%                                                 vals = cS(trialMask);
%                                                 roiTracesPre(z) = mean(vals);
%                                             end
% 
%                                             % process mid-stim data
%                                             for z = 1:size(midFrames,3)
%                                                 cS = midFrames(:,:,z);
%                                                 vals = cS(trialMask);
%                                                 roiTracesMid(z) = mean(vals);
%                                             end
% 
%                                             % process post-stim data
%                                             for z = 1:size(postFrames,3)
%                                                 cS = postFrames(:,:,z);
%                                                 vals = cS(trialMask);
%                                                 roiTracesPost(z) = mean(vals);
%                                             end
%                                             
%                                             allTraces{trialNum,1}=roiTracesPre;
%                                             allTraces{trialNum,2}=roiTracesMid;
%                                             allTraces{trialNum,3}=roiTracesPost;
%                                         end
%                                     end
%                                 end
%                             end           
                            

                            preAveM=0;
                            midPostAveM=0;
                            if(size(runningPre,2)>size(preAveM,2))
                                preAveM(size(preAveM,1),size(runningPre,2))=0;
                            else
                                runningPre(size(runningPre,1),size(preAveM,2))=0;
                            end
                            if(size(runningMid,2)>size(midPostAveM,2))
                                midPostAveM(size(midPostAveM,1),size(runningMid,2))=0;
                            else
                                runningMid(size(runningMid,1),size(midPostAveM,2))=0;
                            end

                            runningPre = [runningPre; preAveM];
                            runningMid = [runningMid; midPostAveM];
                        end
                    end
                    dff = [flip(runningPre(:,1:min([30,size(runningPre,2)])),2),runningMid(:,1:min([121,size(runningMid,2)]))];
                    h = heatmap(dff,'ColorLimits',[0 3]);
                    h.GridVisible = 'off';
                    colormap(parula)
                    h.XDisplayLabels = nan(size(h.XDisplayData));
                    h.YDisplayLabels = nan(size(h.YDisplayData));
                    if(s ~= numel(files2Load))
                        h.ColorbarVisible = 'off';
                    end
                    traceplotData{s} = mean(dff,1);
                end

                % Save figure
                txt = strcat('BaselineLongitudinalHeatmaps-Chan',num2str(chans(ch)),'Curr',num2str(allCurrs(cu))); 
                combineFileName = strcat(srcFld,txt,'.tif');
                % saveas(curFig,combineFileName);
                combineFileName = strcat(srcFld,txt,'.fig');
                % saveas(curFig,combineFileName);
                close(curFig)
                
                
                
                % Save figure
                curFig = figure();
                for s = 1:numel(files2Load)
                    subplot(1,numel(files2Load),s)
                    if(~isempty(traceplotData{s}))
                        plot(traceplotData{s})
                    end
                end
                txt = strcat('BaselineLongitudinalHeatmapPlot-Chan',num2str(chans(ch)),'Curr',num2str(allCurrs(cu))); 
                combineFileName = strcat(srcFld,txt,'.tif');
                % saveas(curFig,combineFileName);
                combineFileName = strcat(srcFld,txt,'.fig');
                % saveas(curFig,combineFileName);
                close(curFig)
            end
        end
        end

    %% Calculate longitudinal trends for all metrics measured
    allTrendsOnset = NaN(1,numel(files2Load));
    allTrendsDur = NaN(1,numel(files2Load));
    allTrendsVol = NaN(1,numel(files2Load));
    allTrendsMax = NaN(1,numel(files2Load));
    allTrendsSpikes = NaN(1,numel(files2Load));
    allTrendsSpikesItime = NaN(1,numel(files2Load));
    allTrendsSpikesDuration = NaN(1,numel(files2Load));
    allTrendsSpikesPeaks = NaN(1,numel(files2Load));
    allTrendsSpikesMtime = NaN(1,numel(files2Load));
    allTrendsN20 = NaN(1,numel(files2Load));
    allTrendsRoiCh = NaN(1,3);
    tN = 0;
    for curR = 1:size(ROIs2Compare,1)
        for ch = 1:numel(chans)
            for cu = 1:numel(allCurrs)
               curOnset = NaN(1,numel(files2Load));
               curDur = NaN(1,numel(files2Load));
               curVol = NaN(1,numel(files2Load));
               curMax = NaN(1,numel(files2Load));
               nT20 = NaN(1,numel(files2Load));
               curSpikes =  NaN(1,numel(files2Load));
               curSpikesItime =  NaN(1,numel(files2Load));
               curSpikesMtime =  NaN(1,numel(files2Load));
               curSpikesDuration =  NaN(1,numel(files2Load));
               curSpikesPeaks =  NaN(1,numel(files2Load));
               cnt=0;
               for s = 1:numel(files2Load)
                   if(~isempty(midTraces{curR,s,ch,cu}))
                       curOnset(s) = spikeOnset(curR,s,ch,cu);
                       curDur(s) = spikeDur(curR,s,ch,cu);
                       curVol(s) = spikeVol(curR,s,ch,cu);
                       curMax(s) = spikeMax(curR,s,ch,cu);
                       curSpikes(s) = spikeCount(curR,s,ch,cu);
                       curSpikesItime(s) = spikeInitialTime(curR,s,ch,cu);
                       curSpikesMtime(s) = spikeMaxTime(curR,s,ch,cu);
                       curSpikesDuration(s) = spikeDeconDur(curR,s,ch,cu);
                       curSpikesPeaks(s) = spikePeaks(curR,s,ch,cu);
                       nT20(s) = spikeN20(curR,s,ch,cu);
                       cnt = cnt+1;
                   end
               end
               if(cnt>minSetCnt)
                   tN = tN+1;
                   allTrendsOnset(tN,numel(files2Load))=NaN; % expand array
                   allTrendsDur(tN,numel(files2Load))=NaN;
                   allTrendsVol(tN,numel(files2Load))=NaN;
                   allTrendsMax(tN,numel(files2Load))=NaN;
                   allTrendsN20(tN,numel(files2Load))=NaN;
                   allTrendsOnset(tN,:)=curOnset; % save valid trendline
                   allTrendsDur(tN,:)=curDur;
                   allTrendsVol(tN,:)=curVol;
                   allTrendsMax(tN,:)=curMax;
                   allTrendsSpikes(tN,:)=curSpikes;
                   allTrendsSpikesItime(tN,:)=curSpikesItime;
                   allTrendsSpikesMtime(tN,:)=curSpikesMtime;
                   allTrendsSpikesDuration(tN,:)=curSpikesDuration;
                   allTrendsSpikesPeaks(tN,:)=curSpikesPeaks;
                   allTrendsN20(tN,:) = nT20;
                   allTrendsRoiCh(tN,1) = curR;
                   allTrendsRoiCh(tN,2) = chans(ch);
                   allTrendsRoiCh(tN,3) = allCurrs(cu);
                   
               end
            end
        end
    end

    weeksTrained = ceil(daysTrained./7);





%     %% Plot active region means and STDs for all ROIS
%     for r = 1:numel(roiSets)
%         curFig = figure();
%         hold on
%         curROIInds = roiSets{r};
%         targMeans = NaN(numel(files2Load),1);
%         targSTDs = NaN(numel(files2Load),1);
%         for f = 1:numel(files2Load) 
%             curRegSTDs = regionSTDs{f};
%             curRegMeans = regionMeans{f};
% 
%             curSessROIind = curROIInds(f);
%             if(curSessROIind>0)
%                 targMeans(f) = curRegMeans(curSessROIind);
%                 targSTDs(f) = curRegSTDs(curSessROIind);
%             end
%         end
%         plot(targMeans,'-o')
%         yyaxis right
%         plot(targSTDs,'--o')
%         legend({'Mean';'STD'})
% 
%         title(strcat('Baseline Trends ROI',num2str(r)));
% %         combineFileName = strcat(curFld,'BaselineROIs-',num2str(r),'.fig');
% %         saveas(curFig,combineFileName);  
%         combineFileName = strcat(curFld,'BaselineROIs-',num2str(r),'.tif');
%         saveas(curFig,combineFileName);
%         close(curFig)
%     end
%     
    
%% Genrate traces for ROIs non-consistently active
% 
%     nonConROIs2Compare = zeros(numel(roiEndSets),numel(files2Load));
%     for r = 1:numel(roiEndSets) % convert keys to single matrix
%         nonConROIs2Compare(r,:) = roiEndSets{r};
%     end
% 
%     % Calculate longitudinal trends for all metrics measured
%     nonConTrendsOnset = NaN(1,numel(files2Load));
%     nonConTrendsDur = NaN(1,numel(files2Load));
%     nonConTrendsVol = NaN(1,numel(files2Load));
%     nonConTrendsMax = NaN(1,numel(files2Load));
%     nonConTrendsSpikes = NaN(1,numel(files2Load));
%     nonConTrendsSpikesItime = NaN(1,numel(files2Load));
%     nonConTrendsSpikesDuration = NaN(1,numel(files2Load));
%     nonConTrendsSpikesPeaks = NaN(1,numel(files2Load));
%     nonConTrendsSpikesMtime = NaN(1,numel(files2Load));
%     nonConTrendsN20 = NaN(1,numel(files2Load));
%     nonConTrendsRoiCh = NaN(1,3);
%     
%     tN = 0;
%     for curR = 1:size(nonConROIs2Compare,1)
%         for ch = 1:numel(chans)
%             for cu = 1:numel(allCurrs)
%                curOnset = NaN(1,numel(files2Load));
%                curDur = NaN(1,numel(files2Load));
%                curVol = NaN(1,numel(files2Load));
%                curMax = NaN(1,numel(files2Load));
%                nT20 = NaN(1,numel(files2Load));
%                curSpikes =  NaN(1,numel(files2Load));
%                curSpikesItime =  NaN(1,numel(files2Load));
%                curSpikesMtime =  NaN(1,numel(files2Load));
%                curSpikesDuration =  NaN(1,numel(files2Load));
%                curSpikesPeaks =  NaN(1,numel(files2Load));
%                cnt=0;
%                for s = 1:numel(files2Load)
%                    if(~isempty(endMidTraces{curR,s,ch,cu}))
%                        curOnset(s) = spikeOnset(curR,s,ch,cu);
%                        curDur(s) = spikeDur(curR,s,ch,cu);
%                        curVol(s) = spikeVol(curR,s,ch,cu);
%                        curMax(s) = spikeMax(curR,s,ch,cu);
%                        curSpikes(s) = spikeCount(curR,s,ch,cu);
%                        curSpikesItime(s) = spikeInitialTime(curR,s,ch,cu);
%                        curSpikesMtime(s) = spikeMaxTime(curR,s,ch,cu);
%                        curSpikesDuration(s) = spikeDeconDur(curR,s,ch,cu);
%                        curSpikesPeaks(s) = spikePeaks(curR,s,ch,cu);
%                        nT20(s) = spikeN20(curR,s,ch,cu);
%                        cnt = cnt+1;
%                    end
%                end
%                if(cnt<minSetCnt)
%                    tN = tN+1;
%                    nonConTrendsOnset(tN,numel(files2Load))=NaN; % expand array
%                    nonConTrendsDur(tN,numel(files2Load))=NaN;
%                    nonConTrendsVol(tN,numel(files2Load))=NaN;
%                    nonConTrendsMax(tN,numel(files2Load))=NaN;
%                    nonConTrendsOnset(tN,:)=curOnset; % save valid trendline
%                    nonConTrendsRoiCh(tN,1) = curR;
%                    nonConTrendsRoiCh(tN,2) = chans(ch);
%                    nonConTrendsRoiCh(tN,3) = allCurrs(cu);
%                    nonConTrendsN20(tN,numel(files2Load))=NaN;
%                    nonConTrendsDur(tN,:)=curDur;
%                    nonConTrendsVol(tN,:)=curVol;
%                    nonConTrendsMax(tN,:)=curMax;
%                    nonConTrendsSpikes(tN,:)=curSpikes;
%                    nonConTrendsSpikesItime(tN,:)=curSpikesItime;
%                    nonConTrendsSpikesMtime(tN,:)=curSpikesMtime;
%                    nonConTrendsSpikesDuration(tN,:)=curSpikesDuration;
%                    nonConTrendsSpikesPeaks(tN,:)=curSpikesPeaks;
%                    nonConTrendsN20(tN,:) = nT20;
%                end
%             end
%         end
%     end
%     


%% Genrate traces for ROIs active only in later sessions

    endROIs2Compare = zeros(numel(roiEndSets),numel(files2Load));
    for r = 1:numel(roiEndSets) % convert keys to single matrix
        endROIs2Compare(r,:) = roiEndSets{r};
    end

    % calculate trace means and stds
    % TODO - report baseline trace means over each session, i.e. see if there is a drifting baseline in each session or for each roi
    allCurrs = [];
    endRegionMeans = cell(numel(files2Load),1);
    endRegionSTDs = cell(numel(files2Load),1);
    endRegionPreStimMeans = cell(numel(files2Load),1);
    endRegionPreStimSTDs = cell(numel(files2Load),1);
    for s = 1:numel(files2Load)
        curData = AllData{s};
        trialInfo = curData.TrialInfo;
        currs = unique(trialInfo(:,3));
        regions = curData.regions;
        allCurrs = [allCurrs; currs];
        cN = 0;
        cPN = 0;
        
        curAllTraces = curData.allTraces;
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
            endRegionMeans{s} = mean(completeTrace,2);
            endRegionSTDs{s} = std(completeTrace,[],2);
            
            endRegionPreStimMeans{s} = mean(completePreTrace,2);
            endRegionPreStimSTDs{s} = std(completePreTrace,[],2);
        end
    end
    allCurrs = unique(allCurrs);
    allCurrs(allCurrs==0)=[];


    endMidTraces = cell(size(endROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    endPreTraces = cell(size(endROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    endSpikeTraces = cell(size(endROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    endPreSpikeTraces = cell(size(endROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    midSTD = cell(size(endROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeOnset = zeros(size(endROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeVol = zeros(size(endROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeDur = zeros(size(endROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeMax = zeros(size(endROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeN20 = zeros(size(endROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeCount = zeros(size(endROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeMaxTime = zeros(size(endROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeInitialTime = zeros(size(endROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeDeconDur = zeros(size(endROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikePeaks = zeros(size(endROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    for curR = 1:size(endROIs2Compare,1)
        currRois = endROIs2Compare(curR,:);
        for s = 1:numel(files2Load)
            if(currRois(s)~=0)
                curData = AllData{s};
                curAllTraces = curData.allTraces;
                curAllTracesClean = curData.cleanTraces;
                regions = curData.regions;
                trialInfo = curData.TrialInfo;
                spikes = curData.deconvolute;
                curRegSTDs = regionSTDs{s};
                curRegMeans = regionMeans{s};

                for ch = 1:numel(chans)
                    for cu = 1:numel(allCurrs)
                        m1 = trialInfo(:,1) == chans(ch);
                        m2 = trialInfo(:,3) == allCurrs(cu);
                        inds = find(m1 & m2);
                        if(numel(inds)>0)
                            preAve = NaN(numel(inds),10);
                            midPostAve = NaN(numel(inds),10);
                            midOnlyAve = NaN(numel(inds),10);
                            preAveSpikes = NaN(numel(inds),10);
                            midPostAveSpikes = NaN(numel(inds),10);
                            midOnlyAveSpikes = NaN(numel(inds),10);
                            for i = 1:numel(inds)
                               pre = (curAllTraces{inds(i),1}-regionMeans{s})./regionMeans{s};
                               midPost = ([curAllTraces{inds(i),2},curAllTraces{inds(i),3}]-regionMeans{s})./regionMeans{s};

                               preAve(i,size(pre,2))=NaN;
                               preAve(i,1:size(pre,2))=smooth(flip(pre(currRois(s),:)),5);
                               midPostAve(i,size(midPost,2))=NaN;
                               midPostAve(i,1:size(midPost,2))=smooth(midPost(currRois(s),:),5);
                               midOnly = (curAllTraces{inds(i),2}-regionMeans{s})./regionMeans{s};
                               midOnlyAve(i,size(midOnly,2))=NaN;
                               midOnlyAve(i,1:size(midOnly,2))=smooth(midOnly(currRois(s),:),5);
                               
                               
                               preSpikes = spikes{inds(i),1};
                               midOnlySpikes = spikes{inds(i),2};
                               midPostSpikes = [spikes{inds(i),2},spikes{inds(i),3}];
                               
                               preAveSpikes(i,size(preSpikes,2))=NaN;
                               preAveSpikes(i,1:size(preSpikes,2))=flip(preSpikes(currRois(s),:));
                               midPostAveSpikes(i,size(midPostSpikes,2))=NaN;
                               midPostAveSpikes(i,1:size(midPostSpikes,2))=midPostSpikes(currRois(s),:);
                               midOnlyAveSpikes(i,size(midOnlySpikes,2))=NaN;
                               midOnlyAveSpikes(i,1:size(midOnlySpikes,2))=midOnlySpikes(currRois(s),:);
                            end

                            % check if current trace has sufficently amplitude
                            % above threshold to be considered a spike.  This
                            % is to reduce the amount of high noise being
                            % processed.
                            curMidTrace = mean(midOnlyAve,1,'omitnan');
                            curSTD = curRegSTDs(currRois(s))./curRegMeans(currRois(s));
                            % activeTP = curMidTrace > (1*curSTD);
                            activeTP = curMidTrace > (2*curSTD);
                            activeTP2 = curMidTrace > (curSTD);

                            % find largest continous region in mid trace
                            d = [true, diff(activeTP) ~= 0];   % TRUE if values change
                            b = activeTP(d);                   % Elements without repetitions
                            n = diff(find([d, true]));  % Number of repetitions
                            n = n(b == 1);   % Care for runs of 1's only
                            m = max(n);
                            if(n > minSetCnt)
                                sTime = find(activeTP2,1,'first')*33;
                                if(sTime < 700) % start time must be before end of stimulation

                                    midPostTrace = mean(midPostAve,1,'omitnan');
                                    preTrace = mean(preAve,1,'omitnan');

                                    endPreTraces{curR,s,ch,cu} = flip(preTrace);
                                    endMidTraces{curR,s,ch,cu} = midPostTrace; % mid and post trace
                                    endSpikeTraces{curR,s,ch,cu} = mean(midPostAveSpikes,1,'omitnan');

                                    % find initial spike time from Deconvolved data
                                    iSpikes = zeros(size(midPostAveSpikes,1),1);
                                    eSpikes = zeros(size(midPostAveSpikes,1),1);
                                    for g = 1:size(midPostAveSpikes,1)
                                        if(isempty(find(midPostAveSpikes(g,:),1)))
                                            iSpikes(g) = NaN;
                                            eSpikes(g) = NaN;
                                        else
                                            iSpikes(g) = find(midPostAveSpikes(g,:),1);
                                            eSpikes(g) = find(midPostAveSpikes(g,:),1,'last');
                                        end
                                    end
                                    spikeInitialTime(curR,s,ch,cu) = mean(iSpikes,'omitnan')*33;


                                    spikeDeconDur(curR,s,ch,cu) =  mean(eSpikes - iSpikes,'omitnan')*33;
                                    temp = smooth( mean(midPostAveSpikes,1,'omitnan'),5);
                                    spikePeaks(curR,s,ch,cu) = numel(findpeaks(temp,'MinPeakProminence',0.05));

                                    endMidSpikeTraces = mean(midOnlyAveSpikes,1,'omitnan');
                                    endPreSpikeTraces{curR,s,ch,cu} = flip(mean(preAveSpikes,1,'omitnan'));
                                    % Determine onset and end times for current
                                    % trace

                                    % save metrics of trace for analysis
                                    spikeOnset(curR,s,ch,cu) = find(activeTP2,1,'first')*33; 
                                    spikeDur(curR,s,ch,cu) = (find(activeTP2,1,'last') - find(activeTP2,1,'first'))*33;
                                    spikeVol(curR,s,ch,cu) = sum(midPostTrace(find(activeTP2,1,'first') : find(activeTP2,1,'last')));
                                    spikeMax(curR,s,ch,cu) = max(midPostTrace);
                                    spikeCount(curR,s,ch,cu) = sum(endMidSpikeTraces); % number of spikes during stimulation period


    %                                 % find max spike time
    %                                 [~,mInd] = max(midPostTrace);
    %                                 spikeMaxTime(curR,s,ch,cu) = mInd*33;
    %                                 midPostTraceN = midPostTrace./max(midPostTrace);
    %                                 iData = interp1(1:numel(midPostTraceN),midPostTraceN,1:0.1:numel(midPostTraceN)/2);
    %                                 iX = 1:0.1:numel(midPostTraceN)/2;
    %                                 iDataR = round(iData*10)./10;
    %                                 tInd = find(iDataR==0.2,1);
    %                                 spikeN20(curR,s,ch,cu) = iX(tInd)*33;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
        
        

    
    
    
    
    endMidTraces = cell(size(endROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    endPreTraces = cell(size(endROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeOnset = zeros(size(endROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeVol = zeros(size(endROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeDur = zeros(size(endROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    spikeMax = zeros(size(endROIs2Compare,1),numel(files2Load),numel(chans),numel(allCurrs));
    for curR = 1:size(endROIs2Compare,1)
        currRois = endROIs2Compare(curR,:);
        for s = 1:numel(files2Load)
            if(currRois(s)~=0)
                curData = AllData{s};
                curAllTraces = curData.allTraces;
                regions = curData.regions;
                trialInfo = curData.TrialInfo;
                curRegSTDs = endRegionSTDs{s};
                curRegMeans = endRegionMeans{s};

                for ch = 1:numel(chans)
                    for cu = 1:numel(allCurrs)
                        m1 = trialInfo(:,1) == chans(ch);
                        m2 = trialInfo(:,3) == allCurrs(cu);
                        inds = find(m1 & m2);
                        if(numel(inds)>0)
                            preAve = NaN(numel(inds),10);
                            midPostAve = NaN(numel(inds),10);
                            midOnlyAve = NaN(numel(inds),10);
                            for i = 1:numel(inds)
                               pre = (curAllTraces{inds(i),1}-endRegionMeans{s})./endRegionMeans{s};
                               midPost = ([curAllTraces{inds(i),2},curAllTraces{inds(i),3}]-endRegionMeans{s})./endRegionMeans{s};

                               preAve(i,size(pre,2))=NaN;
                               preAve(i,1:size(pre,2))=smooth(flip(pre(currRois(s),:)),5);
                               midPostAve(i,size(midPost,2))=NaN;
                               midPostAve(i,1:size(midPost,2))=smooth(midPost(currRois(s),:),5);
                               midOnly = (curAllTraces{inds(i),2}-endRegionMeans{s})./endRegionMeans{s};
                               midOnlyAve(i,size(midOnly,2))=NaN;
                               midOnlyAve(i,1:size(midOnly,2))=smooth(midOnly(currRois(s),:),5);
                            end

                            % check if current trace has sufficently amplitude
                            % above threshold to be considered a spike.  This
                            % is to reduce the amount of high noise being
                            % processed.
                            curMidTrace = mean(midOnlyAve,1,'omitnan');
                            curSTD = curRegSTDs(currRois(s))./curRegMeans(currRois(s));
                            % activeTP = curMidTrace > (1*curSTD);
                            activeTP = curMidTrace > (2*curSTD);%-curRegMeans(currRois(s));
                            activeTP2 = curMidTrace > (curSTD);%-curRegMeans(currRois(s));

                            % find largest continous region in mid trace
                            d = [true, diff(activeTP) ~= 0];   % TRUE if values change
                            b = activeTP(d);                   % Elements without repetitions
                            n = diff(find([d, true]));  % Number of repetitions
                            n = n(b == 1);   % Care for runs of 1's only
                            m = max(n);
                            if(n > minSetCnt)
                                sTime = find(activeTP2,1,'first')*33;
                                if(sTime < 700) % start time must be before end of stimulation
                                    endMidPostTrace = mean(midPostAve,1,'omitnan');
                                    endPreTraces{curR,s,ch,cu} = flip(mean(preAve,1,'omitnan'));
                                    endMidTraces{curR,s,ch,cu} = endMidPostTrace; % mid and post trace

                                    % Determine onset and end times for current
                                    % trace
                                    % activeTP2 = endMidPostTrace > curSTD;

                                    % save metrics of trace for analysis
                                    spikeOnset(curR,s,ch,cu) = find(activeTP2,1,'first')*33; 
                                    spikeDur(curR,s,ch,cu) = (find(activeTP2,1,'last') - find(activeTP2,1,'first'))*33;
                                    spikeVol(curR,s,ch,cu) = sum(endMidPostTrace(find(activeTP2,1,'first') : find(activeTP,1,'last')));
                                    spikeMax(curR,s,ch,cu) = max(endMidPostTrace);
                                    midPostTraceN = endMidPostTrace./max(endMidPostTrace);
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    % Calculate longitudinal trends for all metrics measured
    endTrendsOnset = NaN(1,numel(files2Load));
    endTrendsDur = NaN(1,numel(files2Load));
    endTrendsVol = NaN(1,numel(files2Load));
    endTrendsMax = NaN(1,numel(files2Load));
    endTrendsSpikes = NaN(1,numel(files2Load));
    endTrendsSpikesItime = NaN(1,numel(files2Load));
    endTrendsSpikesDuration = NaN(1,numel(files2Load));
    endTrendsSpikesPeaks = NaN(1,numel(files2Load));
    endTrendsSpikesMtime = NaN(1,numel(files2Load));
%     endTrendsN20 = NaN(1,numel(files2Load));
    endTrendsRoiCh = NaN(1,3);
    
    tN = 0;
    for curR = 1:size(endROIs2Compare,1)
        for ch = 1:numel(chans)
            for cu = 1:numel(allCurrs)
               curOnset = NaN(1,numel(files2Load));
               curDur = NaN(1,numel(files2Load));
               curVol = NaN(1,numel(files2Load));
               curMax = NaN(1,numel(files2Load));
               nT20 = NaN(1,numel(files2Load));
               curSpikes =  NaN(1,numel(files2Load));
               curSpikesItime =  NaN(1,numel(files2Load));
               curSpikesMtime =  NaN(1,numel(files2Load));
               curSpikesDuration =  NaN(1,numel(files2Load));
               curSpikesPeaks =  NaN(1,numel(files2Load));
               cnt=0;
               for s = 1:numel(files2Load)
                   if(~isempty(endMidTraces{curR,s,ch,cu}))
                       curOnset(s) = spikeOnset(curR,s,ch,cu);
                       curDur(s) = spikeDur(curR,s,ch,cu);
                       curVol(s) = spikeVol(curR,s,ch,cu);
                       curMax(s) = spikeMax(curR,s,ch,cu);
                       curSpikes(s) = spikeCount(curR,s,ch,cu);
                       curSpikesItime(s) = spikeInitialTime(curR,s,ch,cu);
                       curSpikesMtime(s) = spikeMaxTime(curR,s,ch,cu);
                       curSpikesDuration(s) = spikeDeconDur(curR,s,ch,cu);
                       curSpikesPeaks(s) = spikePeaks(curR,s,ch,cu);
%                        nT20(s) = spikeN20(curR,s,ch,cu);
                       cnt = cnt+1;
                   end
               end
               if(cnt==1)
                   tN = tN+1;
                   endTrendsOnset(tN,numel(files2Load))=NaN; % expand array
                   endTrendsDur(tN,numel(files2Load))=NaN;
                   endTrendsVol(tN,numel(files2Load))=NaN;
                   endTrendsMax(tN,numel(files2Load))=NaN;
                   endTrendsOnset(tN,:)=curOnset; % save valid trendline
                   endTrendsDur(tN,:)=curDur;
                   endTrendsVol(tN,:)=curVol;
                   endTrendsMax(tN,:)=curMax;
                   endTrendsRoiCh(tN,1) = curR;
                   endTrendsRoiCh(tN,2) = chans(ch);
                   endTrendsRoiCh(tN,3) = allCurrs(cu);
%                    endTrendsN20(tN,numel(files2Load))=NaN;
                   endTrendsSpikes(tN,:)=curSpikes;
                   endTrendsSpikesItime(tN,:)=curSpikesItime;
                   endTrendsSpikesMtime(tN,:)=curSpikesMtime;
                   endTrendsSpikesDuration(tN,:)=curSpikesDuration;
                   endTrendsSpikesPeaks(tN,:)=curSpikesPeaks;
%                    endTrendsN20(tN,:) = nT20;
  
               end
            end
        end
    end
    
    
    
    %% plot traces of end ROIS
    if(figGen==1)

     for curR = 1:size(endMidTraces,1)
        for ch = 1:numel(chans)
            for cu = 1:numel(allCurrs)
                curFig = figure();
                title(strcat('End-ROI:',num2str(curR),{' '},'Chan:',num2str(chans(ch)),{' '},'Curr:',num2str(allCurrs(cu)),'\muA'))
                hold on
                legSess = [];
                for s = 1:numel(files2Load)
                   if(~isempty(endMidTraces{curR,s,ch,cu}))
                       legSess = [legSess, s];
                       pre = endPreTraces{curR,s,ch,cu};
                       midPost = endMidTraces{curR,s,ch,cu};
                       dff = [pre,midPost];
                       x = 33*((1:numel(dff))-numel(pre));
                       plot(x,dff,'linewidth',2)
                   end
                end
                xlim([-1000 4000])
                xline(0,'--');
                xlabel('ms')
                ylabel('\DeltaF/F')
                if(~isempty(legSess))
                   legend(num2str(legSess'))
                   % Save figure
                    txt = strcat('End-ROI',num2str(curR),'Chan',num2str(chans(ch)),'Curr',num2str(allCurrs(cu)));
                    combineFileName = strcat(curFld,txt,'.tif');
                    % saveas(curFig,combineFileName);

                    close(curFig)
                else
                    close(curFig)
                end
            end
        end
     end
    end
     
     
     
   
%% Plot grouped traces
if(figGen==1)
    for ch = 1:numel(chans)
        for cu = 1:numel(allCurrs)
            curFig = figure();
            hold on
            aveDFFPre = NaN(1,1);
            aveDFFMid = NaN(1,1);
            aDC = 0;
            for curR = 1:size(endMidTraces,1)
                for s = 1:numel(files2Load)
                   if(~isempty(endMidTraces{curR,s,ch,cu}))
                       pre = endPreTraces{curR,s,ch,cu};
                       midPost = endMidTraces{curR,s,ch,cu};
                       aDC = aDC+1;
                       aveDFFPre(aDC,numel(pre)) = NaN;
                       aveDFFMid(aDC,numel(midPost)) = NaN;
                       aveDFFPre(aDC,1:numel(pre)) = pre;
                       aveDFFMid(aDC,1:numel(midPost)) = midPost;
                       dff = [pre,midPost];
                       x = 33*((1:numel(dff))-numel(pre));
                       plot(x,dff,'color',[0 0 0 0.15],'linewidth',2)
                       
                   end
                end
            end
            
           ia = endTrendsRoiCh(:,3) == allCurrs(cu);
           ib = endTrendsRoiCh(:,2) == chans(ch);
           ic = ia & ib;
           
           aveOnset = reshape(endTrendsOnset(ic,:),1,[]);
           aveDur = reshape(endTrendsDur(ic,:),1,[]);
            
            
            title(strcat('End ROIs Grouped Chan:',num2str(chans(ch)),{' '},'Curr:',num2str(allCurrs(cu)),'\muA'))

            anstr = strcat('Onset',{' '}, num2str(mean(aveOnset,'omitnan')),'ms',...
            {'   '},'Duration',{' '}, num2str(mean(aveDur,'omitnan')),'ms');
%             {'   '},'Intensity',{' '}, num2str(mean(aveMax,'omitnan')),...
%             {'   '},'Volume',{' '}, num2str(mean(aveVol,'omitnan')));
            dim = [.2 .62 .3 .3];
            annotation('textbox',dim,'String',anstr{1},'FitBoxToText','on');

            meanPre = mean(aveDFFPre,1,'omitnan');
            meanMid = mean(aveDFFMid,1,'omitnan');
            dff = [meanPre,meanMid];
            x = 33*((1:numel(dff))-numel(meanPre));
            plot(x,dff,'color',[0 0 0],'linewidth',3)
            
            xlim([-1000 4000])
            xline(0,'--');
            xlabel('ms')
            ylabel('\DeltaF/F')
            

            

            % Save figure
            txt = strcat('End-Grouped-Chan',num2str(chans(ch)),'Curr',num2str(allCurrs(cu)));
            combineFileName = strcat(curFld,txt,'.tif');
            % saveas(curFig,combineFileName);

% %             % zoom in on rising edge
% %             xlim([-100 400])
% %             combineFileName = strcat(curFld,'RisingEdge',txt,'.tif');
% %             saveas(curFig,combineFileName);
% %             combineFileName = strcat(curFld,'RisingEdge',txt,'.fig');
% %             saveas(curFig,combineFileName);
            close(curFig)
        end
    end
    
    
    
    
    
    
    
    
       %% Plot grouped traces
   if(figGen==1)
    for cu = 1:numel(allCurrs)
        curFig = figure();
            hold on
            aveDFFPre = NaN(1,1);
            aveDFFMid = NaN(1,1);
            aDC = 0;
        for ch = 1:numel(chans)
            
            for curR = 1:size(endMidTraces,1)
                for s = 1:numel(files2Load)
                   if(~isempty(endMidTraces{curR,s,ch,cu}))
                       pre = endPreTraces{curR,s,ch,cu};
                       midPost = endMidTraces{curR,s,ch,cu};
                       aDC = aDC+1;
                       aveDFFPre(aDC,numel(pre)) = NaN;
                       aveDFFMid(aDC,numel(midPost)) = NaN;
                       aveDFFPre(aDC,1:numel(pre)) = pre;
                       aveDFFMid(aDC,1:numel(midPost)) = midPost;
                       dff = [pre,midPost];
                       x = 33*((1:numel(dff))-numel(pre));
                       plot(x,dff,'color',[0 0 0 0.05],'linewidth',2)
                       
                   end
                end
            end
        end
            
       ia = endTrendsRoiCh(:,3) == allCurrs(cu);
       ib = endTrendsRoiCh(:,2) == chans(ch);
       ic = ia & ib;

       aveOnset = reshape(endTrendsOnset(ic,:),1,[]);
       aveDur = reshape(endTrendsDur(ic,:),1,[]);


        title(strcat('End ROIs Grouped Chan:',num2str(chans(ch)),{' '},'Curr:',num2str(allCurrs(cu)),'\muA'))

        anstr = strcat('Onset',{' '}, num2str(mean(aveOnset,'omitnan')),'ms',...
        {'   '},'Duration',{' '}, num2str(mean(aveDur,'omitnan')),'ms');
        dim = [.2 .62 .3 .3];
        annotation('textbox',dim,'String',anstr{1},'FitBoxToText','on');

        meanPre = mean(aveDFFPre,1,'omitnan');
        meanMid = mean(aveDFFMid,1,'omitnan');
        dff = [meanPre,meanMid];
        x = 33*((1:numel(dff))-numel(meanPre));
        plot(x,dff,'color',[0 0 0],'linewidth',3)

        xlim([-1000 4000])
        xline(0,'--');
        xlabel('ms')
        ylabel('\DeltaF/F')




        % Save figure
        txt = strcat('End-Grouped-Curr',num2str(allCurrs(cu)));
        combineFileName = strcat(curFld,txt,'.tif');
        % saveas(curFig,combineFileName);


        close(curFig)
    end
    end
end
    
    
    
    
    
    %% Plot all metric trends

    overallTrends = zeros(8,2);
    for t = 1:10
        if(t==1)
            % plot max DFF for all traces analyzed
            metricTrend = allTrendsMax;
            ytext = '\DeltaF/F';
            titleTxt = 'Maximum \DeltaF/F over Training';
            svTxt = 'MaxDFF';
            edges = 0:0.5:5;
            thresh = 0.2;
        elseif(t==2)
            % plot onset time for all traces analyzed
            metricTrend = allTrendsOnset;
            metricTrend(metricTrend>700)=NaN;
            ytext = 'ms';
            titleTxt = 'Spike Onset Time over Training';
            svTxt = 'Onset';
            edges = 0:50:1000;
            thresh = 15;
        elseif(t==3)
            % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsDur;
            ytext = 'ms';
            titleTxt = 'Spike Duration over Training';
            svTxt = 'Duration';
            edges = 0:200:3000;
            
        elseif(t==4)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsN20;
            ytext = 'ms';
            titleTxt = 'Normalized Onset over Training';
            svTxt = 'NormalizedOnset';
            edges = 0:10:500;
        elseif(t==5)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsSpikes;
            ytext = 'Spikes';
            titleTxt = 'Deconvoluted Spiking over Training';
            svTxt = 'SpikeSum';
        elseif(t==6)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsSpikesItime;
            metricTrend(metricTrend>600)=NaN;
            ytext = 'ms';
            titleTxt = 'Deconvoluted Spiking Onset over Training';
            svTxt = 'SpikeOnset';
        elseif(t==7)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsSpikesMtime;
            ytext = 'ms';
            titleTxt = 'Deconvoluted Spiking Max Time over Training';
            svTxt = 'SpikeMaxTime';
            
        elseif(t==8)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsSpikesDuration;
            ytext = 'ms';
            titleTxt = 'Deconvoluted Spiking Duration over Training';
            svTxt = 'SpikeDeconDur';
        elseif(t==9)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsSpikesPeaks;
            ytext = '# of Peaks';
            titleTxt = 'Deconvoluted number of midstim peaks over Training';
            svTxt = 'SpikePeaks';            
        else
            % plot volume under stimulation curve for all traces analyzed
            metricTrend = allTrendsVol;
            ytext = '\DeltaF/F x ms';
            titleTxt = 'Spike Volume over Training';
            svTxt = 'Volume';
            edges = 0:10:200;
        end
        fitM = zeros(size(metricTrend,1),1);
        allP = zeros(size(metricTrend,1),2);
        fGain = zeros(size(metricTrend,1),1);
        fMax = zeros(size(metricTrend,1),1);
        bMax = zeros(size(metricTrend,1),1);

        % curFig = figure();
        % hold on
        % aX = [];
        % aY = [];

        for i = 1:size(metricTrend,1) % Add check if metricTrend is all NaNs
            % aX = [aX, daysTrained];
            % aY = [aY, metricTrend(i,:)];
            % plot(daysTrained,metricTrend(i,:),'color',[0 0 0 0.15],'linewidth',1)
            % scatter1 = scatter(daysTrained,metricTrend(i,:),'k');
            % scatter1.MarkerFaceAlpha = .15;
            % scatter1.MarkerEdgeAlpha = .15;

%             d = 1:numel(daysTrained);
            d = daysTrained;
            m = metricTrend(i,:);

            if(numel(m)==sum(isnan(m)))
                continue
            end
            d(isnan(m)) = [];
            m(isnan(m)) = [];
            p = polyfit(d, m, 1);
            fitM(i) = p(1);
            allP(i,:) = p;
            fGain(i) = m(end)-m(1);
            fMax(i) = m(end);
            bMax(i) = m(1);
        end

        % aX(isnan(aY))=[];
        % aY(isnan(aY))=[];
        % p = polyfit(aX, aY, 1);
        % overallTrends(t,:) = p;
        % px = [min(aX) max(aX)];
        % py = polyval(p, px);
        % plot(px,py,'color',[0 0 0],'linewidth',3)
        % xlim([min(aX)-1, max(aX)+1])
        % 
        % xlabel('Days of Training')
        % ylabel(ytext)
        % title(titleTxt)
        % combineFileName = strcat(curFld,'OverallTrendsLines',svTxt,'.fig');
        % % saveas(curFig,combineFileName);  
        % combineFileName = strcat(curFld,'OverallTrendsLines',svTxt,'.tif');
        % % saveas(curFig,combineFileName);
        % close(curFig)
        % 
        % 
%         curFig = figure();
%         histogram(fitM,15)
%         hold on
%         xline(mean(fitM),'-r',{'Mean'},'linewidth',2)
% %         xline(mean(fitM) + std(fitM),'--',{strcat('Thresh:',num2str(mean(fitM) + std(fitM)))},'linewidth',2)
%         xline(mean(fitM) + std(fitM),'--r',{'+1 STD'},'linewidth',2)
% 
% 
%         xlabel('Slope')
%         ylabel('Count')
%         title(titleTxt)
%         combineFileName = strcat(curFld,'OverallTrends',svTxt,'.fig');
%         % saveas(curFig,combineFileName);  
%         combineFileName = strcat(curFld,'OverallTrends',svTxt,'.tif');
%         % saveas(curFig,combineFileName);
%         % close(curFig)
        
        
        
        % 
        % % Plot distance of ROIs relative to metric being examined
        % curFig = figure();
        % metricTrendROI = allTrendsRoiCh;
        % cents = NaN(size(metricTrendROI,1),2);
        % for r = 1:size(metricTrendROI,1)
        %     curROIs = ROIs2Compare(metricTrendROI(r),:);
        %     datSet = find(curROIs,1);
        %     curData = AllData{datSet};
        %     regions = curData.regions;
        %     curPx = regions(curROIs(datSet)).PixelList;
        %     cents(r,:) = mean(curPx,1);
        % end
        % curDists = pdist2(mean(chanPositions,1),cents);
        % scatter(curDists,fitM)
        % xlabel('Neuron Distance (\mum)')
        % ylabel(strcat(ytext,{' '},'Slope'))
        % title(strcat(titleTxt,{' '},'Activation Distance'))
        % combineFileName = strcat(curFld,'DistanceTrends',svTxt,'.fig');
        % % saveas(curFig,combineFileName);  
        % combineFileName = strcat(curFld,'DistanceTrends',svTxt,'.tif');
        % % saveas(curFig,combineFileName);
        % close(curFig)
        
        
        
        
        
        
        % split trends into ascending, descending, and non-monotonic
        nonmon = true(size(metricTrend,1),1);
        ascend = fitM > (mean(fitM) + std(fitM));
        decend = fitM < (mean(fitM) - std(fitM));
%         ascend = fitM > thresh;
%         decend = fitM < -thresh;
        nonmon(ascend==1) = 0;
        nonmon(decend==1) = 0;

        % Split groups into those with the highest and lowest flourescence
        % gain over training
        nonmonGain = true(size(metricTrend,1),1);
        ascendGain = fGain > (mean(fGain) + std(fGain));
        decendGain = fGain < (mean(fGain) - std(fGain));
        nonmonGain(ascendGain==1) = 0;
        nonmonGain(decendGain==1) = 0;
        
        % Split groups into those with the highest and lowest flourescence
        % at the end of training
        nonmonEnd = true(size(metricTrend,1),1);
        ascendEnd = fMax > (mean(fMax) + std(fMax));
        decendEnd = fMax < (mean(fMax) - std(fMax));
        nonmonEnd(ascendEnd==1) = 0;
        nonmonEnd(decendEnd==1) = 0;

        if(t==1)
            % subset = ascend; % slope trend
            % subset2 = nonmon | decend; 

            fMaxAll = [fMaxAll; fMax];
            bMaxAll = [bMaxAll; bMax];
            % subset = ascendGain; % gain trend
            % subset2 = nonmonGain | decendGain; 

            subset = ascendEnd; % max end trend
            subset2 = nonmonEnd | decendEnd; 
        end
        
        ascendTrend = metricTrend(ascend,:);
        descendTrend = metricTrend(decend,:);
        nomonoTrend = metricTrend(nonmon,:);
        
% 
%         curFig = figure('Name', titleTxt);
%         subplot(1,3,1) % Ascending ---------
%         hold on
%         for i = 1:size(ascendTrend,1)
%             plot(daysTrained,ascendTrend(i,:),'color',[1 0 0 0.15],'linewidth',2)
%             scatter1 = scatter(daysTrained,ascendTrend(i,:),[],[1 0 0]);
%             scatter1.MarkerFaceAlpha = .15;
%             scatter1.MarkerEdgeAlpha = .15;
%         end
% %         pAve = mean(allP(ascend,:),1);
% %         px = [min(daysTrained) max(daysTrained)];
% %         py = polyval(pAve, px);
% %         plot(px,py,'color',[1 0 0],'linewidth',3)
%         aveAll = mean(ascendTrend,1,'omitnan');
%         plot(daysTrained,aveAll,'color',[1 0 0],'linewidth',3)
%         xlim([min(daysTrained)-1, max(daysTrained)+1])
%         xlabel('Days of Training')
%         ylabel(ytext)
%         title('Trending Ascending')
% 
% 
%         subplot(1,3,2) % Descending -----------
%         hold on
%         for i = 1:size(descendTrend,1)            
%             plot(daysTrained,descendTrend(i,:),'color',[0 0 1 0.15],'linewidth',2)
%             scatter1 = scatter(daysTrained,descendTrend(i,:),[],[0 0 1]);
%             scatter1.MarkerFaceAlpha = .15;
%             scatter1.MarkerEdgeAlpha = .15;
%         end
%         aveAll = mean(descendTrend,1,'omitnan');
%         plot(daysTrained,aveAll,'color',[0 0 1],'linewidth',3)
% %         pAve = mean(allP(decend,:),1);
% %         px = [min(daysTrained) max(daysTrained)];
% %         py = polyval(pAve, px);
% %         plot(px,py,'color',[0 0 1],'linewidth',3)
%         xlim([min(daysTrained)-1, max(daysTrained)+1])
%         xlabel('Days of Training')
%         ylabel(ytext)
%         title('Trending Descending')
% 
%         subplot(1,3,3) % Non-Monotonic -----------
%         hold on
%         for i = 1:size(nomonoTrend,1)
%             plot(daysTrained,nomonoTrend(i,:),'color',[0 1 0 0.15],'linewidth',2)
%             scatter1 = scatter(daysTrained,nomonoTrend(i,:),[],[0 1 0]);
%             scatter1.MarkerFaceAlpha = .15;
%             scatter1.MarkerEdgeAlpha = .15;
%         end
%         aveAll = mean(nomonoTrend,1,'omitnan');
%         plot(daysTrained,aveAll,'color',[0.0549    0.4314    0.0314],'linewidth',3)
% %         pAve = mean(allP(nonmon,:),1);
% %         px = [min(daysTrained) max(daysTrained)];
% %         py = polyval(pAve, px);
% %         plot(px,py,'color',[0.0549    0.4314    0.0314],'linewidth',3)
%         xlim([min(daysTrained)-1, max(daysTrained)+1])
%         xlabel('Days of Training')
%         ylabel(ytext)
%         title('Trending Non-Monotonic')
%         combineFileName = strcat(curFld,'DividedTrends',svTxt,'.fig');
%         % saveas(curFig,combineFileName);  
%         combineFileName = strcat(curFld,'DividedTrends',svTxt,'.tif');
%         % saveas(curFig,combineFileName);
%         close(curFig)

    end

    %% Plot spatial distribution of subset and non-subset populations   
    if(figGen==1)

    metricTrendInvROI = allTrendsRoiCh(subset2,:);
    metricTrendROI = allTrendsRoiCh(subset,:);
    
    cents = NaN(size(metricTrendROI,1),2);
    for r = 1:size(metricTrendROI,1)
        curROIs = ROIs2Compare(metricTrendROI(r),:);
        datSet = find(curROIs,1);
        curData = AllData{datSet};
        regions = curData.regions;
        curPx = regions(curROIs(datSet)).PixelList;
        cents(r,:) = mean(curPx,1);
    end
%     meanCent = mean(cents,1);
    curDistsA = pdist2(mean(chanPositions,1),cents);
    distanceSub{animal} = curDistsA;
    
    cents = NaN(size(metricTrendInvROI,1),2);
    for r = 1:size(metricTrendInvROI,1)
        curROIs = ROIs2Compare(metricTrendInvROI(r),:);
        datSet = find(curROIs,1);
        curData = AllData{datSet};
        regions = curData.regions;
        curPx = regions(curROIs(datSet)).PixelList;
        cents(r,:) = mean(curPx,1);
    end
    meanCent = mean(cents,1);
    curDists = pdist2(mean(chanPositions,1),cents);
    distanceNon{animal} = curDists;
    
    % Plot difference between two populations
    curFig = figure();
    curDistVals = curDists';
    curDistG = ones(size(curDists',1),1)*2;
    curDistVals = [curDistVals;curDistsA'];
    curDistG = [curDistG;ones(size(curDistsA',1),1)*1];
    boxplot(curDistVals,curDistG)
    title('Comparision of Subgroup Vs Main Population Activation Distance')
    ylabel('Activation Distance (\mum)')
    combineFileName = strcat(curFld,'ActivatedDistance',svTxt,'.fig');
    % saveas(curFig,combineFileName);  
    combineFileName = strcat(curFld,'ActivatedDistance',svTxt,'.tif');
    % saveas(curFig,combineFileName);
    close(curFig)
    
    
    
    allRasters = cell(numel(files2Load),1);
    for s = 1:numel(files2Load)
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
        allRasters{s} = sessData.allRasterH;
    end
    
%     % determine which ROIs are active for each electrode, for each scenario
%     % calculate distance
%     curDistsA_elec = [];
%     cents = NaN(size(metricTrendROI,1),2);
%     for r = 1:size(metricTrendROI,1)
%         curROIs = ROIs2Compare(metricTrendROI(r),:);
%         datSet = find(curROIs,1);
%         curData = AllData{datSet};
%         regions = curData.regions;
%         curPx = regions(curROIs(datSet)).PixelList;
%         curCent = mean(curPx,1);
%         
%         % check what channels it was activated by, measure distance from
%         % each
%         sessRasterCnts = zeros(3,size(ROIs2Compare,2));
%         for s = 1:size(ROIs2Compare,2)
%             if(curROIs(s)~=0)
%                 allRasterH = allRasters{s};
%                 sessRasterCnts(:,s) = sum(allRasterH(:,:,curROIs(s)),2);
%             end
%         end
%         sessRasterBin = sum(sessRasterCnts,2)>1;
%         for sr = 1:3
%             if(sessRasterBin(sr)==1)
%                 % ROI was active in multiple sessions for given channel
%                 curDistsA_elec = [curDistsA_elec, pdist2(chanPositions(sr,:),curCent)];
%             end
%         end
%     end
%     distanceSub_ver2{animal} = curDistsA_elec;
    
%     
%         % determine which ROIs are active for each electrode, for each scenario
%     % calculate distance
%     curDists_elec = [];
%     cents = NaN(size(metricTrendInvROI,1),2);
%     for r = 1:size(metricTrendInvROI,1)
%         curROIs = ROIs2Compare(metricTrendInvROI(r),:);
%         datSet = find(curROIs,1);
%         curData = AllData{datSet};
%         regions = curData.regions;
%         curPx = regions(curROIs(datSet)).PixelList;
%         curCent = mean(curPx,1);
%         
%         % check what channels it was activated by, measure distance from
%         % each
%         sessRasterCnts = zeros(3,size(ROIs2Compare,2));
%         for s = 1:size(ROIs2Compare,2)
%             if(curROIs(s)~=0)
%                 allRasterH = allRasters{s};
%                 sessRasterCnts(:,s) = sum(allRasterH(:,:,curROIs(s)),2)>0;
%             end
%         end
%         sessRasterBin = sum(sessRasterCnts,2)>1;
%         for sr = 1:3
%             if(sessRasterBin(sr)==1)
%                 % ROI was active in multiple sessions for given channel
%                 curDists_elec = [curDists_elec, pdist2(chanPositions(sr,:),curCent)];
%             end
%         end
%     end
%     distanceNon_ver2{animal} = curDists_elec;
%     
%     % Plot difference between two populations
%     curFig = figure();
%     curDistVals = curDists_elec';
%     curDistG = ones(size(curDists_elec',1),1)*2;
%     curDistVals = [curDistVals;curDistsA_elec'];
%     curDistG = [curDistG;ones(size(curDistsA_elec',1),1)*1];
%     boxplot(curDistVals,curDistG)
%     title('Comparision of Subgroup Vs Main Population Activation Distance')
%     ylabel('Activation Distance (\mum)')
%     combineFileName = strcat(curFld,'ActivatedDistance_Ver2',svTxt,'.fig');
%     saveas(curFig,combineFileName);  
%     combineFileName = strcat(curFld,'ActivatedDistance_Ver2',svTxt,'.tif');
%     saveas(curFig,combineFileName);
%     close(curFig)
%     
    
    end
    
    
      %% Plot all metric trends
  if(figGen==1)

    overallTrends = zeros(8,2);
    for t = 1:8
        if(t==1)
            % plot max DFF for all traces analyzed
            metricTrend = allTrendsMax;
            ytext = '\DeltaF/F';
            titleTxt = 'Maximum \DeltaF/F over Training';
            svTxt = 'MaxDFF';
            edges = 0:0.5:5;
            thresh = 0.2;
        elseif(t==2)
            % plot onset time for all traces analyzed
            metricTrend = allTrendsOnset;
            metricTrend(metricTrend>700)=NaN;
            ytext = 'ms';
            titleTxt = 'Spike Onset Time over Training';
            svTxt = 'Onset';
            edges = 0:50:1000;
            thresh = 15;
        elseif(t==3)
            % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsDur;
            ytext = 'ms';
            titleTxt = 'Spike Duration over Training';
            svTxt = 'Duration';
            edges = 0:200:3000;
            
        elseif(t==4)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsN20;
            ytext = 'ms';
            titleTxt = 'Normalized Onset over Training';
            svTxt = 'NormalizedOnset';
            edges = 0:10:500;
        elseif(t==5)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsSpikes;
            ytext = 'Spikes';
            titleTxt = 'Deconvoluted Spiking over Training';
            svTxt = 'SpikeSum';
        elseif(t==6)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsSpikesItime;
            metricTrend(metricTrend>600)=NaN;
            ytext = 'ms';
            titleTxt = 'Deconvoluted Spiking Onset over Training';
            svTxt = 'SpikeOnset';
        elseif(t==7)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsSpikesMtime;
            ytext = 'ms';
            titleTxt = 'Deconvoluted Spiking Max Time over Training';
            svTxt = 'SpikeMaxTime';
            
        elseif(t==8)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsSpikesDuration;
            ytext = 'ms';
            titleTxt = 'Deconvoluted Spiking Duration over Training';
            svTxt = 'SpikeDeconDur';
        elseif(t==9)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsSpikesPeaks;
            ytext = '# of Peaks';
            titleTxt = 'Deconvoluted number of midstim peaks over Training';
            svTxt = 'SpikePeaks'; 
            
        else
            % plot volume under stimulation curve for all traces analyzed
            metricTrend = allTrendsVol;
            ytext = '\DeltaF/F x ms';
            titleTxt = 'Spike Volume over Training';
            svTxt = 'Volume';
            edges = 0:10:200;
        end
        fitM = zeros(size(metricTrend,1),1);

        curFig = figure();
        hold on
        aX = [];
        aY = [];
        for i = 1:size(metricTrend,1)
            aX = [aX, daysTrained];
            aY = [aY, metricTrend(i,:)];
            plot(daysTrained,metricTrend(i,:),'color',[0 0 0 0.15],'linewidth',1)
            scatter1 = scatter(daysTrained,metricTrend(i,:),'k');
            scatter1.MarkerFaceAlpha = .15;
            scatter1.MarkerEdgeAlpha = .15;

            d = 1:numel(daysTrained);
            m = metricTrend(i,:);
            d(isnan(m)) = [];
            m(isnan(m)) = [];
            p = polyfit(d, m, 1);
            fitM(i) = p(1);
        end

        aX(isnan(aY))=[];
        aY(isnan(aY))=[];
        p = polyfit(aX, aY, 1);
        overallTrends(t,:) = p;
        px = [min(aX) max(aX)];
        py = polyval(p, px);
        plot(px,py,'color',[0 0 0],'linewidth',3)
        xlim([min(aX)-1, max(aX)+1])

        xlabel('Days of Training')
        ylabel(ytext)
        title(titleTxt)
        combineFileName = strcat(curFld,'OverallTrends',svTxt,'.fig');
        % saveas(curFig,combineFileName);  
        combineFileName = strcat(curFld,'OverallTrends',svTxt,'.tif');
        % saveas(curFig,combineFileName);
        close(curFig)

        % split trends into ascending, descending, and non-monotonic
        nonmon = true(size(metricTrend,1),1);
        ascend = fitM > thresh;
        decend = fitM < -thresh;
        nonmon(ascend==1) = 0;
        nonmon(decend==1) = 0;
        ascendTrend = metricTrend(ascend,:);
        descendTrend = metricTrend(decend,:);
        nomonoTrend = metricTrend(nonmon,:);
        

        curFig = figure('Name', titleTxt);
        subplot(1,3,1) % Ascending ---------
        hold on
        for i = 1:size(ascendTrend,1)
            plot(daysTrained,ascendTrend(i,:),'color',[1 0 0 0.15],'linewidth',2)
            scatter1 = scatter(daysTrained,ascendTrend(i,:),[],[1 0 0]);
            scatter1.MarkerFaceAlpha = .15;
            scatter1.MarkerEdgeAlpha = .15;
        end
        aveAll = mean(ascendTrend,1,'omitnan');
        plot(daysTrained,aveAll,'color',[1 0 0],'linewidth',3)
        xlim([min(daysTrained)-1, max(daysTrained)+1])
        xlabel('Days of Training')
        ylabel(ytext)
        title('Trending Ascending')
        
        
        subplot(1,3,2) % Descending -----------
        hold on
        for i = 1:size(descendTrend,1)            
            plot(daysTrained,descendTrend(i,:),'color',[0 0 1 0.15],'linewidth',2)
            scatter1 = scatter(daysTrained,descendTrend(i,:),[],[0 0 1]);
            scatter1.MarkerFaceAlpha = .15;
            scatter1.MarkerEdgeAlpha = .15;
        end
        aveAll = mean(descendTrend,1,'omitnan');
        plot(daysTrained,aveAll,'color',[0 0 1],'linewidth',3)
        xlim([min(daysTrained)-1, max(daysTrained)+1])
        xlabel('Days of Training')
        ylabel(ytext)
        title('Trending Descending')
        
        subplot(1,3,3) % Non-Monotonic -----------
        hold on
        for i = 1:size(nomonoTrend,1)
            plot(daysTrained,nomonoTrend(i,:),'color',[0 1 0 0.15],'linewidth',2)
            scatter1 = scatter(daysTrained,nomonoTrend(i,:),[],[0 1 0]);
            scatter1.MarkerFaceAlpha = .15;
            scatter1.MarkerEdgeAlpha = .15;
        end
        aveAll = mean(nomonoTrend,1,'omitnan');
        plot(daysTrained,aveAll,'color',[0.0549    0.4314    0.0314],'linewidth',3)
        xlim([min(daysTrained)-1, max(daysTrained)+1])
        xlabel('Days of Training')
        ylabel(ytext)
        title('Trending Non-Monotonic')
        combineFileName = strcat(curFld,'DividedTrends',svTxt,'.fig');
        % saveas(curFig,combineFileName);  
        combineFileName = strcat(curFld,'DividedTrends',svTxt,'.tif');
        % saveas(curFig,combineFileName);
        close(curFig)
    end

  end

    %% Perform comparisions at channels over sessions/weeks at individual
    if(figGen==1)
    % currents and plot
    for curR = 1:size(ROIs2Compare,1)
        for ch = 1:numel(chans)
            for cu = 1:numel(allCurrs)
               curFig = figure();
               title(strcat('ROI:',num2str(curR),{' '},'Chan:',num2str(chans(ch)),{' '},'Curr:',num2str(allCurrs(cu)),'\muA'))
               hold on
               legSess = [];
               for s = 1:numel(files2Load)
                   if(~isempty(midTraces{curR,s,ch,cu}))
                       legSess = [legSess, s];
                       pre = preTraces{curR,s,ch,cu};
                       midPost = midTraces{curR,s,ch,cu};
                       dff = [pre,midPost];
                       x = 33*((1:numel(dff))-numel(pre));
                       plot(x,dff,'linewidth',2)
                   end
               end
               xlim([-1000 4000])
               xline(0,'--');
               xlabel('ms')
               ylabel('\DeltaF/F')
               if(~isempty(legSess))
%                    legSess = weeksTrained(legSess);
                   legend(num2str(legSess'))
               end
               if(numel(legSess)>minSetCnt)
                   % Save figure
                   txt = strcat('ROI',num2str(curR),'Chan',num2str(chans(ch)),'Curr',num2str(allCurrs(cu)));
    %                combineFileName = strcat(curFld,txt,'.fig');
    %                 saveas(curFig,combineFileName);  
                    combineFileName = strcat(curFld,txt,'.tif');
                    % saveas(curFig,combineFileName);
                    
% %                     % zoom in on rising edge
% %                     xlim([-100 400])
% %                     combineFileName = strcat(curFld,'RisingEdge',txt,'.tif');
% %                     saveas(curFig,combineFileName);
                    close(curFig)
               else
                   % dont save, just close figure
                   close(curFig)
               end
           end
        end
    end
    end
    
    

    
% %     %% Plot Spiking Trends traces over training
% %     for curR = 1:size(ROIs2Compare,1)
% %         for ch = 1:numel(chans)
% %             for cu = 1:numel(allCurrs)
% %                curFig = figure();
% %                hold on
% %                sessCnt = 0;
% %                
% %                for s = 1:numel(files2Load)
% %                    if(~isempty(spikeTraces{curR,s,ch,cu}))
% %                        sessCnt = sessCnt+1;
% %                        subplot(numel(files2Load),1,s)
% %                        pre = preSpikeTraces{curR,s,ch,cu};
% %                        midPost = spikeTraces{curR,s,ch,cu};
% %                        dff = [pre,midPost];
% %                        x = 33*((1:numel(dff))-numel(pre));
% %                        bar(x,dff)
% %                        if(sessCnt==1)
% %                           title(strcat('ROI:',num2str(curR),{' '},'Chan:',num2str(chans(ch)),{' '},'Curr:',num2str(allCurrs(cu)),'\muA'))
% %                        end
% %                        xlim([-500 2000])
% %                        ylim([0 1])
% %                        xline(0,'--');
% %                        xlabel('ms')
% %                        ylabel('Spiking Probability')
% %                    end
% %                end
% % 
% %                if(sessCnt>minSetCnt)
% %                    % Save figure
% %                    txt = strcat('Spiking-ROI',num2str(curR),'Chan',num2str(chans(ch)),'Curr',num2str(allCurrs(cu)));
% %     %                combineFileName = strcat(curFld,txt,'.fig');
% %     %                 saveas(curFig,combineFileName);  
% %                     combineFileName = strcat(curFld,txt,'.tif');
% %                     saveas(curFig,combineFileName);                   
% %                end
% %                close(curFig)
% %            end
% %         end
% %     end
% %     
% %     
    
    
    
    %%
    if(figGen==1)

	for s = 1:numel(files2Load)
        curFig = figure();
        hold on
        runningPre = [];
        runningMid = [];
        rN = 0;
        for curR = 1:size(ROIs2Compare,1)
            for ch = 1:numel(chans)
                for cu = 1:numel(allCurrs)
                   
                    sessVld = false(numel(files2Load),1);
                    for s2 = 1:numel(files2Load)
                        sessVld(s2) = ~isempty(spikeTraces{curR,s2,ch,cu});
                    end
                    if(sum(sessVld)>minSetCnt)
                       if(~isempty(spikeTraces{curR,s,ch,cu}))
                           pre = preSpikeTraces{curR,s,ch,cu};
                           midPost = spikeTraces{curR,s,ch,cu};
                           dff = [pre,midPost];
                           x = 33*((1:numel(dff))-numel(pre));
                           plot(x,smoothdata(dff,"gaussian",10),'color',[0 0 0 0.1])
                           
                           rN = rN+1;
                           runningPre(rN,size(pre,2))=NaN;
                           runningPre(rN,1:size(pre,2))=flip(pre);
                           runningMid(rN,size(midPost,2))=NaN;
                           runningMid(rN,1:size(midPost,2))=midPost;
                       end
                   end
                end
           end
        end
        rPM = flip(mean(runningPre,1,'omitnan'));
        dff = [rPM,mean(runningMid,1,'omitnan')];
        x = 33*((1:numel(dff))-numel(rPM));
        plot(x,smoothdata(dff,"gaussian",10),'color',[0 0 0],'linewidth',3)
        
        title(strcat('Session:',num2str(s)))
        xlim([-500 2000])
        ylim([0 1])
        xline(0,'--');
        xlabel('ms')
        ylabel('Spiking Probability')
        
        % Save figure
        txt = strcat('Spiking-SessionAve_Session',num2str(s));
        combineFileName = strcat(curFld,txt,'.tif');
        % saveas(curFig,combineFileName);   
        close(curFig)
	end               
    end
    
    
    % Run analysis and see if I can segment out double peaks
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
% %     %% 
% %      % currents and plot
% %     for curR = 1:size(ROIs2Compare,1)
% %         for ch = 1:numel(chans)
% %             for cu = 1:numel(allCurrs)
% %                curFig = figure();
% %                legSess = [];
% %                for s = 1:numel(files2Load)
% %                    subplot(numel(files2Load),1,s)
% %                    hold on
% %                    if(~isempty(midTraces{curR,s,ch,cu}))
% %                        legSess = [legSess, s];
% %                        pre = preTraces{curR,s,ch,cu};
% %                        midPost = midTraces{curR,s,ch,cu};
% %                        dff = [pre,midPost];
% %                        x = 33*((1:numel(dff))-numel(pre));
% %                        plot(x,dff,'linewidth',2)
% %                        xlim([-1000 4000])
% %                        xline(0,'--');
% %                        xlabel('ms')
% %                        ylabel('\DeltaF/F')
% %                        
% %                        yyaxis right
% %                        midSpiking = spikeTraces{curR,s,ch,cu};
% %                        conActiveS = midSpiking>0.5;
% %                        casInds = find(conActiveS);
% %                        if(~isempty(casInds))
% %                            for m = 1:numel(casInds)
% %                                xline(casInds(m)*33,'r');
% %                            end
% %                        end
% %                    end
% %                end
% %               
% %                if(numel(legSess)>minSetCnt)
% %                    % Save figure
% %                    txt = strcat('TraceAndSpike-ROI',num2str(curR),'Chan',num2str(chans(ch)),'Curr',num2str(allCurrs(cu)));
% %                    combineFileName = strcat(curFld,txt,'.tif');
% %                    saveas(curFig,combineFileName);
% %                end
% %                close(curFig)
% %            end
% %         end
% %     end
    
    
    
    
    
    
    
    
% % % % %     %% Perform comparisions at channels over sessions/weeks at individual
% % % % %     % currents and plot
% % % % %     for curR = 1:size(ROIs2Compare,1)
% % % % %         for ch = 1:numel(chans)
% % % % %             for cu = 1:numel(allCurrs)
% % % % %                curFig = figure();
% % % % %                title(strcat('Normalized---ROI:',num2str(curR),{' '},'Chan:',num2str(chans(ch)),{' '},'Curr:',num2str(allCurrs(cu)),'\muA'))
% % % % %                hold on
% % % % %                legSess = [];
% % % % %                for s = 1:numel(files2Load)
% % % % %                    if(~isempty(midTraces{curR,s,ch,cu}))
% % % % %                        legSess = [legSess, s];
% % % % %                        pre = preTraces{curR,s,ch,cu};
% % % % %                        midPost = midTraces{curR,s,ch,cu};
% % % % %                        dff = [pre,midPost];
% % % % %                        x = 33*((1:numel(dff))-numel(pre));
% % % % %                        dff = dff./max(dff);
% % % % %                        plot(x,dff,'linewidth',2)
% % % % %                    end
% % % % %                end
% % % % %                xlim([-1000 4000])
% % % % %                xline(0,'--');
% % % % %                xlabel('ms')
% % % % %                ylabel('\DeltaF/F')
% % % % %                if(~isempty(legSess))
% % % % %                    legend(num2str(legSess'))
% % % % %                end
% % % % %                if(numel(legSess)>minSetCnt)
% % % % %                    % Save figure
% % % % %                    txt = strcat('Normalized---ROI',num2str(curR),'Chan',num2str(chans(ch)),'Curr',num2str(allCurrs(cu)));
% % % % %                     combineFileName = strcat(curFld,txt,'.tif');
% % % % %                     saveas(curFig,combineFileName);
% % % % %                     
% % % % %                     % zoom in on rising edge
% % % % %                     xlim([-100 400])
% % % % %                     combineFileName = strcat(curFld,'RisingEdge',txt,'.tif');
% % % % %                     saveas(curFig,combineFileName);
% % % % %                     combineFileName = strcat(curFld,'RisingEdge',txt,'.fig');
% % % % %                     saveas(curFig,combineFileName);
% % % % %                     close(curFig)                    
% % % % %                else
% % % % %                    % dont save, just close figure
% % % % %                    close(curFig)
% % % % %                end
% % % % %            end
% % % % %         end
% % % % %     end
% % % % %     
% % % % %     
    
% 
%     % Plot all traces for all sessions and all currents
%     for curR = 1:size(ROIs2Compare,1)
%         for ch = 1:numel(chans)
%             txt = strcat('ICMS 98 - ROI',{' '},num2str(curR),{' '},'Channel',{' '},num2str(chans(ch)));
%            figure('Name',txt{1},'WindowState','maximized') 
%            for s = 1:numel(files2Load)
%                for cu = 1:numel(allCurrs)
%                    plotInd = (cu-1)*numel(files2Load) + s;
%                    h2 = subplot(numel(allCurrs),numel(files2Load),plotInd);
%                    if(~isempty(midTraces{curR,s,ch,cu}))
%                        pre = preTraces{curR,s,ch,cu};
%                        midPost = midTraces{curR,s,ch,cu};
%                        dff = [pre,midPost];
%                        x = 33*((1:numel(dff))-numel(pre));
%                        plot(x,dff,'linewidth',2)
%                        xlim([-1000 4000])
%                    end
%                    if(cu==1)
%                        % Report session
%                        title(strcat('Sess',num2str(s)))
%                    end
%                    if(s==1)
%                        % report current
%                        ylabel(strcat(num2str(allCurrs(cu)),'\muA'))
%                    end
%                end
%            end
%         end
%     end



allDFF{animal} = allTrendsMax;
allOnset{animal} = allTrendsOnset;
allDur{animal} = allTrendsDur;
allSpikeNum{animal} = allTrendsSpikes;
allSpikeOnset{animal} = allTrendsSpikesItime;
allVol{animal} = allTrendsVol;
allSubset{animal} = subset;
allDays{animal} = daysTrained;
allROI2Comp{animal} = ROIs2Compare;
allChanData{animal} = allTrendsRoiCh;




 %% Plot all metric trends - manual selected
%     subset = subset1;
%     subset2 = nonmon;

if(figGen==1)
    for t = 1:8
        if(t==1)
            % plot max DFF for all traces analyzed
            metricTrend = allTrendsMax;
            ytext = '\DeltaF/F';
            titleTxt = 'Maximum \DeltaF/F over Training';
            svTxt = 'MaxDFF';
            edges = 0:0.5:5;
        elseif(t==2)
            % plot onset time for all traces analyzed
            metricTrend = allTrendsOnset;
            ytext = 'ms';
            titleTxt = 'Spike Onset Time over Training';
            svTxt = 'Onset';
            edges = 0:50:1000;
        elseif(t==3)
            % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsDur;
            ytext = 'ms';
            titleTxt = 'Spike Duration over Training';
            svTxt = 'Duration';
            edges = 0:200:3000;
        elseif(t==4)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsN20;
            ytext = 'ms';
            titleTxt = 'Normalized Onset over Training';
            svTxt = 'NormalizedOnset';
            edges = 0:10:500;
        elseif(t==5)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsSpikes;
            ytext = 'Spikes';
            titleTxt = 'Deconvoluted Spiking over Training';
            svTxt = 'SpikeSum';
        elseif(t==6)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsSpikesItime;
            ytext = 'ms';
            titleTxt = 'Deconvoluted Spiking Onset over Training';
            svTxt = 'SpikeOnset';
        elseif(t==7)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsSpikesMtime;
            ytext = 'ms';
            titleTxt = 'Deconvoluted Spiking Max Time over Training';
            svTxt = 'SpikeMaxTime';
        elseif(t==8)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsSpikesDuration;
            ytext = 'ms';
            titleTxt = 'Deconvoluted Spiking Duration over Training';
            svTxt = 'SpikeDeconDur';
        elseif(t==9)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsSpikesPeaks;
            ytext = '# of Peaks';
            titleTxt = 'Deconvoluted number of midstim peaks over Training';
            svTxt = 'SpikePeaks'; 
        else
            % plot volume under stimulation curve for all traces analyzed
            metricTrend = allTrendsVol;
            ytext = '\DeltaF/F x ms';
            titleTxt = 'Spike Volume over Training';
            svTxt = 'Volume';
            edges = 0:10:200;
        end
        metricTrendInv = metricTrend(subset2,:);
        metricTrend = metricTrend(subset,:);

        curFig = figure();
        hold on
        for i = 1:size(metricTrend,1)
            plot(daysTrained,metricTrend(i,:),'color',[0 0 0 0.15],'linewidth',1)
            scatter1 = scatter(daysTrained,metricTrend(i,:),'k');
            scatter1.MarkerFaceAlpha = .15;
            scatter1.MarkerEdgeAlpha = .15;
        end
        
        
        
%         px = [min(daysTrained) max(daysTrained)];
%         py = polyval(overallTrends(t,:), px);
%         plot(px,py,'--','color','red','linewidth',2)
        
        aveAll = mean(metricTrendInv,1,'omitnan');
        plot(daysTrained,aveAll,'--','color','red','linewidth',2)
        
        aveAll = mean(metricTrend,1,'omitnan');
        plot(daysTrained,aveAll,'color',[0 0 0],'linewidth',3)
        xlim([min(daysTrained)-1, max(daysTrained)+1])

        
        
        xlabel('Days of Training')
        ylabel(ytext)
        title(titleTxt)
        
        combineFileName = strcat(curFld,'ManOverallTrends',svTxt,'.fig');
        % saveas(curFig,combineFileName);  
        combineFileName = strcat(curFld,'ManOverallTrends',svTxt,'.tif');
        % saveas(curFig,combineFileName);
        close(curFig)
    end
end

    
    
     %% Plot all metric trends - manual selected
%     subset = subset1;
%     subset2 = nonmon;
if(figGen==1)

    for t = 1:10
        if(t==1)
            % plot max DFF for all traces analyzed
            metricTrend = allTrendsMax;
            ytext = '\DeltaF/F';
            titleTxt = 'Maximum \DeltaF/F over Training';
            svTxt = 'MaxDFF';
            edges = 0:0.5:5;
        elseif(t==2)
            % plot onset time for all traces analyzed
            metricTrend = allTrendsOnset;
            ytext = 'ms';
            titleTxt = 'Spike Onset Time over Training';
            svTxt = 'Onset';
            edges = 0:50:1000;
        elseif(t==3)
            % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsDur;
            ytext = 'ms';
            titleTxt = 'Spike Duration over Training';
            svTxt = 'Duration';
            edges = 0:200:3000;
        elseif(t==4)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsN20;
            ytext = 'ms';
            titleTxt = 'Normalized Onset over Training';
            svTxt = 'NormalizedOnset';
            edges = 0:10:500;
        elseif(t==5)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsSpikes;
            ytext = 'Spikes';
            titleTxt = 'Deconvoluted Spiking over Training';
            svTxt = 'SpikeSum';
        elseif(t==6)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsSpikesItime;
            ytext = 'ms';
            titleTxt = 'Deconvoluted Spiking Onset over Training';
            svTxt = 'SpikeOnset';
        elseif(t==7)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsSpikesMtime;
            ytext = 'ms';
            titleTxt = 'Deconvoluted Spiking Max Time over Training';
            svTxt = 'SpikeMaxTime';
        elseif(t==8)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsSpikesDuration;
            ytext = 'ms';
            titleTxt = 'Deconvoluted Spiking Duration over Training';
            svTxt = 'SpikeDeconDur';
        elseif(t==9)
             % plot stimulation duration for all traces analyzed
            metricTrend = allTrendsSpikesPeaks;
            ytext = '# of Peaks';
            titleTxt = 'Deconvoluted number of midstim peaks over Training';
            svTxt = 'SpikePeaks'; 
        else
            % plot volume under stimulation curve for all traces analyzed
            metricTrend = allTrendsVol;
            ytext = '\DeltaF/F x ms';
            titleTxt = 'Spike Volume over Training';
            svTxt = 'Volume';
            edges = 0:10:200;
        end
        metricTrendInv = metricTrend(subset2,:);
        metricTrend = metricTrend(subset,:);

        curFig = figure();
        hold on
        for s = 1:size(metricTrendInv,2)
            h = ttest2(metricTrendInv(:,s),metricTrend(:,s));
            if(h==1)
                scatter(daysTrained(s),1,'*k')
            end
        end
        
        
        SEM = std(metricTrendInv,[],1,'omitnan')/sqrt(length(metricTrendInv));               % Standard Error
        ts = tinv(0.975  ,length(metricTrendInv)-1);      % T-Score
        CI_inv =  ts*SEM;                      % Confidence Intervals
        
        SEM = std(metricTrend,'omitnan')/sqrt(length(metricTrend));               % Standard Error
        ts = tinv( 0.975,length(metricTrend)-1);      % T-Score
        CI = ts*SEM;                      % Confidence Intervals
        
        errorbar(daysTrained,mean(metricTrendInv,1,'omitnan'),CI_inv,'blue','linewidth',2)
        errorbar(daysTrained+0.3,mean(metricTrend,1,'omitnan'),CI,'red','linewidth',2)
        xlim([min(daysTrained)-1, max(daysTrained)+1])
        
        
        xlabel('Days of Training')
        ylabel(ytext)
        title(titleTxt)
        combineFileName = strcat(curFld,'Box_ManOverallTrends',svTxt,'.fig');
        % saveas(curFig,combineFileName);  
        combineFileName = strcat(curFld,'Box_ManOverallTrends',svTxt,'.tif');
        % saveas(curFig,combineFileName);
        close(curFig)
    end
end

%% Plot metrics for non-active ROIs from session to session.  Perform a 
%  pusdo population analysis for comparision to the trends were observing
 % in our subpopulation
 
 if(figGen==1)
    
     for t = 1:10
            if(t==1)
                % plot max DFF for all traces analyzed
                metricTrend = endTrendsMax;
                metricTrendO = allTrendsMax;
                ytext = '\DeltaF/F';
                titleTxt = 'Maximum \DeltaF/F over Training';
                svTxt = 'MaxDFF';
                edges = 0:0.5:5;
            elseif(t==2)
                % plot onset time for all traces analyzed
                metricTrend = endTrendsOnset;
                metricTrendO = allTrendsOnset;
                ytext = 'ms';
                titleTxt = 'Spike Onset Time over Training';
                svTxt = 'Onset';
                edges = 0:50:1000;
            elseif(t==3)
                % plot stimulation duration for all traces analyzed
                metricTrend = endTrendsDur;
                metricTrendO = allTrendsDur;
                ytext = 'ms';
                titleTxt = 'Spike Duration over Training';
                svTxt = 'Duration';
                edges = 0:200:3000;
    
            elseif(t==5)
                 % plot stimulation duration for all traces analyzed
                metricTrend = endTrendsSpikes;
                metricTrendO = allTrendsSpikes;
                ytext = 'Spikes';
                titleTxt = 'Deconvoluted Spiking over Training';
                svTxt = 'SpikeSum';
            elseif(t==6)
                 % plot stimulation duration for all traces analyzed
                metricTrend = endTrendsSpikesItime;
                metricTrendO = allTrendsSpikesItime;
                ytext = 'ms';
                titleTxt = 'Deconvoluted Spiking Onset over Training';
                svTxt = 'SpikeOnset';
            elseif(t==7)
                 % plot stimulation duration for all traces analyzed
                metricTrend = endTrendsSpikesMtime;
                metricTrendO = allTrendsSpikesMtime;
                ytext = 'ms';
                titleTxt = 'Deconvoluted Spiking Max Time over Training';
                svTxt = 'SpikeMaxTime';
            elseif(t==8)
                 % plot stimulation duration for all traces analyzed
                metricTrend = endTrendsSpikesDuration;
                metricTrendO = allTrendsSpikesDuration;
                ytext = 'ms';
                titleTxt = 'Deconvoluted Spiking Duration over Training';
                svTxt = 'SpikeDeconDur';
            elseif(t==9)
                 % plot stimulation duration for all traces analyzed
                metricTrend = endTrendsSpikesPeaks;
                metricTrendO = allTrendsSpikesPeaks;
                ytext = '# of Peaks';
                titleTxt = 'Deconvoluted number of midstim peaks over Training';
                svTxt = 'SpikePeaks'; 
            else
                % plot volume under stimulation curve for all traces analyzed
                metricTrend = endTrendsVol;
                metricTrendO = allTrendsVol;
                ytext = '\DeltaF/F x ms';
                titleTxt = 'Spike Volume over Training';
                svTxt = 'Volume';
                edges = 0:10:200;
            end
                    
            curFig = figure();
            hold on
    
            endAveAll = mean(metricTrend,'all','omitnan');
            endValAll = reshape(metricTrend,1,[]);
            endValAll(isnan(endValAll))=[];
    %         endStdAll = std(metricTrend,[],'all','omitnan');
            SEM = std(endValAll,'omitnan')/sqrt(length(endValAll));               % Standard Error
            ts = tinv( 0.975,length(endValAll)-1);      % T-Score
            CIE = ts*SEM;                      % Confidence Intervals
            
            
            
            SEM = std(metricTrendO,'omitnan')/sqrt(length(metricTrendO));               % Standard Error
            ts = tinv( 0.975,length(metricTrendO)-1);      % T-Score
            CI = ts*SEM;                      % Confidence Intervals
            
            errorbar(daysTrained,mean(metricTrendO,1,'omitnan'),CI,'red','linewidth',2);
            xlim([min(daysTrained)-1, max(daysTrained)+1])
            
            errorbar(mean(daysTrained(numel(daysTrained)-1:end)),endAveAll,CIE,'Blue','linewidth',2);
    
            
            
            xlabel('Days of Training')
            ylabel(ytext)
            title(titleTxt)
            
            combineFileName = strcat(curFld,'EndOverallTrends',svTxt,'.fig');
            % saveas(curFig,combineFileName);  
            combineFileName = strcat(curFld,'EndOverallTrends',svTxt,'.tif');
            % saveas(curFig,combineFileName);
            close(curFig)
     end
 end
endDFF{animal} = endTrendsMax;
endOnset{animal} = endTrendsOnset;
endDur{animal} = endTrendsDur;
endSpikeNum{animal} = endTrendsSpikes;
endSpikeOnset{animal} = endTrendsSpikesItime;
endVol{animal} = endTrendsVol;
allEndROIs{animal} = endROIs2Compare;
endNeuroDat{animal} = endTrendsRoiCh;

 
%   %% Plot all metric trends - noncontinous selected
% %     subset = subset1;
% %     subset2 = nonmon;
%     for t = 1:8
%         if(t==1)
%             % plot max DFF for all traces analyzed
%             metricTrend = nonConTrendsMax;
%             ytext = '\DeltaF/F';
%             titleTxt = 'Maximum \DeltaF/F over Training';
%             svTxt = 'MaxDFF';
%             edges = 0:0.5:5;
%         elseif(t==2)
%             % plot onset time for all traces analyzed
%             metricTrend = nonConTrendsOnset;
%             ytext = 'ms';
%             titleTxt = 'Spike Onset Time over Training';
%             svTxt = 'Onset';
%             edges = 0:50:1000;
%         elseif(t==3)
%             % plot stimulation duration for all traces analyzed
%             metricTrend = nonConTrendsDur;
%             ytext = 'ms';
%             titleTxt = 'Spike Duration over Training';
%             svTxt = 'Duration';
%             edges = 0:200:3000;
%         elseif(t==4)
%              % plot stimulation duration for all traces analyzed
%             metricTrend = nonConTrendsN20;
%             ytext = 'ms';
%             titleTxt = 'Normalized Onset over Training';
%             svTxt = 'NormalizedOnset';
%             edges = 0:10:500;
%         elseif(t==5)
%              % plot stimulation duration for all traces analyzed
%             metricTrend = nonConTrendsSpikes;
%             ytext = 'Spikes';
%             titleTxt = 'Deconvoluted Spiking over Training';
%             svTxt = 'SpikeSum';
%         elseif(t==6)
%              % plot stimulation duration for all traces analyzed
%             metricTrend = nonConTrendsSpikesItime;
%             ytext = 'ms';
%             titleTxt = 'Deconvoluted Spiking Onset over Training';
%             svTxt = 'SpikeOnset';
%         elseif(t==7)
%              % plot stimulation duration for all traces analyzed
%             metricTrend = nonConTrendsSpikesMtime;
%             ytext = 'ms';
%             titleTxt = 'Deconvoluted Spiking Max Time over Training';
%             svTxt = 'SpikeMaxTime';
%         elseif(t==8)
%              % plot stimulation duration for all traces analyzed
%             metricTrend = nonConTrendsSpikesDuration;
%             ytext = 'ms';
%             titleTxt = 'Deconvoluted Spiking Duration over Training';
%             svTxt = 'SpikeDeconDur';
%         elseif(t==9)
%              % plot stimulation duration for all traces analyzed
%             metricTrend = nonConTrendsSpikesPeaks;
%             ytext = '# of Peaks';
%             titleTxt = 'Deconvoluted number of midstim peaks over Training';
%             svTxt = 'SpikePeaks'; 
%         else
%             % plot volume under stimulation curve for all traces analyzed
%             metricTrend = nonConTrendsVol;
%             ytext = '\DeltaF/F x ms';
%             titleTxt = 'Spike Volume over Training';
%             svTxt = 'Volume';
%             edges = 0:10:200;
%         end
% 
%         curFig = figure();
%         hold on
%         for i = 1:size(metricTrend,1)
%             plot(daysTrained,metricTrend(i,:),'color',[0 0 0 0.15],'linewidth',1)
%             scatter1 = scatter(daysTrained,metricTrend(i,:),'k');
%             scatter1.MarkerFaceAlpha = .15;
%             scatter1.MarkerEdgeAlpha = .15;
%         end
%         
%         
%         
% %         px = [min(daysTrained) max(daysTrained)];
% %         py = polyval(overallTrends(t,:), px);
% %         plot(px,py,'--','color','red','linewidth',2)
% %         
% %         aveAll = mean(metricTrendInv,1,'omitnan');
% %         plot(daysTrained,aveAll,'--','color','red','linewidth',2)
%         
%         aveAll = mean(metricTrend,1,'omitnan');
%         plot(daysTrained,aveAll,'color',[0 0 0],'linewidth',3)
%         xlim([min(daysTrained)-1, max(daysTrained)+1])
% 
%         
%         
%         xlabel('Days of Training')
%         ylabel(ytext)
%         title(titleTxt)
%         
%         combineFileName = strcat(curFld,'NonConOverallTrends',svTxt,'.fig');
%         saveas(curFig,combineFileName);  
%         combineFileName = strcat(curFld,'NonConnOverallTrends',svTxt,'.tif');
%         saveas(curFig,combineFileName);
%         close(curFig)
%     end



    %% Plot active region means and for subset 1
% %     subset1 = [8,9,12,20,22,35,39,40,44,53,71,80,81]; % increasing
% %         subset2 = [1,7,19,23,27,29,33,36,43,53,58,61,64,66]; % increasing
    if(figGen==1)

    aX = [];
    aY = [];
    curFig = figure();
    hold on
    roiSubSets = roiSets(subset2);
    overallMeans = NaN(1,numel(files2Load));
    for r = 1:numel(roiSubSets)
        curROIInds = roiSubSets{r};
        targMeans = NaN(numel(files2Load),1);
        targSTDs = NaN(numel(files2Load),1);
        for f = 1:numel(files2Load)
            curRegSTDs = regionPreStimSTDs{f};
            curRegMeans = regionPreStimMeans{f};

            curSessROIind = curROIInds(f);
            if(curSessROIind>0)
                targMeans(f) = curRegMeans(curSessROIind);
                targSTDs(f) = curRegSTDs(curSessROIind);
            end
        end
        aX = [aX, 1:numel(targMeans)];
        aY = [aY, targMeans];
        plot(1:numel(targMeans),targMeans,'color',[0 0 0 0.2],'linewidth',2)
        overallMeans(r,:)=targMeans;
    end
%     aX(isnan(aY))=[];
%     aY(isnan(aY))=[];
%     p = polyfit(aX, aY, 1);
%     px = [min(aX) max(aX)];
%     py = polyval(p, px);
    px = 1:numel(files2Load);
    py = mean(overallMeans,1,'omitnan');
    plot(px,py,'color',[0 0 0],'linewidth',3)
    title('Baseline Mean Trends Subset 2');
    xlabel('session')
    xticks(px)
    ylabel('Mean Baseline \DeltaF/F')
    combineFileName = strcat(curFld,'BaselineROIs-Subset2.tif');
    % saveas(curFig,combineFileName);
    close(curFig)
    end
    
    
    
    
    
    
    
%     %% Plot intensity normalized traces over training
%     subset1 = [8,9,12,20,22,35,39,40,44,53,71,80,81]; % increasing
%     roiSubSets = roiSets(subset1);
%     for r = 1:numel(roiSubSets)
%         curROIInds = roiSubSets{r};
%         curROI = allTrendsRoiCh(subset1(r),1);
%         curChan = allTrendsRoiCh(subset1(r),2);
%         curCurr = allTrendsRoiCh(subset1(r),3);
%         curChanInd = find(chans==curChan);
%         curCurrInd = find(allCurrs==curCurr);
%         
%         curFig = figure();
%         title(strcat('Normalized---ROI:',num2str(curROI),{' '},'Chan:',num2str(curChan),{' '},'Curr:',num2str(curCurr),'\muA'))
%         hold on
%         legSess = [];
%         for s = 1:numel(files2Load)
%             curSessROIind = curROIInds(s);
%             if(curSessROIind>0)
%                 legSess = [legSess, s];
%                 pre = preTraces{subset1(r),s,curChanInd,curCurrInd};
%                 midPost = midTraces{subset1(r),s,curChanInd,curCurrInd};
%                 dff = [pre,midPost];
%                 x = 33*((1:numel(dff))-numel(pre));
%                 dff = dff./max(dff);
%                 plot(x,dff,'linewidth',2)   
%             end
%         end
%         xlim([-1000 4000])
%         xline(0,'--');
%         xlabel('ms')
%         ylabel('\DeltaF/F')
%         legend(num2str(legSess'))
%                        
%         % Save figure
%         txt = strcat('Normalized---ROI',num2str(curROI),'Chan',num2str(curChan),'Curr',num2str(curCurr));
%         combineFileName = strcat(curFld,'Subset1----',txt,'.tif');
%         saveas(curFig,combineFileName);
% 
%         % zoom in on rising edge
%         xlim([-100 500])
%         combineFileName = strcat(curFld,'Subset1----RisingEdge',txt,'.tif');
%         saveas(curFig,combineFileName);
%         combineFileName = strcat(curFld,'Subset1----RisingEdge',txt,'.fig');
%         saveas(curFig,combineFileName);
%         close(curFig)    
%     end
% 
%     
%     
%     
    
    
    
    
    
% %     
% %     
% %     
% %     %% Plot active region means and for subset 2
% % %     subset2 = [1,7,19,23,27,29,33,36,43,53,58,61,64,66]; % increasing
% %     aX = [];
% %     aY = [];
% %     curFig = figure();
% %     hold on
% %     roiSubSets = roiSets(subset2);
% %     for r = 1:numel(roiSubSets)
% %         curROIInds = roiSubSets{r};
% %         targMeans = NaN(numel(files2Load),1);
% %         targSTDs = NaN(numel(files2Load),1);
% %         for f = 1:numel(files2Load) 
% %             curRegSTDs = regionSTDs{f};
% %             curRegMeans = regionPreStimMeans{f};
% % 
% %             curSessROIind = curROIInds(f);
% %             if(curSessROIind>0)
% %                 targMeans(f) = curRegMeans(curSessROIind);
% %                 targSTDs(f) = curRegSTDs(curSessROIind);
% %             end
% %         end
% %         aX = [aX, 1:numel(targMeans)];
% %         aY = [aY, targMeans];
% %         plot(1:numel(targMeans),targMeans,'-o','color',[0 0 0 0.2],'linewidth',2)
% %     end
% %     aX(isnan(aY))=[];
% %     aY(isnan(aY))=[];
% %     p = polyfit(aX, aY, 1);
% %     px = [min(aX) max(aX)];
% %     py = polyval(p, px);
% %     plot(px,py,'color',[0 0 0],'linewidth',3)
% %     title('Baseline Mean Trends Subset 2');
% %     combineFileName = strcat(curFld,'BaselineROIs-',num2str(r),'.tif');
% %     saveas(curFig,combineFileName);
% %     close(curFig)
% %       
    


     %% Manual plotting and labeling of datasets for roi selection
    %         
    % for f = 1:numel(files2Load)
    %     load(files2Load{f});
    %     
    %     curFig = figure('WindowState','maximized');
    %     imshow(meanPre./max(meanPre,[],'all'))
    %     hold on
    %     for r = 1:numel(regions)
    %         zProjectedBW = false(512,512);
    %         curPixels = regions(r).PixelList;
    %         for m = 1:size(curPixels,1)
    %             zProjectedBW(curPixels(m,2),curPixels(m,1))=1;
    %         end
    %         curCent = mean(curPixels,1);
    %         boundaries = bwboundaries(zProjectedBW);
    %         for k=1:numel(boundaries)
    %             b = boundaries{k};
    %             plot(b(:,2),b(:,1),'r','LineWidth',1);
    %         end
    %         text(curCent(1),curCent(2),num2str(r),'Color','y')
    %     end
    %     combineFileName = strcat('ICMS98_Sess',num2str(f),'.png');
    %     saveas(curFig,combineFileName);
    % end
end



%% merge all data based on weeks and just append all the data together
mergedDFF = NaN(1,5);
mergedTrendsOnset = NaN(1,5);
mergedTrendsDur =  NaN(1,5);
mergedTrendsSpikes =  NaN(1,5);
mergedTrendsSpikesItime =  NaN(1,5);
mergedTrendsVol = NaN(1,5);
mergedROI2Comp = NaN(1,5);
mergedChanData = NaN(1,3);
subsetMerged = [];
weekTimepoints = 0:4;
sAll=0;
sAll2=0;
for aInd = 1:3
    curWeeks = ceil(allDays{aInd}./7);
    subsetMerged = [subsetMerged; allSubset{aInd}];

    for w = 0:4
        tw = find(curWeeks==w);
        if(~isempty(tw))
            temp = allDFF{aInd};
            mergedDFF(sAll+1: sAll+size(allSubset{aInd},1), w+1) =  mean(temp(:,tw),2,'omitnan');
            temp = allOnset{aInd};
            mergedTrendsOnset(sAll+1: sAll+size(allSubset{aInd},1), w+1) = mean(temp(:,tw),2,'omitnan');
            temp =  allDur{aInd};
            mergedTrendsDur(sAll+1: sAll+size(allSubset{aInd},1), w+1) = mean(temp(:,tw),2,'omitnan');
            temp = allSpikeNum{aInd};
            mergedTrendsSpikes(sAll+1: sAll+size(allSubset{aInd},1), w+1) =  mean(temp(:,tw),2,'omitnan');
            temp = allSpikeOnset{aInd};
            mergedTrendsSpikesItime(sAll+1: sAll+size(allSubset{aInd},1), w+1) = mean(temp(:,tw),2,'omitnan');
            temp = allVol{aInd};
            mergedTrendsVol(sAll+1: sAll+size(allSubset{aInd},1), w+1) =  mean(temp(:,tw),2,'omitnan');
            temp = allROI2Comp{aInd};
            mergedROI2Comp(sAll+1: sAll+size(allROI2Comp{aInd},1), w+1) =  mean(temp(:,tw),2,'omitnan');


            % update index ref to unique 
            temp = allChanData{aInd};
            temp(:,1) = temp(:,1)+10000*aInd;
            mergedChanData(sAll+1: sAll+size(allSubset{aInd},1), :) = temp;
        end
    end
    sAll = sAll+size(allSubset{aInd},1);
    sAll2 = sAll2+size(allROI2Comp{aInd},1);
end

mergedTrendsVol(mergedTrendsVol==0)=NaN;
mergedTrendsSpikesItime(mergedTrendsSpikesItime==0)=NaN;
mergedTrendsSpikes(mergedTrendsSpikes==0)=NaN;
mergedTrendsDur(mergedTrendsDur==0)=NaN;
mergedTrendsOnset(mergedTrendsOnset==0)=NaN;
mergedDFF(mergedDFF==0)=NaN;
mergedROI2Comp(mergedROI2Comp==0)=NaN;

% exclude longitudinal tracks with will currents outside of 4 - 6 uA
invalidCurr = mergedChanData(:,3)<4 | mergedChanData(:,3)>6;
mergedTrendsVol(invalidCurr,:)=NaN;
mergedTrendsSpikesItime(invalidCurr,:)=NaN;
mergedTrendsSpikes(invalidCurr,:)=NaN;
mergedTrendsDur(invalidCurr,:)=NaN;
mergedTrendsOnset(invalidCurr,:)=NaN;
mergedDFF(invalidCurr,:)=NaN;

validLongNeurons = mergedChanData(:,1);
validLongNeurons(invalidCurr==1)=[];



%% Plot longitudinal activation raster of all ROIs

% Count of all weeks that had activation induced by stimulation -
% longitudinal responses
binROIRaster = mergedROI2Comp>0;
roiRastSum = sum(binROIRaster,2);

binROIRaster = binROIRaster(roiRastSum>0,:);
roiRastSum = roiRastSum(roiRastSum>0);

startInds = zeros(size(roiRastSum,1),1);
for i = 1:size(binROIRaster,1)
    startInds(i) = find(binROIRaster(i,:)>0,1,'first');
end
[~,inds] = sort(startInds);
binROIRaster = binROIRaster(inds,:);
roiRastSum = roiRastSum(inds);
heatRaster =zeros(size(binROIRaster));
for i = 1:size(binROIRaster,1)
    heatRaster(i,binROIRaster(i,:)==1) = roiRastSum(i);
end
heatRaster(4500,5)=0; % Pad raster to clean edge

figure()
h = heatmap(double(heatRaster));
h.GridVisible = 'off';




% % % for each neuron determine which weeks it was active <-- new B plot
% % binROIRaster = ~isnan(mergedDFF);
% % roiRastSum = sum(binROIRaster,2);
% % neuronIDs = unique(mergedChanData(:,1));
% % neuronROIRaster = zeros(numel(neuronIDs),5);
% % for n = 1:numel(neuronIDs)
% %     for i = 1:size(binROIRaster,1)
% %         if(mergedChanData(i,1)==neuronIDs(n))
% %             % Add activation to row
% %             neuronROIRaster(n,binROIRaster(i,:)==1)=1;
% %         end
% %     end
% % end
% % uniqueRoiRastSum = sum(neuronROIRaster,2);
% % 
% % 
% % singleNeurons = neuronIDs(uniqueRoiRastSum==1);
% % 
% % neuronROIRaster = neuronROIRaster(uniqueRoiRastSum>0,:);
% % uniqueRoiRastSum = uniqueRoiRastSum(uniqueRoiRastSum>0);
% % 
% % startInds = zeros(size(uniqueRoiRastSum,1),1);
% % for i = 1:size(neuronROIRaster,1)
% %     startInds(i) = find(neuronROIRaster(i,:)>0,1,'first');
% % end
% % [~,inds] = sort(startInds);
% % neuronROIRaster = neuronROIRaster(inds,:);
% % uniqueRoiRastSum = uniqueRoiRastSum(inds);
% % heatRaster =zeros(size(neuronROIRaster));
% % for i = 1:size(neuronROIRaster,1)
% %     heatRaster(i,neuronROIRaster(i,:)==1) = uniqueRoiRastSum(i);
% % end
% % figure()
% % h = heatmap(double(heatRaster));
% % h.GridVisible = 'off';

% % pair data with which neuron (unique) was for each longitudinal response
% neuronActiveWeeks = [mergedChanData(:,1),roiRastSum];






%% Merge end population dataset
mergedEndDFF = [];
mergedEndTrendsOnset = [];
mergedEndTrendsDur =  [];
mergedEndTrendsSpikes =  [];
mergedEndTrendsSpikesItime =  [];
mergedEndTrendsVol = [];
mergedEndChanData = [];
for aInd = 1:3
    % mergedEndDFF = [mergedEndDFF,reshape(endDFF{aInd},1,[])];
    % mergedEndTrendsOnset = [mergedEndTrendsOnset,reshape(endSpikeOnset{aInd},1,[])];
    % mergedEndTrendsDur =  [mergedEndTrendsDur,reshape(endDur{aInd},1,[])];
    % mergedEndTrendsSpikes =  [mergedEndTrendsSpikes,reshape(endSpikeNum{aInd},1,[])];
    % mergedEndTrendsSpikesItime =  [mergedEndTrendsSpikesItime,reshape(endSpikeOnset{aInd},1,[])];
    % mergedEndTrendsVol = [mergedEndTrendsVol,reshape(endVol{aInd},1,[])]; 
    % mergedEndChanData = [mergedEndChanData;endNeuroDat{aInd}];

    mergedEndDFF = [mergedEndDFF; mean(endDFF{aInd},2,'omitnan')];
    mergedEndTrendsOnset = [mergedEndTrendsOnset; mean(endSpikeOnset{aInd},2,'omitnan')];
    mergedEndTrendsDur =  [mergedEndTrendsDur;mean(endDur{aInd},2,'omitnan')];
    mergedEndTrendsSpikes =  [mergedEndTrendsSpikes; mean(endSpikeNum{aInd},2,'omitnan')];
    mergedEndTrendsSpikesItime =  [mergedEndTrendsSpikesItime; mean(endSpikeOnset{aInd},2,'omitnan')];
    mergedEndTrendsVol = [mergedEndTrendsVol; mean(endVol{aInd},2,'omitnan')]; 

end
mergedEndDFF(mergedEndDFF==0)=NaN;
mergedEndTrendsOnset(mergedEndTrendsOnset==0)=NaN;
mergedEndTrendsDur(mergedEndTrendsDur==0)=NaN;
mergedEndTrendsSpikes(mergedEndTrendsSpikes==0)=NaN;
mergedEndTrendsSpikesItime(mergedEndTrendsSpikesItime==0)=NaN;
mergedEndTrendsVol(mergedEndTrendsVol==0)=NaN;

% check through all rois and see which ones are present in longitudinal
% dataset
for aInd = 1:3
    temp = endNeuroDat{aInd};
    sharedRoiInd = zeros(size(temp,1),1);
    endROIs = allEndROIs{aInd};
    longROIs = allROI2Comp{aInd};
    allTemp = allChanData{aInd};
    newROICnt = max(allTemp(:,1));
    for s = 1:size(temp,1)
        searchROI = endROIs(temp(s),:);
        searchInd = find(searchROI);

        % find index
        longInd = find(longROIs(:,searchInd)==searchROI(searchInd));
        
        if(~isempty(longInd))
            sharedRoiInd(s) = longInd+10000*aInd;
        else
            newROICnt = newROICnt+1;
            sharedRoiInd(s) = newROICnt + 10000*aInd;
        end
    end

    temp(:,1) = sharedRoiInd;
    mergedEndChanData = [mergedEndChanData; temp];
end


% Remove entries outside of current range
invalidCurrEnd = mergedEndChanData(:,3)<4 | mergedEndChanData(:,3)>6;
mergedEndDFF(invalidCurrEnd,:)=NaN;
mergedEndTrendsOnset(invalidCurrEnd,:)=NaN;
mergedEndTrendsDur(invalidCurrEnd,:)=NaN;
mergedEndTrendsSpikes(invalidCurrEnd,:)=NaN;
mergedEndTrendsSpikesItime(invalidCurrEnd,:)=NaN;
mergedEndTrendsVol(invalidCurrEnd,:)=NaN;


% Remove all entries where the ROI was present in longitudinal analysis
longROIs = unique(mergedChanData(:,1));
invalidNeuronEnd = ismember(mergedEndChanData(:,1),longROIs);
mergedEndDFF(invalidNeuronEnd,:)=NaN;
mergedEndTrendsOnset(invalidNeuronEnd,:)=NaN;
mergedEndTrendsDur(invalidNeuronEnd,:)=NaN;
mergedEndTrendsSpikes(invalidNeuronEnd,:)=NaN;
mergedEndTrendsSpikesItime(invalidNeuronEnd,:)=NaN;
mergedEndTrendsVol(invalidNeuronEnd,:)=NaN;


% mergedEndChanData2 = [mergedEndChanData, isnan(mergedEndDFF)];
validEndNeurons = mergedEndChanData(:,1);
validEndNeurons(isnan(mergedEndDFF))=[];
validEndNeurons = unique(validEndNeurons);




 %% Plot all metric trends
    for t = 1:3
        if(t==1)
            % plot max DFF for all traces analyzed
            metricTrend = mergedDFF;
            ytext = '\DeltaF/F';
            titleTxt = 'Maximum \DeltaF/F over Training';
            svTxt = 'MaxDFF';
            edges = 0:0.5:5;
        elseif(t==2)
            % plot onset time for all traces analyzed
            metricTrend = mergedTrendsOnset;
            ytext = 'ms';
            titleTxt = 'Spike Onset Time over Training';
            svTxt = 'Onset';
            edges = 0:50:1000;
        elseif(t==3)
            % plot stimulation duration for all traces analyzed
            metricTrend = mergedTrendsDur;
            ytext = 'ms';
            titleTxt = 'Spike Duration over Training';
            svTxt = 'Duration';
            edges = 0:200:3000;
        end

        eGrp = zeros(size(metricTrend,1),1);
        lGrP = zeros(size(metricTrend,1),1);

        curFig = figure();
        hold on
        for i = 1:size(metricTrend,1)
            plot(weekTimepoints,metricTrend(i,:),'color',[0 0 0 0.15],'linewidth',1)
            scatter1 = scatter(weekTimepoints,metricTrend(i,:),'k');
            scatter1.MarkerFaceAlpha = .15;
            scatter1.MarkerEdgeAlpha = .15;
            eGrp(i) = mean(metricTrend(i,1:2),'omitnan');
            lGrP(i) = mean(metricTrend(i,3:5),'omitnan');
        end

        metSkip = isnan(eGrp);
        lGrP(isnan(eGrp))=[];
        eGrp(isnan(eGrp))=[];
        [p,~,stats] = ranksum(eGrp,lGrP);
        n1 = numel(eGrp);
        n2 = numel(lGrP);
        Z = stats.zval;
        p
        r = abs(Z) / sqrt(n1 + n2)

        curFig = figure();
        hold on
        for i = 1:size(metricTrend,1)
            if(metSkip(i)==0)
                x = weekTimepoints;
                y = metricTrend(i,:);
                x(isnan(y))=[];
                y(isnan(y))=[];
                scatter1 = scatter(x,y,'k');
                scatter1.MarkerFaceAlpha = .15;
                scatter1.MarkerEdgeAlpha = .15;
                P = polyfit(x,y,1);
                yfit = P(1)*x+P(2);
                plot(x,yfit,'color',[0 0 0 0.15],'linewidth',1);
                % plot(weekTimepoints,metricTrend(i,:),'color',[0 0 0 0.15],'linewidth',1)
            end
        end


        aveAll = mean(metricTrend,1,'omitnan');

        P = polyfit(weekTimepoints,aveAll,1);
        yfit = P(1)*weekTimepoints+P(2);
        plot(weekTimepoints,yfit,'color',[0 0 0],'linewidth',3)

        % plot(weekTimepoints,aveAll,'color',[0 0 0],'linewidth',3)
        xlim([min(weekTimepoints)-0.5, max(weekTimepoints)+0.5])
        xticks([0:4])

        xlabel('Weeks of Training')
        ylabel(ytext)
        title(titleTxt)        
        combineFileName = strcat(curFld,'MergedOverallTrends',svTxt,'.fig');
        % saveas(curFig,combineFileName);  
        combineFileName = strcat(curFld,'MergedOverallTrends',svTxt,'.tif');
        % saveas(curFig,combineFileName);
        close(curFig)

        dataTrend = [];
        gTrend = [];
        for g = 1:5
            temp = metricTrend(:,g);
            dataTrend = [dataTrend; temp];
            gTrend = [gTrend; ones(numel(temp),1)*g];
        end
        gTrend(isnan(dataTrend))=[];
        dataTrend(isnan(dataTrend))=[];


        % gTrend = double(gTrend>2);

        % [p,~,stats] = ranksum(dataTrend(gTrend<3),dataTrend(gTrend>2));
       
        % p
        % [p,tbl,stats] = kruskalwallis(dataTrend,gTrend);
        % c = multcompare(stats);
    end




%% Plot longituidnal comparision of early and late rois
figure()
histogram(reshape(mergedDFF(:,1:2),1,[]),40)
hold on
xline(mean(reshape(mergedDFF(:,1:2),1,[]),'omitnan'),'r')
xlim([0.2 3.7])
ylim([0 30])
xticks([0.5,1.5,2.5,3.5])
yticks([0,15,30])

figure()
hold on
histogram(reshape(mergedDFF(:,4:5),1,[]),40)
xline(mean(reshape(mergedDFF(:,4:5),1,[]),'omitnan'),'r')
threshSTD = mean(reshape(mergedDFF(:,4:5),1,[]),'omitnan') + std(reshape(mergedDFF(:,4:5),1,[]),[],'omitnan');
xline(threshSTD,'r--')
xlim([0.2 3.7])
ylim([0 30])
xticks([0.5,1.5,2.5,3.5])
yticks([0,15,30])


 %% Plot all metric trends - manual selected
%     subset = subset1;
%     subset2 = nonmon;
% % % % %     for t = 1:6
% % % % %         if(t==1)
% % % % %             % plot max DFF for all traces analyzed
% % % % %             metricTrend = mergedDFF;
% % % % %             ytext = '\DeltaF/F';
% % % % %             titleTxt = 'Maximum \DeltaF/F over Training';
% % % % %             svTxt = 'MaxDFF';
% % % % %             edges = 0:0.5:5;
% % % % %         elseif(t==2)
% % % % %             % plot onset time for all traces analyzed
% % % % %             metricTrend = mergedTrendsOnset;
% % % % %             ytext = 'ms';
% % % % %             titleTxt = 'Spike Onset Time over Training';
% % % % %             svTxt = 'Onset';
% % % % %             edges = 0:50:1000;
% % % % %         elseif(t==3)
% % % % %             % plot stimulation duration for all traces analyzed
% % % % %             metricTrend = mergedTrendsDur;
% % % % %             ytext = 'ms';
% % % % %             titleTxt = 'Spike Duration over Training';
% % % % %             svTxt = 'Duration';
% % % % %             edges = 0:200:3000;
% % % % %         elseif(t==4)
% % % % %             % plot stimulation duration for all traces analyzed
% % % % %             metricTrend = mergedTrendsSpikes;
% % % % %             ytext = 'Spikes';
% % % % %             titleTxt = 'Deconvoluted Spiking over Training';
% % % % %             svTxt = 'SpikeSum';
% % % % %         elseif(t==5)
% % % % %              % plot stimulation duration for all traces analyzed
% % % % %             metricTrend = mergedTrendsSpikesItime;
% % % % %             ytext = 'ms';
% % % % %             titleTxt = 'Deconvoluted Spiking Onset over Training';
% % % % %             svTxt = 'SpikeOnset';
% % % % %         else
% % % % %             % plot volume under stimulation curve for all traces analyzed
% % % % %             metricTrend = mergedTrendsVol;
% % % % %             ytext = '\DeltaF/F x ms';
% % % % %             titleTxt = 'Spike Volume over Training';
% % % % %             svTxt = 'Volume';
% % % % %             edges = 0:10:200;
% % % % %         end
% % % % %         metricTrendInv = metricTrend(subsetMerged==0,:);
% % % % %         metricTrend = metricTrend(subsetMerged==1,:);
% % % % % 
% % % % %         curFig = figure();
% % % % %         hold on
% % % % %         for i = 1:size(metricTrend,1)
% % % % %             plot(weekTimepoints,metricTrend(i,:),'color',[0 0 0 0.15],'linewidth',1)
% % % % %             scatter1 = scatter(weekTimepoints,metricTrend(i,:),'k');
% % % % %             scatter1.MarkerFaceAlpha = .15;
% % % % %             scatter1.MarkerEdgeAlpha = .15;
% % % % %         end
% % % % % 
% % % % %         aveAll = mean(metricTrendInv,1,'omitnan');
% % % % %         plot(weekTimepoints,aveAll,'--','color','red','linewidth',2)
% % % % % 
% % % % %         aveAll = mean(metricTrend,1,'omitnan');
% % % % %         plot(weekTimepoints,aveAll,'color',[0 0 0],'linewidth',3)
% % % % %         xlim([min(weekTimepoints)-0.5, max(weekTimepoints)+0.5])
% % % % %         xticks([0:4])
% % % % % 
% % % % %         xlabel('Weeks of Training')
% % % % %         ylabel(ytext)
% % % % %         title(titleTxt)        
% % % % %         combineFileName = strcat(curFld,'MergedOverallTrends',svTxt,'.fig');
% % % % %         saveas(curFig,combineFileName);  
% % % % %         combineFileName = strcat(curFld,'MergedOverallTrends',svTxt,'.tif');
% % % % %         saveas(curFig,combineFileName);
% % % % % %         close(curFig)
% % % % %     end
% % % % % 

    for t = 1:6
        % if(t==1)
        %     % plot max DFF for all traces analyzed
        %     metricTrend = mergedDFF;
        %     ytext = '\DeltaF/F';
        %     titleTxt = 'Maximum \DeltaF/F over Training';
        %     svTxt = 'MaxDFF';
        %     edges = 0:0.5:5;
        % elseif(t==2)
        %     % plot onset time for all traces analyzed
        %     metricTrend = mergedTrendsOnset;
        %     ytext = 'ms';
        %     titleTxt = 'Spike Onset Time over Training';
        %     svTxt = 'Onset';
        %     edges = 0:50:1000;
        % elseif(t==3)
        %     % plot stimulation duration for all traces analyzed
        %     metricTrend = mergedTrendsDur;
        %     ytext = 'ms';
        %     titleTxt = 'Spike Duration over Training';
        %     svTxt = 'Duration';
        %     edges = 0:200:3000;
        % elseif(t==4)
        %      % plot stimulation duration for all traces analyzed
        %     metricTrend = mergedTrendsSpikes;
        %     ytext = 'Spikes';
        %     titleTxt = 'Deconvoluted Spiking over Training';
        %     svTxt = 'SpikeSum';
        % elseif(t==5)
        %      % plot stimulation duration for all traces analyzed
        %     metricTrend = mergedTrendsSpikesItime;
        %     ytext = 'ms';
        %     titleTxt = 'Deconvoluted Spiking Onset over Training';
        %     svTxt = 'SpikeOnset';
        % else
        %     % plot volume under stimulation curve for all traces analyzed
        %     metricTrend = mergedTrendsVol;
        %     ytext = '\DeltaF/F x ms';
        %     titleTxt = 'Spike Volume over Training';
        %     svTxt = 'Volume';
        %     edges = 0:10:200;
        % end
         if(t==1)
            % plot max DFF for all traces analyzed
            metricTrend = mergedDFF;
            ytext = '\DeltaF/F';
            titleTxt = 'Maximum \DeltaF/F over Training';
            svTxt = 'MaxDFF';
            edges = 0:0.5:5;

        elseif(t==2)
            % Deconvolved spiking
            metricTrend = mergedTrendsSpikes;
            ytext = 'Spikes';
            titleTxt = 'Deconvoluted Spiking over Training';
            svTxt = 'SpikeSum';
        elseif(t==3)
             % plot stimulation duration for all traces analyzed
            metricTrend = mergedTrendsSpikesItime;
            ytext = 'ms';
            titleTxt = 'Deconvoluted Spiking Onset over Training';
            svTxt = 'SpikeOnset';
        elseif(t==4)
             % plot onset time for all traces analyzed
            metricTrend = mergedTrendsOnset;
            ytext = 'ms';
            titleTxt = 'Spike Onset Time over Training';
            svTxt = 'Onset';
        elseif(t==5)
            % plot volume under stimulation curve for all traces analyzed
            metricTrend = mergedTrendsVol;
            ytext = '\DeltaF/F x ms';
            titleTxt = 'Spike Volume over Training';
            svTxt = 'Volume';
        elseif(t==6)
            % plot stimulation duration for all traces analyzed
            metricTrend = mergedTrendsDur;
            ytext = 'ms';
            titleTxt = 'Spike Duration over Training';
            svTxt = 'Duration';
            edges = 0:200:3000;
        end
        metricTrendInv = metricTrend(subsetMerged==0,:);
        metricTrend = metricTrend(subsetMerged==1,:);
svTxt
        % data =[];
        % g = [];
        % for k = 1:size(metricTrend,2)
        %     t = metricTrend(:,k);
        %     data = [data;t];
        %     gA = ones(1,numel(t))*k;
        %     g = [g,gA];
        % end
        % g(isnan(data))=[];
        % data(isnan(data))=[];
        data = zeros(size(metricTrend,1),2);
        for k = 1:size(metricTrend,1)
            data(k,1) = mean(metricTrend(k,1:2),'omitnan');
            data(k,2) = mean(metricTrend(k,3:5),'omitnan');
        end
        % [~,~,stats] = anovan(data,g');
        % c1 = multcompare(stats)
        data(isnan(data(:,1)),:)=[];
        data(isnan(data(:,2)),:)=[];
        % p = ranksum(data(:,1),data(:,2));
        [p,~,stats] = ranksum(data(:,1),data(:,2));
        n1 = numel(data(:,1));
        n2 = numel(data(:,2));
        Z = stats.zval;
        p
        r = abs(Z) / sqrt(n1 + n2)

        % % % data =[];
        % % % g = [];
        % % % for k = 1:size(metricTrendInv,2)
        % % %     t = metricTrendInv(:,k);
        % % %     data = [data;t];
        % % %     gA = ones(1,numel(t))*k;
        % % %     g = [g,gA];
        % % % end
        % % % g(isnan(data))=[];
        % % % data(isnan(data))=[];
        % % % % [~,~,stats] = anovan(data,g');
        % % % % c1 = multcompare(stats)
        % % % ranksum(data(g<3),dlata(g>2))

        data = zeros(size(metricTrendInv,1),2);
        for k = 1:size(metricTrendInv,1)
            data(k,1) = mean(metricTrendInv(k,1:2),'omitnan');
            data(k,2) = mean(metricTrendInv(k,3:5),'omitnan');
        end
        % [~,~,stats] = anovan(data,g');
        % c1 = multcompare(stats)
        data(isnan(data(:,1)),:)=[];
        data(isnan(data(:,2)),:)=[];
        % ranksum(data(:,1),data(:,2))
        [p,~,stats] = ranksum(data(:,1),data(:,2));
        n1 = numel(data(:,1));
        n2 = numel(data(:,2));
        Z = stats.zval;
        p
        r = abs(Z) / sqrt(n1 + n2)

        svTxt
        curFig = figure();
        hold on
        % for s = 1:size(metricTrendInv,2)
        %     h = ttest2(metricTrendInv(:,s),metricTrend(:,s));
        %     if(h==1)
        %         scatter(weekTimepoints(s),1,'*k')
        %     end
        % end
        for s = 1:size(metricTrendInv,2)
            [p,~,stats] = ranksum(metricTrendInv(:,s),metricTrend(:,s))
            n1 = numel(data(:,1));
            n2 = numel(data(:,2));
            Z = stats.zval;
            p
            r = abs(Z) / sqrt(n1 + n2)
            if(p<0.001)
                scatter(weekTimepoints(s),-1,'*r')
            elseif(p<0.01)
                scatter(weekTimepoints(s),-1,'*b')
            elseif(p<0.05)
                scatter(weekTimepoints(s),-1,'*k')
            end
        end

        % SEM = std(metricTrendInv,[],1,'omitnan')/sqrt(length(metricTrendInv));               % Standard Error
        % ts = tinv(0.975  ,length(metricTrendInv)-1);      % T-Score
        % CI_inv =  ts*SEM;                      % Confidence Intervals
        % 
        % SEM = std(metricTrend,'omitnan')/sqrt(length(metricTrend));               % Standard Error
        % ts = tinv( 0.975,length(metricTrend)-1);      % T-Score
        % CI = ts*SEM;                      % Confidence Intervals
        % 
        % errorbar(weekTimepoints,mean(metricTrendInv,1,'omitnan'),CI_inv,'blue','linewidth',2)
        % errorbar(weekTimepoints,mean(metricTrend,1,'omitnan'),CI,'red','linewidth',2)

        errorbar(weekTimepoints,mean(metricTrendInv,1,'omitnan'),std(metricTrendInv,[],1,'omitnan'),'blue','linewidth',2)
        errorbar(weekTimepoints+0.2,mean(metricTrend,1,'omitnan'),std(metricTrend,[],1,'omitnan'),'red','linewidth',2)
        xlim([min(weekTimepoints)-0.5, max(weekTimepoints)+0.5])
        xticks([0:4])


        xlabel('Weeks of Training')
        ylabel(ytext)
        title(titleTxt)
        combineFileName = strcat(curFld,'MergedBox_ManOverallTrends',svTxt,'.fig');
        % saveas(curFig,combineFileName);  
        combineFileName = strcat(curFld,'MergedBox_ManOverallTrends',svTxt,'.tif');
        % saveas(curFig,combineFileName);
        close(curFig)
    end
    
    
   %% Plot select group from end neurons
     for t = 1:4
        if(t==1)
            % plot max DFF for all traces analyzed
            metricTrend2 = mergedEndDFF;
            metricTrend = mergedDFF;
            ytext = '\DeltaF/F';
            titleTxt = 'Maximum \DeltaF/F over Training';
            svTxt = 'MaxDFF';
            edges = 0:0.5:5;
        elseif(t==2)
        %     % plot onset time for all traces analyzed
        %     metricTrend2 = mergedEndTrendsOnset;
        %     metricTrend = mergedTrendsOnset;
        %     ytext = 'ms';
        %     titleTxt = 'Spike Onset Time over Training';
        %     svTxt = 'Onset';
        %     edges = 0:50:1000;
        % elseif(t==3)
        %     % plot stimulation duration for all traces analyzed
        %     metricTrend2 = mergedEndTrendsDur;
        %     metricTrend = mergedTrendsDur;
        %     ytext = 'ms';
        %     titleTxt = 'Spike Duration over Training';
        %     svTxt = 'Duration';
        %     edges = 0:200:3000;
        % elseif(t==4)
            metricTrend2 = mergedEndTrendsSpikes;
            metricTrend = mergedTrendsSpikes;
            ytext = 'Spikes';
            titleTxt = 'Deconvolved Spiking over Training';
            svTxt = 'SpikeSum';
        elseif(t==3)
            metricTrend2 = mergedEndTrendsSpikesItime;
            metricTrend = mergedTrendsSpikesItime;
            ytext = 'ms';
            titleTxt = 'Deconvolved Spiking Onset over Training';
            svTxt = 'SpikeOnset';
        else
            % plot volume under stimulation curve for all traces analyzed
            metricTrend2 = mergedEndTrendsVol;
            metricTrend = mergedTrendsVol;
            ytext = '\DeltaF/F x ms';
            titleTxt = 'Spike Volume over Training';
            svTxt = 'Volume';
            edges = 0:10:200;
        end
%         metricTrendInv = metricTrend(subsetMerged==0,:);
        metricTrend = metricTrend(subsetMerged==1,:);

        curFig = figure();
        hold on
        for s = 1:size(metricTrend,2)
            [p,~,stats] = ranksum(metricTrend2,metricTrend(:,s));
            p
            totN = [metricTrend2;metricTrend(:,s)];
%             h = ttest2(metricTrend2,metricTrend(:,s));
            if(p<0.001)
                scatter(weekTimepoints(s),1,'*r')
            elseif(p<0.01)
                scatter(weekTimepoints(s),1,'*b')
            elseif(p<0.05)
                scatter(weekTimepoints(s),1,'*k')
            end
        end


        % y = reshape(metricTrend,1,[])';
        % g = [];
        % for s = 1:size(metricTrend,2)
        %     g = [g; ones(size(metricTrend,1),1)*(s-1)];
        % end
        % 
        % [~,~,stats] = anovan(y,double(g>1));
        % c = multcompare(stats);

        
        
        % SEM = std(metricTrend2,[],'omitnan')/sqrt(length(metricTrend2));               % Standard Error
        % ts = tinv(0.975  ,length(metricTrend2)-1);      % T-Score
        % CI_inv =  ts*SEM;                      % Confidence Intervals
        % 
        % SEM = std(metricTrend,'omitnan')/sqrt(length(metricTrend));               % Standard Error
        % ts = tinv( 0.975,length(metricTrend)-1);      % T-Score
        % CI = ts*SEM;                      % Confidence Intervals
        % 
        % 
        errorbar(4.3,mean(metricTrend2,'omitnan'),std(metricTrend2,'omitnan'),'o',"MarkerEdgeColor",'green',"MarkerFaceColor",'green','linewidth',2)
        errorbar(weekTimepoints,mean(metricTrend,1,'omitnan'),std(metricTrend,[],1,'omitnan'),'red','linewidth',2)
        % errorbar(4.3,mean(metricTrend2,'omitnan'),CI_inv,'blue','linewidth',2)
        % errorbar(weekTimepoints,mean(metricTrend,1,'omitnan'),CI,'red','linewidth',2)
        xlim([min(weekTimepoints)-0.5, max(weekTimepoints)+0.5])
        xticks([0:4])
        
        
        xlabel('Weeks of Training')
        ylabel(ytext)
        title(titleTxt)
        combineFileName = strcat(curFld,'MergedBox_Endcompare_ManOverallTrends',svTxt,'.fig');
        % saveas(curFig,combineFileName);  
        combineFileName = strcat(curFld,'MergedBox_Endcompare_ManOverallTrends',svTxt,'.tif');
        % saveas(curFig,combineFileName);
        close(curFig)
    end
    
    
    




     %% Plot select group from end neurons
     for t = 1:6
        if(t==1)
            % plot max DFF for all traces analyzed
            metricTrend2 = mergedEndDFF;
            metricTrend = mergedDFF;
            ytext = '\DeltaF/F';
            titleTxt = 'Maximum \DeltaF/F over Training';
            svTxt = 'MaxDFF';
            edges = 0:0.5:5;
        elseif(t==2)
            % plot onset time for all traces analyzed
            metricTrend2 = mergedEndTrendsOnset;
            metricTrend = mergedTrendsOnset;
            ytext = 'ms';
            titleTxt = 'Spike Onset Time over Training';
            svTxt = 'Onset';
            edges = 0:50:1000;
        elseif(t==3)
            % plot stimulation duration for all traces analyzed
            metricTrend2 = mergedEndTrendsDur;
            metricTrend = mergedTrendsDur;
            ytext = 'ms';
            titleTxt = 'Spike Duration over Training';
            svTxt = 'Duration';
            edges = 0:200:3000;
        elseif(t==4)
            metricTrend2 = mergedEndTrendsSpikes;
            metricTrend = mergedTrendsSpikes;
            ytext = 'Spikes';
            titleTxt = 'Deconvolved Spiking over Training';
            svTxt = 'SpikeSum';
        elseif(t==5)
            metricTrend2 = mergedEndTrendsSpikesItime;
            metricTrend = mergedTrendsSpikesItime;
            ytext = 'ms';
            titleTxt = 'Deconvolved Spiking Onset over Training';
            svTxt = 'SpikeOnset';
        else
            % plot volume under stimulation curve for all traces analyzed
            metricTrend2 = mergedEndTrendsVol;
            metricTrend = mergedTrendsVol;
            ytext = '\DeltaF/F x ms';
            titleTxt = 'Spike Volume over Training';
            svTxt = 'Volume';
            edges = 0:10:200;
        end
        metricTrendAll = metricTrend;
        metricTrend = metricTrend(subsetMerged==1,:);
        earDat = mean(metricTrend(:,1),2,'omitnan');
        latDat = mean(metricTrend(:,3:5),2,'omitnan');
        latDat(isnan(latDat))=[];
        earDat(isnan(earDat))=[];
        metricTrend2(isnan(metricTrend2))=[];

        barDat = [earDat',latDat',metricTrend2'];
        barG = [ones(1,numel(earDat))*1, ones(1,numel(latDat))*2, ones(1,numel(metricTrend2))*3];
        data = [mean(earDat,'omitnan'),mean(latDat,'omitnan'),mean(metricTrend2,'omitnan')];
        stDevs = [std(earDat,'omitnan'),std(latDat,'omitnan'),std(metricTrend2,'omitnan')];
        errHigh = data ;
        errLow = data - stDevs;

        curFig = figure();
        hold on
        bar(1:3,data)
        er = errorbar(1:3,data,stDevs);    
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';  
        xlabel('Weeks of Training')
        ylabel(ytext)
        title(titleTxt)


        titleTxt
        [p,tbl,stats] = kruskalwallis(barDat,barG);
        c = multcompare(stats)
        % ranksum(metricTrend2,latDat)


        % validTest = zeros(1000,3);
        % for tr = 1:1000
        %     metricTrendS = randsample(metricTrend2,50);
        %     validTest(tr,1) = ranksum(earDat,latDat);
        %     validTest(tr,2) = ranksum(earDat,metricTrendS);
        %     validTest(tr,3) = ranksum(metricTrendS,latDat);
        % end
        % vres = mean(validTest,1);

    end
    
    
    

 %%

     curDists = [];
    curDistG = [];
    for aInd = 1:3
        temp = distanceNon{aInd}';
        curDists = [curDists;temp];
        curDistG = [curDistG;ones(size(temp,1),1)*2];
        
        temp = distanceSub{aInd}';
        curDists = [curDists;temp];
        curDistG = [curDistG;ones(size(temp,1),1)*1];
    end

 figure()
boxplot(curDists,curDistG)
title('Comparision of Subgroup Vs Main Population Activation Distance')
ylabel('Activation Distance (\mum)')

[~,~,stats] = anovan(curDists,curDistG)
c1 = multcompare(stats);