% Requires raw session data (not in repo).
% Internal paths must be configured manually for your local environment.

cd('S:\Roy\robin - Copy\planar_analysis');
addpath(genpath('S:\Roy\ryhattori-PatchWarp-1.3.3.0'));

% calculate trace means and stds
% TODO - report baseline trace means over each session, i.e. see if there is a drifting baseline in each session or for each roi
allROI2Comp = cell(3,1);
allDays = cell(3,1);


for animal =  1:3
    if(animal==1)
    icms = 'ICMS92';
    files2Load = {'S:\ICMS92\9-6-23\PlanarTrialsG8\SessionMetrics.mat';...
                'S:\ICMS92\9-8-23\PlanarTrialsG8\SessionMetrics.mat';...
                'S:\ICMS92\9-12-23\PlanarTrialsG8\SessionMetrics.mat';...
                'S:\ICMS92\9-14-23\PlanarTrialsG8\SessionMetrics.mat';...
                'S:\ICMS92\9-19-23\PlanarTrialsG8\SessionMetrics.mat';...
                'S:\ICMS92\9-21-23\PlanarTrialsG8\SessionMetrics.mat';...
                'S:\ICMS92\9-25-23\PlanarTrialsG8\SessionMetrics.mat'};
            
    % minSetCnt = 4;
    minSetCnt = 1;

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
        % minSetCnt = 3;
        minSetCnt = 1;

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
        % minSetCnt = 3;
        minSetCnt = 1;

        chans = [9,11,12]; % hardcoded for this animal - update for each animal
        daysTrained = [0,7,15,23,27];     
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
        end

        %save alignment
        save(alignFile,'ptForms');

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



    % find overlapping regions across sessions
    roiSets = cell(1,1);
    roiSetsPx = cell(1,1);
    rsC = 0;
    bob=0;
    for f = 1:numel(files2Load)
        curRegions = updatedRegions{f};
        curProc = processedRegions{f};
        for r1 = 1:numel(curRegions) % iterate through all regions and find subsequent regions that overlap
            compCnt =0;
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
                                compCnt = compCnt+1;
                            end
                        end
                    end

                    if(optInd>0) % if these was an overlap add to current set and make region as processed
                        curSet(f2) = optInd;
                        % if(compCnt>1)
                        %     bob = bob+1
                        % end
                        searchProc(optInd) = 1;
                        processedRegions{f2} = searchProc;
                    end
                end

                rsC = rsC+1;
                roiSets{rsC} = curSet;
                roiSetsPx{rsC} = curRegions{r1};
            end
        end
    end


    
    
    AllData = cell(numel(files2Load),1);
    for s = 1:numel(files2Load)
        AllData{s} = load(files2Load{s});
    end


    %% Generate traces for all ROIs selected
    ROIs2Compare = zeros(numel(roiSets),numel(files2Load));
    for r = 1:numel(roiSets) % convert keys to single matrix
        ROIs2Compare(r,:) = roiSets{r};
    end
    allROI2Comp{animal} = ROIs2Compare;
    allDays{animal} = daysTrained;
end




%% merge all data based on weeks and just append all the data together
mergedROI2Comp = zeros(1,5);
weekTimepoints = 0:4;
sAll=0;
for aInd = 1:3
    curWeeks = ceil(allDays{aInd}./7);

    for w = 0:4
        tw = find(curWeeks==w);
        if(~isempty(tw))
            temp = allROI2Comp{aInd};
            mergedROI2Comp(sAll+1: sAll+size(allROI2Comp{aInd},1), w+1) =  mean(temp(:,tw),2,'omitnan');

        end
    end
    sAll = sAll+size(allROI2Comp{aInd},1);
end
binROIRaster = mergedROI2Comp>0;
roiRastSum = sum(binROIRaster,2);
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
h = heatmap(double(heatRaster));
h.GridVisible = 'off';
% mergedROI2Comp(mergedROI2Comp==0)=NaN;

