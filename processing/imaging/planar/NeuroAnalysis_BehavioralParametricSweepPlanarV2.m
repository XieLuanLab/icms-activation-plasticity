function [] = NeuroAnalysis_BehavioralParametricSweepPlanarV2(sourceFolders)
% Program which processes the data aquired during behavioral paired imaging
% and stimulation experiments.
     
    warning('off','all')
    warning('off','MATLAB:MKDIR:DirectoryExists');
    elecPos = zeros(32,1); % Possible electrode position file, will be populated by external file (if it exists)
    
    % Load supporting files
    addpath(genpath('./rasterplot'));
    addpath(genpath('./natsortfiles'));
    addpath(genpath('./imshow3D'));
    addpath(genpath('./bfmatlab'));
    addpath(genpath('./strel3d'));
    addpath(genpath('./mij'));
    javaaddpath 'C:\Program Files\MATLAB\R2019a\java\mij.jar'
    javaaddpath 'C:\Program Files\MATLAB\R2019a\java\ij.jar'

    % Start ImageJ interface
    %    If you need to load the libraries look here: imagej.net/plugins/miji
    for curFile = 1:numel(sourceFolders) %Process each experiment folder/dataset

        selpath = sourceFolders{curFile}; % select folder containing raw files
        fprintf('\nProcessing Dataset: ');
        pathDispName = erase(selpath,'\\10.129.151.108\nei1\xieluanlabs\xl_stimulation\');
        pathDispName = erase(pathDispName,'\RAW');
        disp(pathDispName)

        idcs = strfind(selpath,'\');
        pullPath = selpath(1:idcs(end)-1);
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

        % Extract key files with timing of all stimulation and  imaging 
        keyListing = dir(keypathO);
        keyFile = keyListing(3).name; % There should only be one key per session, grab it
        
        


        %% Process raw images into individual stacks 
        procPath = strcat(pullPath,'\PROCCESSED-PLANAR');
        if(~exist(procPath, 'dir'))
            mkdir(procPath) % if directory doesn't exist make it
            fprintf('Performing initial file conversion  ');
            
            % Open raw data folder - only one for planar data
            curDir = strcat(selpath,'\',segListing(1).name);
            rawListing = dir(curDir);
            rawNames = cell(numel(rawListing),1);
            for k = 1:numel(rawListing)
                rawNames{k} = rawListing(k).name;
            end
            % Discard files not from Ch3 (this excludes 
            validInds = find(contains(rawNames,'Ch3')); % check to make sure additional imaging channels are not processed
            
            
            
            % Sort into cycles (using a lazy loop)
            cyclePos = zeros(numel(validInds),2);
            for k = 1:numel(validInds)
                curStr = rawNames{validInds(k)};
                newStr1 = split(curStr,'Cycle');
                newStr2 = split(newStr1{2},'_');
                newStr3 = split(newStr2{3},'.');
                cyclePos(k,1)= str2double(newStr2{1});
                cyclePos(k,2)= str2double(newStr3{1});
            end
            [imgNum,imgInd] = sortrows(cyclePos,1);
            curCycles = unique(imgNum(:,1));

            % %---------- MANUAL OPTION FOR ICMS 101 11/1/23 ---------%
            % % find groups of cycles
            % overrideCycles = cell(1);
            % oCnt = 0;
            % priorCycle = curCycles(1,1);
            % cycleSet = curCycles(1,1);
            % for cN = 2:size(curCycles,1)
            %     if(curCycles(cN,1) == priorCycle+1)
            %       % continue cycle
            %       cycleSet = [cycleSet,curCycles(cN,1)];
            %     else
            %       % new cycle
            %       oCnt = oCnt+1;
            %       overrideCycles{oCnt} = cycleSet;
            %       cycleSet = curCycles(cN,1);
            %     end
            %     priorCycle = curCycles(cN,1);
            % end
            % % save last cycle
            % oCnt = oCnt+1;
            % overrideCycles{oCnt} = cycleSet;
            % 
            % % update cycle data
            % newImgNum = imgNum;
            % for h = 1:oCnt
            %     curCycleSet = overrideCycles{h};
            %     for j = 1:numel(curCycleSet)
            %         newImgNum(newImgNum(:,1)==curCycleSet(j), 1) = h;
            %     end
            % end
            % curCycles = 1:oCnt;
            % imgNum = newImgNum;
            
            %-------------------------------------------------------%
            
            
            
            
            f = waitbar(0,'Converting Raw Scan Data');
            for k = 1:numel(curCycles) % For each cycle load images and append to image block
                curInds = imgInd(imgNum(:,1)==curCycles(k));
                waitbar(k/numel(curCycles),f,'Converting Raw Scan Data');

                % Temporal image sequences can have any range of
                % frames, simply collect into a shared matrix
                imgBlock = uint16(zeros(512,512,numel(curInds)));
                for m = 1:numel(curInds)
                    curTifFile = strcat(curDir,'\',rawListing(validInds(curInds(m))).name);
                    imgBlock(:,:,m) = uint16(imgaussfilt(imread(curTifFile,'tif'),1));
                end
                
                % Save Processed images for each trial of stimulation channel imaged
                destFile = strcat(procPath,'\Data_PlanarProcessed',num2str(k),'.mat');
                save(destFile,'imgBlock');
            end
            close(f)
        end
        

        %% Perform trial analysis
        trialPath = strcat(pullPath,'\PlanarTrialsG8'); %construct folder for trial specific data
        if(~exist(trialPath, 'dir'))
        
           
            % Trial folder does not exist, process individual trial data and
            % save to trial folder
            fprintf('Processing trial data  ');
            mkdir(trialPath) % if output directory doesn't exist make it
            
            
            
            % run EXCTRACT algorithm on all individual trials
            % Prepare config parameters
            config=[];
            config = get_defaults(config);
            config.use_gpu=0;
            config.avg_cell_radius=4;
            config.trace_output_option='nonneg'; 
            config.num_partitions_x=4;
            config.num_partitions_y=4;
            config.cellfind_min_snr=1; 

            nDiv = 2;
            roiMask = false(512,512,nDiv);
            nPM = 0;
            sumVol = zeros(512,512);
            sumSqVol = zeros(512,512);
            regions = [];
            rN = 0;
            rCutOff = 0;
            for div = 1:nDiv
                s = (div-1)*(100/nDiv) + 1; % assuming 100 trials
                e = div*(100/nDiv);

                M = [];
                for k = s:e % For each cycle load images and append to image block
                    % Save Processed images for each trial of stimulation channel imaged
                    destFile = strcat(procPath,'\Data_PlanarProcessed',num2str(k),'.mat');
                    curBlock = load(destFile,'imgBlock');
                    curBlock = curBlock.imgBlock;


                    % extract metrics over entire sessionm frame by frame
                    for z = 1:size(curBlock,3)
                        nPM = nPM+1;
                        sumVol = sumVol + double(curBlock(:,:,z));
                        sumSqVol = sumSqVol + double((curBlock(:,:,z) .* curBlock(:,:,z)));
                    end

                    if(isempty(M))
                        M = curBlock;
                    else
                        M = cat(3, M, curBlock);
                    end
                end

                %Perform the extraction
                output=extractor(M,config); 
                
                for r = 1:size(output.spatial_weights,3)
                    curReg = output.spatial_weights(:,:,r)>0;
                    stats = regionprops(curReg,'PixelList');
                    if(size(stats(1).PixelList,1)>10)
                        rN = rN+1;
                        regions(rN).PixelList = stats.PixelList;
                    end
                end
                regions = regions';
                if(div==1)
                    rCutOff = numel(regions);
                end

                % sliceMask = max(output.spatial_weights,[],3);
                % bwMask = false(size(sliceMask));
                % if(size(sliceMask,3)>0)
                %     bwMask(sliceMask>0.5)=1;
                %     bwMask = bwareaopen(bwMask,10);
                %     roiMask(:,:,div)=bwMask;
                % end 
            end

            % remove overlapping regions (repeate from divided dataset)
            toRmv = false(numel(regions),1);
            for r = 1:numel(regions)
                extractMask = false(512,512);
                curPixels = regions(r).PixelList;
                for m = 1:size(curPixels,1)
                    extractMask(curPixels(m,2),curPixels(m,1))=1;
                end 
                preSum = sum(extractMask,'all');
                for r2 = r+1:numel(regions)
                    extractMask2 = extractMask;
                    curPixels2 = regions(r2).PixelList;
                    for m = 1:size(curPixels2,1)
                        extractMask2(curPixels2(m,2),curPixels2(m,1))=0;
                    end 

                    if(preSum > 1.75*sum(extractMask2,'all'))
                        toRmv(r2)=1;
                    end
                end
            end
            regions(toRmv) = [];

            % Flatten ROIMask
            extractMask = false(512,512);
            for r = 1:numel(regions)
                curPixels = regions(r).PixelList;
                for m = 1:size(curPixels,1)
                    extractMask(curPixels(m,2),curPixels(m,1))=1;
                end 
            end
            % extractMask = max(roiMask,[],3);

            meanVol = double(sumVol)/nPM;
            varianceAll = (double(sumSqVol)/nPM) - double((meanVol.*meanVol));
            standardDev = sqrt(varianceAll);
            % regions=regionprops(extractMask,'PixelList');

            combineFileName = strcat(trialPath,'\Regions.mat');
            save(combineFileName,'regions');  

            % Plot detected regions
            curFig = figure();
            zProjectedBW = max(extractMask,[],3,'omitnan');
            imshow(zProjectedBW/max(zProjectedBW,[],'all'))
            combineFileName = strcat(trialPath,'\Overall_Seg.fig');
            saveas(curFig,combineFileName);  
            combineFileName = strcat(trialPath,'\Overall_Seg.tif');
            saveas(curFig,combineFileName);
            close(curFig)
            combineFileName = strcat(trialPath,'\Overall_Seg.mat');
            save(combineFileName,'roiMask');  


            % Plot outlines of each detected region
            curFig = figure();
            combinedActive = max(standardDev,[],3)/max(standardDev,[],'all');
            imshow(combinedActive)
            hold on
            for r = 1:numel(regions)
                zProjectedBW = false(512,512);
                curPixels = regions(r).PixelList;
                for m = 1:size(curPixels,1)
                    zProjectedBW(curPixels(m,2),curPixels(m,1))=1;
                end 
                boundaries = bwboundaries(zProjectedBW);
                for k=1:numel(boundaries)
                   b = boundaries{k};
                   plot(b(:,2),b(:,1),'r','LineWidth',1);
                end
            end
            combineFileName = strcat(trialPath,'\Overlay_Seg.fig');
            saveas(curFig,combineFileName);  
            combineFileName = strcat(trialPath,'\Overlay_Seg.tif');
            saveas(curFig,combineFileName);
            close(curFig)


            
            trialTot = 100;
            TrialKeys = cell(trialTot,1);
            TrialInfo = zeros(trialTot,1);
            trialCount = 0;
            deltaAll = zeros(512,512);
            prestimAve = zeros(100,1);
            prestimSTD = zeros(100,1);
            
            % Load timing key data
            dataCellArr = load(strcat(keypathO,'\',keyFile));
            dataCellArr = dataCellArr.dataCellArr;
            eCell = zeros(numel(dataCellArr),1);
            for k = 1:numel(dataCellArr)
                if(isempty(dataCellArr{k}))
                    eCell(k)=1;
                end
            end
            dataCellArr(eCell==1)=[];
            allTraces = cell(100,3);

            nPr = 0;
            nM = 0;
            nPo = 0;
            sumPr = zeros(512,512);
            sumM = zeros(512,512);
            sumPo = zeros(512,512);
            sumSqPr = zeros(512,512);
            sumSqM = zeros(512,512);
            sumSqPo = zeros(512,512);
            f = waitbar(0,'Processing Trial Data');

            for trialNum = 1:numel(dataCellArr) % process each trial
                
                % Valid trial, process data. Trials with error flags
                % will be ignored/discarded
                trial_i = dataCellArr{trialNum};
                validTrial = trial_i{7};
                waitbar(trialNum/numel(dataCellArr),f,'Processing Trial Data');
                
                if(validTrial)
                    trialCount = trialCount+1; % increment trial count
                    
                    % load imaged volume set for trial
                    destFile = strcat(procPath,'\Data_PlanarProcessed',num2str(trialNum),'.mat');
                    fileData = load(destFile,'imgBlock');
                    curSeq = fileData.imgBlock;

                    % Construct trial indexing data and data file name
                    TrialInfo(trialNum,1) = trial_i{2}; % Channel #
                    TrialInfo(trialNum,3) = trial_i{1}; % stimulation current
                    TrialInfo(trialNum,4) = trial_i{5}; % Hit(1) of Miss trial (0)
                    TrialInfo(trialNum,5) = trial_i{6}; % Response time
                    TrialInfo(trialNum,6) = trial_i{3}; % trial frequency
                    TrialInfo(trialNum,7) = trial_i{4}; % trial pulse duration

                    % Collect trial specific timing data
                    if(numel(trial_i{10})<2)
                        stim_ts = [seconds(trial_i{10}), seconds(trial_i{10}+0.6)];
                    else
                        stim_ts = seconds(trial_i{10});
                    end
                    imgStarts = seconds(trial_i{9});
                    imgStarts = imgStarts(1:size(curSeq,3)); % Trim Slices to limit of collected sequence
                    stimStart = stim_ts(1);
                    stimEnd = stim_ts(2) + seconds(1); % Add an additional 1000 ms to end of stim period to collect full flourescence activity
                    preSlices = imgStarts<stimStart;
                    postSlices = imgStarts>stimEnd;
                    midSlices = imgStarts>=stimStart;
                    b = imgStarts<=stimEnd;
                    midSlices(b==0)=0;

                    % Assign frames to pre, mid, and post stim regions
                    preFrames = curSeq(:,:,preSlices);
                    midFrames = curSeq(:,:,midSlices);
                    postFrames = curSeq(:,:,postSlices);
                    preStim = mean(double(preFrames),3);
                    midStim = mean(double(midFrames),3);
                    postStim = mean(double(postFrames),3);
                    
                    % Quantify flouresence traces for given trial
                    roiTracesPre = zeros(numel(regions),size(preSlices,3));
                    roiTracesMid = zeros(numel(regions),size(midSlices,3));
                    roiTracesPost = zeros(numel(regions),size(postSlices,3));
                    for r = 1:numel(regions)
                        trialMask = false(512,512);
                        curPixels = regions(r).PixelList;
                        for m = 1:size(curPixels,1)
                            trialMask(curPixels(m,2),curPixels(m,1))=1;
                        end
                        
                        % process pre-stim data
                        for z = 1:size(preFrames,3)
                            cS = preFrames(:,:,z);
                            vals = cS(trialMask);
                            roiTracesPre(r,z) = mean(vals);
                            
                            nPr = nPr+1;
                            sumPr = sumPr + double(cS);
                            sumSqPr = sumSqPr + double(cS .* cS);
                        end
                        
                        % process mid-stim data
                        for z = 1:size(midFrames,3)
                            cS = midFrames(:,:,z);
                            vals = cS(trialMask);
                            roiTracesMid(r,z) = mean(vals);
                            
                            nM = nM+1;
                            sumM = sumM + double(cS);
                            sumSqM = sumSqM + double(cS .* cS);
                        end
                        
                        % process post-stim data
                        for z = 1:size(postFrames,3)
                            cS = postFrames(:,:,z);
                            vals = cS(trialMask);
                            roiTracesPost(r,z) = mean(vals);
                            
                            nPo = nPo+1;
                            sumPo = sumPo + double(cS);
                            sumSqPo = sumSqPo + double(cS .* cS);
                        end
                    end
                    allTraces{trialNum,1}=roiTracesPre;
                    allTraces{trialNum,2}=roiTracesMid;
                    allTraces{trialNum,3}=roiTracesPost;
                    
                    TrialKeys{trialNum} = strcat(trialPath,'\t',num2str(trialNum),'.mat'); % trial key filename
                    save(TrialKeys{trialNum},'preStim','midStim','postStim','preFrames','preFrames','postStim','roiTracesPre','roiTracesMid','roiTracesPost');

                    % Add to running metrics
                    prestimAve(trialCount) = mean(preStim,'all','omitnan');
                    prestimSTD(trialCount) = std(preStim,[],'all','omitnan');
                    deltaAll = deltaAll + double(midStim)-double(preStim);
                end
            end
            prestimAve = prestimAve(1:trialCount);
            prestimSTD = prestimSTD(1:trialCount);
            
            % Calculate session Metrics
            meanPre = double(sumPr)/nPr;
            variancePre = (double(sumSqPr)/nPr) - double((meanPre.*meanPre));
            stdPre = sqrt(variancePre);

            meanMid = double(sumM)/nM;
            varianceMid = (double(sumSqM)/nM) - double((meanMid.*meanMid));
            stdMid = sqrt(varianceMid);

            meanPost = double(sumPo)/nPo;
            variancePost = (double(sumSqPo)/nPo) - double((meanPost.*meanPost));
            stdPost = sqrt(variancePost);

            
            % calculate baselines for each ROI during pre, mid, & post
            % periods
            preReg = zeros(numel(regions),2);
            midReg = zeros(numel(regions),2);
            postReg = zeros(numel(regions),2);
            regPr = [];
            regM = [];
            regPo = [];
            for t = 1:100 % collect data across trails
               regPr = [regPr, allTraces{t,1}];
               regM = [regM, allTraces{t,2}];
               regPo = [regPo, allTraces{t,3}];
            end
            for r = 1:numel(regions)
                preReg(r,:) = [mean(regPr(r,:)), std(regPr(r,:))];
                midReg(r,:) = [mean(regM(r,:)), std(regM(r,:))];
                postReg(r,:) = [mean(regPo(r,:)), std(regPo(r,:))];
            end
            
            % Save session wide metrics
            destPath = strcat(trialPath,'\SessionMetrics.mat'); % trial key filename
            save(destPath,'meanPre','variancePre','stdPre','meanMid','varianceMid','stdMid',...
                'meanPost','variancePost','stdPost','prestimAve','prestimSTD',...
                'preReg','regM','regPo','allTraces','regions','dataCellArr','TrialInfo');
            
            % Save z-projection of baseline imageBlock standard deviation
            curFig = figure();
            zProjected = max(stdPre,[],3,'omitnan');
            imshow(zProjected/max(zProjected,[],'all'))
            combineFileName = strcat(trialPath,'\baselineSTD_Overall.fig');
            saveas(curFig,combineFileName);  
            combineFileName = strcat(trialPath,'\baselineSTD_Overall.tif');
            saveas(curFig,combineFileName);
            close(curFig)

            curFig = figure();
            zProjected = max(deltaAll,[],3,'omitnan');
            imshow(zProjected/max(zProjected,[],'all'))
            combineFileName = strcat(trialPath,'\delta_Overall.fig');
            saveas(curFig,combineFileName);  
            combineFileName = strcat(trialPath,'\delta_Overall.tif');
            saveas(curFig,combineFileName);
            close(curFig)

            % save plot of intensity across session
            curFig = figure();
            errorbar(prestimAve,prestimSTD)
            title('Average Trial Intensity')
            combineFileName = strcat(trialPath,'\trialPreStimIntensity.fig');
            saveas(curFig,combineFileName);  
            combineFileName = strcat(trialPath,'\trialPreStimIntensity.tif');
            saveas(curFig,combineFileName);
            close(curFig)
            close(f)
        end
        
        
        
        
        %% Perform secondary processing aggrigating trials based on channel #, current, hit status, or combinations therein
        fprintf('\nSeconday processing initiated\n');  
        figureFlag=0;
        
        % Create output folder for storing results        
        selpathOut = strcat(pullPath,'\PlanarAnalysis_G8_',datestr(datetime('now'),'mm-dd-yyyy_HH-MM'));
        mkdir(selpathOut); %Create output folder
        
        % load data
        destPath = strcat(trialPath,'\SessionMetrics.mat'); % trial key filename
        sessMetrics = load(destPath);
        
        % Load metrics from session
        allTraces = sessMetrics.allTraces;
        allSpikes = cell(100,3); % genrate and fill spiking dataset
        for t = 1:100
            for g = 1:3
            	allSpikes{t,g} = false(size(allTraces{t,g}));
            end
        end
        
        regions = sessMetrics.regions;
        dataCellArr = sessMetrics.dataCellArr;
        TrialInfo = sessMetrics.TrialInfo;
        validTrials = false(100,1);
        validTrials = TrialInfo(:,6)~=0;
%         for t = 1:100
%             trial_i = dataCellArr{t};
%             validTrials(t) = trial_i{7};
%         end
        
        
        % Calculate mean value for each ROI
        regionMeans = zeros(numel(regions),1);
        regionSTDs = zeros(numel(regions),1);
        spikes = cell(numel(regions),3);
        for r = 1:numel(regions)
            trVals = [];
            for t = 1:100
                if(validTrials(t)) % process only valid trials
                    pr = allTraces{t,1};
                    mi = allTraces{t,2};
                    po = allTraces{t,3};
                    trVals = [trVals,pr(r,:),mi(r,:),po(r,:)];
                end
            end
            regionMeans(r) = mean(trVals);
            regionSTDs(r) = std(trVals);
            
            
            % Perform deconvolution
            curTrace = smooth(trVals,5);
            curTrace = (curTrace-mean(curTrace))./mean(curTrace);
            g = 0.95;
            [~,temp , ~] = deconvolveCa(curTrace, 'ar1', g,'thresholded', 'optimize_pars', 'thresh_factor', 0.95);
            cSpikes = temp > 0.1;
            
%             if(sum(cSpikes)>0)
%                 sum(cSpikes)
%             end
            
            sCnt = 0;
            for t = 1:100
                if(validTrials(t)) % process only valid trials
                    pr = allTraces{t,1};
                    mi = allTraces{t,2};
                    po = allTraces{t,3};
                    
                    spr =allSpikes{t,1};
                    smi =allSpikes{t,2};
                    spo =allSpikes{t,3};
                    spr(r,:) = cSpikes(1+sCnt : sCnt+size(pr,2));
                    smi(r,:) = cSpikes(1+sCnt+size(pr,2) : sCnt+size(pr,2)+size(mi,2));
                    spo(r,:) = cSpikes(1+sCnt+size(pr,2)+size(mi,2) : sCnt+size(pr,2)+size(mi,2)+size(po,2));
                    allSpikes{t,1} = spr;
                    allSpikes{t,2} = smi;
                    allSpikes{t,3} = spo;
                    
                    sCnt = sCnt+size(pr,2)+size(mi,2)+size(po,2);
                end
            end
        end
        
        
        
        % Sweep through all stimulation parameters, create folder to store
        % overall results.  Individual results will be stored in secondary
        % folders.
        sweepResultPath = strcat(selpathOut,'\OverallResults');
        mkdir(sweepResultPath); %Create output folder
        
        
        uChans = unique(TrialInfo(:,1));
        uCurrs = unique(TrialInfo(:,3));
        uChans(uChans==0)=[];% Remove catch trials from those to process
        uCurrs(uCurrs==0)=[];
        trialMaskHvM = TrialInfo(100,4);
        
        % Calculate distance of each ROI to the Stimulating electrodes
%         chanDistances = zeros(numel(uChans),numel(regions));
        %TODO - add comparision to mapped electrode
        
        allRasterH = false(numel(uChans),numel(uCurrs),numel(regions));
        allRasterM = false(numel(uChans),numel(uCurrs),numel(regions));
        allDataH = cell(numel(uChans),numel(uCurrs));
        allDataM = cell(numel(uChans),numel(uCurrs));
        % Analyze the hits, misses, and currents for each channel
        for i = 1:numel(uChans)
            disp(strcat('  Processing channel: ',num2str(uChans(i))));        
            for ci = 1:numel(uCurrs)

                % Analyze each stimulation parameter, relatively
                % quick to perform with flourescence trace data

                % Find indicies for trial with current parameters,
                % if they exist
                trialMask = false(size(TrialInfo,1),1);
                trialMask(TrialInfo(:,1) == uChans(i))=1; % channel
                trialMask(TrialInfo(:,3) ~= uCurrs(ci))=0; % stimulation current
                trialMask(TrialInfo(:,6) ~= 100)=0; % stimulation frequency (just examine standard)
                trialMask(TrialInfo(:,7) ~= 167)=0; % stimulation pulse duration  (just examine standard)

                trHitMask = trialMask; % Hit trials
                trHitMask(trialMaskHvM==0)=0;
                trHitInds = find(trHitMask);
                curFileNameHit = strcat(sweepResultPath,'\CH',num2str(uChans(i)),'-CURR',num2str(uCurrs(ci)),'-FREQ100-PULSE186-HIT\');
                if(~isempty(trHitInds))
                    % Process Hits for current stimulation parameters
                    mkdir(curFileNameHit); %Create output folder
                    [hRaster,hData] = SecondaryProcess(allTraces,curFileNameHit,trHitInds,'HIT',uCurrs(ci),uChans(i),regions,validTrials,regionMeans,regionSTDs,figureFlag,allSpikes);
                    allRasterH(i,ci,:) = hRaster;
                    allDataH{i,ci} = hData;
                end


%                 trMissMask = trialMask; % Miss trials
%                 trMissMask(trialMaskHvM==1)=0;
%                 trMissInds = find(trMissMask);
%                 curFileNameMiss = strcat(sweepResultPath,'\CH',num2str(uChans(i)),'-CURR',num2str(uCurrs(ci)),'-FREQ100-PULSE167-MISS\');
%                 if(~isempty(trMissInds))
%                     % Process Misses for current stimulation parameters
%                     mkdir(curFileNameMiss); %Create output folder
%                     [mRaster,mData] = SecondaryProcess(allTraces,curFileNameMiss,trMissInds,'MISS',uCurrs(ci),uChans(i),regions,validTrials,regionMeans,regionSTDs,figureFlag);
%                     allRasterM(i,ci,:) = mRaster;
%                     allDataM{i,ci} = mData;
%                 end
%             end
            
                if(figureFlag==1) 

                    % Plot traces from all ROIs for each channel
                    if(~isempty(trHitInds))
                        curFig = figure(); 
                        hold on;
                        sTimes = [];
                        for r = 1:numel(regions) % Process each active region

                            allTrace = [];
                            if(r==1)
                                sTimes = zeros(numel(trHitInds),1);
                            end
                            for t = 1:numel(trHitInds)
                                if(validTrials(trHitInds(t)))
                                    % load trial data
                                    curPre = allTraces{trHitInds(t),1};
                                    curPre = curPre(r,:);
                                    curMid = allTraces{trHitInds(t),2};
                                    curMid = curMid(r,:);
                                    curPost = allTraces{trHitInds(t),3};
                                    curPost = curPost(r,:);

                                    % trim pre (last 500 ms) and post (first 800)
                                    trimPre = curPre(numel(curPre)-15:end);
                                    % mid stimulation should be a fixed
                                    % duration -- double check

                                    % String parts together
                                    trTrace = [trimPre,curMid,curPost];
                                    dff = (trTrace-regionMeans(r))./regionMeans(r); % calculate regional DF/F
                                    dff = smooth(dff,5);

                                    sTimes(t) = numel(trimPre) + numel(allTrace);
                                    allTrace = [allTrace;dff;NaN(15,1)];
                                end
                            end

                            ms = 1:numel(allTrace);

                            if(hRaster(r)==1)
                                plot(ms*33,allTrace+r*1.5,'red','linewidth',1);
                            else
                                plot(ms*33,allTrace+r*1.5,'blue','linewidth',1);
                            end

                        end

                        xlabel('Time (ms)')
                        ylabel('Regions (\DeltaF/F)')
                        for s = 1:numel(sTimes)
                            xline(sTimes(s)*33);
                        end
                        tText = 'Active region traces';
                        title(tText)
                        combineFileName = strcat(curFileNameHit,'CH',num2str(uChans(i)),'-CUR',num2str(uCurrs(ci)),'ROI-Traces.fig');
                        saveas(curFig,combineFileName);  
                        combineFileName = strcat(curFileNameHit,'CH',num2str(uChans(i)),'-CUR',num2str(uCurrs(ci)),'ROI-Traces.tif');
                        saveas(curFig,combineFileName);
                        close(curFig)
                    end
                end
            end
        end
        
        % Perform cross metric analyses of planar temporal responses within
        % this session acreoss channels
        for keyVal = 1:7
            PlotMetricsVsCurrents(allDataH,keyVal,uChans,uCurrs,allRasterH,sweepResultPath)
        end

        % Save overall analysis metrics from this dataset
        combineFileName = strcat(sweepResultPath,'\SessionMetrics.mat');
        save(combineFileName,'allRasterH','allRasterM','allDataH','allDataM','regions','allTraces','TrialInfo','validTrials','uChans','uCurrs'); 
    end
end


function [] = PlotMetricsVsCurrents(data,keyVal,uChans,uCurrs,allRaster,sweepResultPath)
    % Selected metric vs Weeks 
    curFig = figure();
    hold on
    
    metricVals = NaN(numel(uCurrs),numel(uChans));
    metricSTDs = NaN(numel(uCurrs),numel(uChans));
    metricCurrVals = cell(numel(uCurrs),numel(uChans));

    for ci = 1:numel(uCurrs) % process each current
        for i = 1:numel(uChans) % process each channel
            tRaster = allRaster(i,ci,:);
            curData =  data{i,ci};
            if(sum(tRaster)>0)
                % group data from active ROIs
                switch keyVal
                    case 1
                        curVals = curData.spikeStart10Times;
                    case 2
                        curVals = curData.spikePeaksTimes;
                    case 3
                        curVals = curData.spikeVol10;
                    case 4
                        curVals = curData.spikeVolSTD;
                    case 5
                        curVals = curData.spikeDuration10;
                    case 6
                        curVals = curData.spikeDurationSTD;
                    case 7
                        curVals = curData.spikeMax;
                end

                tVals = curVals(tRaster);
                metricVals(ci,i) = mean(tVals);
                metricSTDs(ci,i) = std(tVals);
                temp = metricCurrVals{ci};
                if(isempty(temp))
                    temp= tVals;
                else
                    temp = [temp;tVals];
                end
                metricCurrVals{ci} = temp;
            end
        end
    end 
%     curAves = mean(metricVals,2,'omitnan');
%     plot(uCurrs,curAves)
    if(keyVal==7)
        for ci = 1:numel(uCurrs)
            subplot(numel(uCurrs),1,ci)
            histogram(metricCurrVals{ci},0:0.1:3)
            set(gca, 'YScale', 'log')
        end
    else
        errorbar(uCurrs,mean(metricVals,2,'omitnan'),mean(metricSTDs,2,'omitnan'))
    end
    
    switch keyVal
        case 1
            title('Average Spike Start Time [10% Peak] over Currents')
            sTxt = 'SpikeStart10';
            ylabel('Time Post Stimulation Onset (ms)')
            xlabel('Stimulating Current (\muA)')
        case 2
            title('Average Spike Start TIme [STD Threshold] over Currents')
            sTxt = 'SpikeStartThresh';
            ylabel('Time Post Stimulation Onset (ms)')
            xlabel('Stimulating Current (\muA)')
        case 3
            title('Average Spike Volume [10% Peak] over Currents')
            sTxt = 'SpikeVolume10';
            ylabel('Spike Volume (ms)')
            xlabel('Stimulating Current (\muA)')
        case 4
            title('Average Spike Volume [STD Threshold] over Currents')
            sTxt = 'SpikeVolumeThresh';
            ylabel('Spike Volume (ms)')
            xlabel('Stimulating Current (\muA)')
        case 5
            title('Average Spike Duration [10% Peak] over Currents')
            sTxt = 'SpikeDur10';
            ylabel('Spike Duration (ms)')
            xlabel('Stimulating Current (\muA)')

        case 6
            title('Average Spike Duration [STD Threshold] over Currents')
            sTxt = 'SpikeDurThresh';
            ylabel('Spike Duration (ms)')
            xlabel('Stimulating Current (\muA)')

        case 7
            title('Peak flourescense Histogram')
            sTxt = 'PeakFlour';
%             ylabel('\DeltaF/F')
            ylabel('Count')
            xlabel('\DeltaF/F')
    end
    combineFileName = strcat(sweepResultPath,sTxt,'.fig');
    saveas(curFig,combineFileName);  
    combineFileName = strcat(sweepResultPath,sTxt,'.tif');
    saveas(curFig,combineFileName);
    close(curFig)
end
















function [raster,data] = SecondaryProcess(allTraces,fileNameApp,tInds,HS,curr,curCh,regions,validTrials,regionMeans,regionSTDs,figureFlag,spikes)
    % Process data from trials included in selection set
    raster = false(numel(regions),1);
    spikeStart10Times = zeros(numel(regions),1);
    spikePeaksTimes = zeros(numel(regions),1);
    spikeVol10 = zeros(numel(regions),1);
    spikeVolSTD = zeros(numel(regions),1);
    spikeDuration10 = zeros(numel(regions),1);
    spikeDurationSTD = zeros(numel(regions),1);
    spikeMax = zeros(numel(regions),1);
    deconSpikeStartTimes = zeros(numel(regions),1);
    deconSpikeCnt = zeros(numel(regions),1);
    
    
    for r = 1:numel(regions) % Process each region
        sThresh = regionSTDs(r)/regionMeans(r);
        
        % analyze aligned mid traces overlapping across trials
        aveTrace = NaN(numel(tInds),10);
        aveTraceSpikes = NaN(numel(tInds),10);
        for t = 1:numel(tInds)
            if(validTrials(tInds(t)))
                % load trial data
                curPre = allTraces{tInds(t),1};
                curPre = curPre(r,:);
                curMid = allTraces{tInds(t),2};
                curMid = curMid(r,:);
                curPost = allTraces{tInds(t),3};
                curPost = curPost(r,:);
                
                curPreS = spikes{tInds(t),1};
                curPreS = curPreS(r,:);
                curMidS = spikes{tInds(t),2};
                curMidS = curMidS(r,:);
                curPostS = spikes{tInds(t),3};
                curPostS = curPostS(r,:);
                

                % trim pre (last 500 ms) and post (first 800)
                trimPre = curPre(numel(curPre)-15:end);
                trimPreS = curPreS(numel(curPreS)-15:end);
                % mid stimulation should be a fixed
                % duration -- double check

                % String parts together
                trTrace = [trimPre,curMid,curPost];
                dff = (trTrace-regionMeans(r))./regionMeans(r); % calculate regional DF/F
                dff = smooth(dff,5);

                aveTrace(t,1:numel(dff)) = dff(1:numel(dff));
                
                trTraceSpikes = [trimPreS,curMidS,curPostS];
                aveTraceSpikes(t,1:numel(trTraceSpikes)) = trTraceSpikes;
            end
        end

        % calculate average spiking
        grpSpikesAve = mean(double(aveTraceSpikes),1);
        grpSpikes = grpSpikesAve>=0.5;
        
        % calculate  average trace
        aveTrace(aveTrace==0)=NaN;
        grpTrace = mean(aveTrace,1);
        
        % Determine if average trace corresponds to active neuron or not,
        % simple thresholding 
        grpTraceMidPost = grpTrace(16:42); % examine stimulation region and a 100 ms after 
        if(max(grpTraceMidPost)>sThresh)
            raster(r)=1;
            
            % Calculate Spiking metrics for active regions (all
            % calculations made with DF/F
            [maxDFF, mInd] = max(grpTrace);
            thresh10 = maxDFF/10; % 10% of maximum stimulation floor
            
            % find time points above 10% threshold but after pre-stim
            % period
            aInds = find(grpTrace>thresh10);
            aInds(aInds<16)=[]; % Remove any before stimulation
            ind = find([true, aInds(2:end) < aInds(1:end-1), true]);
            len = diff(ind);
            [~, p] = max(len);
            aSeg  = aInds(ind(p):ind(p + 1) - 1); % Isolatee the longest continous segment of activation and quantify
            
            % Sum volume of DF/F in 10% threshold activated case
            vol10 = 0;
            for i = 1:numel(aSeg)
                vol10 = vol10 + grpTrace(aSeg(i));
            end
            
            % Find time points above normalized region STD
            sInds = find(grpTrace>sThresh);
            sInds(sInds<16)=[]; % Remove any before stimulation
            ind = find([true, sInds(2:end) < sInds(1:end-1), true]);
            len = diff(ind);
            [~, p] = max(len);
            sSeg  = sInds(ind(p):ind(p + 1) - 1); % Isolatee the longest continous segment of activation and quantify
            
            % Sum volume of DF/F in STD threshold case
            volSTD = 0;
            for i = 1:numel(sSeg)
                volSTD = volSTD + grpTrace(sSeg(i));
            end
            
            spikeStart10Times(r) = (aInds(1)-16)*33; % Time from stim to an increase of 10% of maximum
            spikePeaksTimes(r) = (mInd-16)*33; % Time to maximum flourescence (in ms)
            spikeVol10(r) = vol10;
            spikeVolSTD(r) = volSTD;
            spikeDuration10(r) = numel(aSeg)*33; % DF/F over 10% (in ms)
            spikeDurationSTD(r) = numel(sSeg)*33; % DF/F over 10% (in ms)
            spikeMax(r) = maxDFF;
            sTime = find(grpSpikes,1,'first')*33;
            if(isempty(sTime))
                deconSpikeStartTimes(r) = NaN;
            else
                deconSpikeStartTimes(r) = sTime;
            end
            deconSpikeCnt(r) = sum(grpSpikes);
        end
    end
        
    if(figureFlag==1) 
        for r = 1:numel(regions) % Process each active region
            if(raster(r)==1)
                % Plot aligned traces overlapping across
                % trails
                txc = strcat('Ch',num2str(curCh),{' '},'ROI',num2str(r),{' '},'Raw Traces',{' '},HS,{' '},'Curr',{' '},num2str(curr));
                curFig = figure('NumberTitle', 'off', 'Name', txc{1},'visible','off'); 
                hold on;
                aveTrace = NaN(numel(tInds),10); % equates to a 2000 ms trace
                for t = 1:numel(tInds)
                    if(validTrials(tInds(t)))
                        % load trial data
                        curPre = allTraces{tInds(t),1};
                        curPre = curPre(r,:);
                        curMid = allTraces{tInds(t),2};
                        curMid = curMid(r,:);
                        curPost = allTraces{tInds(t),3};
                        curPost = curPost(r,:);

                        % trim pre (last 500 ms) and post (first 800)
                        trimPre = curPre(numel(curPre)-15:end);
                        % mid stimulation should be a fixed
                        % duration -- double check

                        % String parts together
                        trTrace = [trimPre,curMid,curPost];
                        dff = (trTrace-regionMeans(r))./regionMeans(r); % calculate regional DF/F
                        dff = smooth(dff,5);
                        ms = -numel(trimPre):(numel(dff)-numel(trimPre)-1);
                        p1 = plot(ms*33,dff,'linewidth',1);
                        p1.Color = [0 0 1 0.25];

                        aveTrace(t,1:numel(dff)) = dff(1:numel(dff));
                    end
                end

                % Plot  average trace
                aveTrace(aveTrace==0)=NaN;
                grpTrace = mean(aveTrace,1);
                ms = -15:size(aveTrace,2)-16;
                p1 = plot(ms*33,grpTrace,'linewidth',2);
                p1.Color = [0 0 1];

                xline(0);
                xline(700);
                ylabel('\DeltaF/F')
                xlabel('Time (\mus)')
                yline(((regionSTDs(r))/regionMeans(r)),'--');

                combineFileName = strcat(fileNameApp,'RawTrace',num2str(r),'.tif');
                saveas(curFig,combineFileName);
                close(curFig);
            end
        end
    end
        
    data.spikeStart10Times = spikeStart10Times;
    data.spikePeaksTimes = spikePeaksTimes;
    data.spikeVol10 = spikeVol10;
    data.spikeVolSTD = spikeVolSTD;
    data.spikeDuration10 = spikeDuration10;
    data.spikeDurationSTD = spikeDurationSTD;
    data.spikeMax = spikeMax;
    data.deconStart = deconSpikeStartTimes;
    data.deconCount = deconSpikeCnt;
    
    
        
%         
%         % Plot all aligned traces overlapping across
%         % trails (Very Slow)
%         txc = strcat('Ch',num2str(curCh),{' '},'ROI',num2str(r),{' '},'Raw Traces',{' '},HS,{' '},'Curr',{' '},num2str(curr));
%         curFig = figure('NumberTitle', 'off', 'Name', txc{1},'visible','off'); 
%         hold on;
%         aveTrace = NaN(numel(tInds),62); % equates to a 2000 ms trace
%         for t = 1:numel(tInds)
%             if(validTrials(t))
%                 % load trial data
%                 curPre = allTraces{tInds(t),1};
%                 curPre = curPre(r,:);
%                 curMid = allTraces{tInds(t),2};
%                 curMid = curMid(r,:);
%                 curPost = allTraces{tInds(t),3};
%                 curPost = curPost(r,:);
% 
%                 % trim pre (last 500 ms) and post (first 800)
%                 trimPre = curPre(numel(curPre)-15:end);
%                 % mid stimulation should be a fixed
%                 % duration -- double check
% 
%                 % String parts together
%                 trTrace = [trimPre,curMid,curPost];
%                 dff = (trTrace-regionMeans(r))./regionMeans(r); % calculate regional DF/F
%                 dff = smooth(dff,5);
%                 ms = -numel(trimPre):(numel(dff)-numel(trimPre)-1);
%                 p1 = plot(ms*33,dff,'linewidth',1);
%                 p1.Color = [0 0 1 0.25];
% 
%                 aveTrace(t,1:numel(dff)) = dff(1:numel(dff));
%                 if(numel(dff)<size(aveTrace,2))
%                     % Extra spaces following current trace, add padding NaNs
%                     aveTrace(t,numel(dff)+1:end)=NaN; 
%                 end
%             end
%         end
% 
%         % Plot  average trace
%         grpTrace = mean(aveTrace,1);
%         ms = -15:size(aveTrace,2)-16;
%         p1 = plot(ms*33,grpTrace,'linewidth',2);
%         p1.Color = [0 0 1];
% 
%         xline(0);
%         xline(700);
%         ylabel('\DeltaF/F')
%         xlabel('Time (\mus)')
%         
%         combineFileName = strcat(fileNameApp,'RawTrace',num2str(r),'.tif');
%         saveas(curFig,combineFileName);
%         close(curFig);
%         
%         % Determine if average trace corresponds to active neuron or not,
%         % simple thresholding 
%         if(max(grpTrace)>regionSTDs(r))
%             raster(r)=1;
%         end
%         
        
        
        
        
    
                            
end










