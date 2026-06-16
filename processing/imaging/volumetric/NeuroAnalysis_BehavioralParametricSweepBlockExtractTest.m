function [] = NeuroAnalysis_BehavioralParametricSweepBlockExtractTest(sourceFolders)
% Program which processes the data aquired during behavioral paired imaging
% and stimulation experiments.
   
    threshOverride = 0;
  
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
        
        keypathO = strcat(pullPath,'\KEY');
        segListing = dir(selpath);
        segListing(1:2) = []; % remove directory entries


        % Check if electrode posistion file exisits
        elecPosFile = strcat(pullPath,'\ElecPos.mat');
        if(isfile(elecPosFile))
            distData = load(elecPosFile);
            elecPos = distData.contactPos;
        end

        % Extract key files with timing of all stimulation and  imaging 
        keyListing = dir(keypathO);
        keyListing(1:2) = []; % remove directory entries
        keyNames = zeros(numel(keyListing),2);
        for i = 1:numel(keyListing)
            temp = keyListing(i).name;
            C = strsplit(temp,'_');
            bInd = contains(C,'block');
            keyNames(i,1) = str2double(erase(C{bInd},'block')); % Extract the block index
        end
        [~,keyInds] = sortrows(keyNames);


        %% Process raw images into individual stacks 
        % Only perform this 
        procPath = strcat(pullPath,'\PROC2');
        if(~exist(procPath, 'dir'))
            mkdir(procPath) % if directory doesn't exist make it

            fprintf('Performing initial file conversion  ');
            strCR= -1;
            
            % Rearrange raw folders into temporal order - Key files must have Rep and
            % then an Underscore in their name
            rawScanCh = zeros(numel(segListing),1);
            for i = 1:numel(segListing)
                
                filep = strsplit(segListing(i).name,'_');
                bInd = contains(filep,'block');       
                filep2 = strsplit(filep{bInd},'-');
                rawScanCh(i) = str2double(erase(filep2{1},'block')); % Extract the block index
                
            end
            [ChanRep,order] = sort(rawScanCh);
            segListing = segListing(order); % Now the seg listing should be ordered temporally

            destFile = strcat(procPath,'\DataOrder.mat');
            save(destFile,'ChanRep');
            ssNum=0;
            
            % Open each raw folder
            for i = 1:numel(segListing)
                
                % Obtain file listing for current channel
                curDir = strcat(selpath,'\',segListing(i).name);
                rawListing = dir(curDir);
                rawNames = cell(numel(rawListing),1);
                for k = 1:numel(rawListing)
                    rawNames{k} = rawListing(k).name;
                end
                % Discard files not from Ch3 (this excludes 
                validInds = find(contains(rawNames,'Ch3'));

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
                numCycles = unique(imgNum(:,1));

                % The cycle numbers directly ascend during collection, but jump
                % between.
                imgGroup = cell(10,1);
                curGroup = 1;
                curCycleGroup = numCycles(1);
                oldCycle = numCycles(1);          
                for k = 2:numel(numCycles)
                    curCycle = numCycles(k);

                    if(oldCycle ~= (curCycle-1))
                        % Cycle broke, save  previous cycle
                        imgGroup{curGroup} = curCycleGroup;
                        curCycleGroup = curCycle;
                        oldCycle = curCycle;
                        curGroup = curGroup+1;
                    else
                        % cycle continues, add to it
                        curCycleGroup = [curCycleGroup,curCycle];
                        oldCycle = curCycle;
                    end
                end
                %Save last cycle
                imgGroup{curGroup} = curCycleGroup;
                
                
                % Assume first volumetric scan has the correct number of
                % slices. Future stacks cannot be assumed to be the correct
                % size due to ending early or somehow late..... still confused
                % about that one.  Assume it is an additional stack that merged
                % with the prior stack.
                if(ssNum==0)
                    ssNum = size(cyclePos(cyclePos(:,1)==1,:),1);
                end
                
                cpUnit = 100*(1/numel(segListing)); % progress percentage variables
                cpMin = 100*((i-1)/numel(segListing));
                for grp = 1:numel(imgGroup)
                    curCycles = imgGroup{grp};
                    trialImages = cell(numel(curCycles),1);
                    
                    % Display processing status
                    pVal = round(cpUnit*(grp/numel(imgGroup))+cpMin);     
                    percentageOut = [num2str(pVal) '%%'];
                    strOut = [percentageOut repmat(' ',1,10-length(percentageOut)-1)];
                    fprintf([strCR strOut]);
                    strCR = repmat('\b',1,length(strOut)-1);
                    
                
                    cs = 0; % collected stack count (likely larger than # of cycles due to appended cycles)
                    for k = 1:numel(curCycles) % For each cycle load images and append to image block
                        curInds = imgInd(imgNum(:,1)==curCycles(k));
                        cs = cs+1;
                        if(numel(curInds)==ssNum)
                            % One full cycle collected in cycle
                            imgBlock = uint16(zeros(512,512,ssNum));
                            for m = 1:numel(curInds)
                                curTifFile = strcat(curDir,'\',rawListing(validInds(curInds(m))).name);
                                imgBlock(:,:,m) = uint16(imgaussfilt(imread(curTifFile,'tif'),1));
                            end
                            trialImages{cs}=imgBlock;
                            
                            
                        elseif(numel(curInds)>ssNum)
                            % an additional stack was collected in this cycle
                            % process initial volume
                            imgBlock = uint16(zeros(512,512,ssNum));
                            for m = 1:ssNum
                                curTifFile = strcat(curDir,'\',rawListing(validInds(curInds(m))).name);
                                imgBlock(:,:,m) = uint16(imgaussfilt(imread(curTifFile,'tif'),1));
                            end
                            trialImages{cs}=imgBlock;

                            %Process additional merged volume
                            cs = cs+1;
                            imgBlock = uint16(zeros(512,512,ssNum));
                            for m = 1:numel(curInds)-ssNum
                                curTifFile = strcat(curDir,'\',rawListing(validInds(curInds(m))).name);
                                imgBlock(:,:,m) = uint16(imgaussfilt(imread(curTifFile,'tif'),1));
                            end
                            if((numel(curInds)-ssNum)<ssNum)
                                for n = (m+1):ssNum % Fill NaNs for all non existing slices
                                    imgBlock(:,:,m) = NaN(512,512);
                                end
                            end
                            trialImages{cs}=imgBlock;
                            
                            
                        elseif(numel(curInds)<ssNum)
                            % Less than a full stack was collected in cycle
                            imgBlock = uint16(zeros(512,512,numel(ssNum)));
                            for m = 1:numel(curInds) % Fill data that does exist
                                curTifFile = strcat(curDir,'\',rawListing(validInds(curInds(m))).name);
                                imgBlock(:,:,m) = uint16(imgaussfilt(imread(curTifFile,'tif'),1));
                            end
                            for n = (m+1):ssNum % Fill NaNs for all non existing slices
                                imgBlock(:,:,m) = NaN(512,512);
                            end
                            trialImages{cs}=imgBlock;
                        end
                    end
                    trialImgNum = cs;
                    
                    % Save Processed images for each trial of stimulation
                    % channel imaged
                    destFile = strcat(procPath,'\Data_BLK',num2str(ChanRep(i)),'TR',num2str(grp),'.mat');
                    save(destFile,'trialImages','trialImgNum');
                end
            end
        else
            %Load processing key
            destFile = strcat(procPath,'\DataOrder');
            temp = load(destFile,'ChanRep');
            ChanRep = temp.ChanRep;
        end

    
        
        
        
        
        
        
       
        
        %% Process each channel (if not already processed)
        % Process all trials across all channels and stair combinations, store
        % data from each individual trial for post processing following this
        % initial sweep and processing
        trialTot = 10;
        TrialKeys = cell(trialTot,1);
        TrialInfo = zeros(trialTot,1);
        trialCount = 1;
        

         trialPath = strcat(pullPath,'\TRIALS2'); %construct folder for trial specific data
        if(~exist(trialPath, 'dir'))
            % Trial folder does not exist, process individual trial data and
            % save to trial folder
            fprintf('Processing trial data  ');
            strCR = -1;

             % load initial dataset for sizing
             destFile = strcat(procPath,'\Data_BLK',num2str(ChanRep(1)),'TR',num2str(1),'.mat');
             iniStack = load(destFile,'trialImages','trialImgNum');
             iniStack = iniStack.trialImages;
             imgSz = size(iniStack{1});
             ssNum = imgSz(3);
             
            
            mkdir(trialPath) % if directory doesn't exist make it

            

            deltaAll = zeros(imgSz);     
            intensAve = zeros(10,1);
            intensSTD = zeros(10,1);
            sumVol = zeros(imgSz);
            sumSqVol = zeros(imgSz);
            n = 0;
            sumVolPM = zeros(imgSz);
            sumSqVolPM = zeros(imgSz);
            nPM = 0;
            for curBlock = 1:numel(keyListing) % Process each channel/staircase/Block dataset

                % calculate processing status steps
                cpUnit = 100*(1/numel(keyListing)); % progress percentage variables
                cpMin = 100*((curBlock-1)/numel(keyListing));
                
                
               % Load timing key data
                dataCellArr = load(strcat(keypathO,'\',keyListing(keyInds(curBlock)).name));
                dataCellArr = dataCellArr.dataCellArr;
                eCell = zeros(numel(dataCellArr),1);
                for k = 1:numel(dataCellArr)
                    if(isempty(dataCellArr{k}))
                        eCell(k)=1;
                    end
                end
                dataCellArr(eCell==1)=[];

                
                
                


                

                for trialNum = 1:numel(dataCellArr) % process each trial
                    
                    % Display processing status
                    pVal = round(cpUnit*(trialNum/numel(dataCellArr))+cpMin);     
                    percentageOut = [num2str(pVal) '%%'];
                    strOut = [percentageOut repmat(' ',1,10-length(percentageOut)-1)];
                    fprintf([strCR strOut]);
                    strCR = repmat('\b',1,length(strOut)-1);
                    
                    
                    % Valid trial, process data. Trials with error flags
                    % will be ignored/discarded
                    trial_i = dataCellArr{trialNum};
                    validTrial = trial_i{7};
                    
                    if(validTrial)
                        
                        % load imaged volume set for trial
                        destFile = strcat(procPath,'\Data_BLK',num2str(ChanRep(curBlock)),'TR',num2str(trialNum),'.mat');
                        fileData = load(destFile,'trialImages','trialImgNum');
                        numScanVols = fileData.trialImgNum;
                        imageStacks = fileData.trialImages;

                        
                        % Construct trial indexing data and data file name
                        TrialInfo(trialCount,1) = trial_i{2}; % Channel #
                        TrialInfo(trialCount,3) = trial_i{1}; % stimulation current
                        TrialInfo(trialCount,4) = trial_i{5}; % Hit(1) of Miss trial (0)
                        TrialInfo(trialCount,5) = trial_i{6}; % Response time
                        TrialInfo(trialCount,6) = trial_i{3}; % trial frequency
                        TrialInfo(trialCount,7) = trial_i{4}; % trial pulse duration

                        TrialKeys{trialCount} = strcat(trialPath,'\t',num2str(trialCount),'.mat'); % trial key filename

                        % Collect trial specific timing data
                         if(numel(trial_i{10})<2)
                            stim_ts = [seconds(trial_i{10}), seconds(trial_i{10}+0.6)];
                        else
                            stim_ts = seconds(trial_i{10});
                        end
                        imgStarts = seconds(trial_i{9});
                        stimStart = stim_ts(1);
                        stimEnd = stim_ts(2) + seconds(1); % Add an additional 500 ms to end of stim period to collect full flourescence activity
                        preSlices = imgStarts<stimStart;
                        postSlices = imgStarts>stimEnd;
                        midSlices = imgStarts>=stimStart;
                        b = imgStarts<=stimEnd;
                        midSlices(b==0)=0;

                        % Assign frames to pre, mid, and post stim regions
                        preBlocks = zeros(numScanVols,1);
                        for k = 1:numScanVols
                            if(ssNum*k<numel(preSlices))
                                preBlocks(k) = sum(preSlices(1+ssNum*(k-1):ssNum*k));
                            else
                                preBlocks(k) = sum(preSlices(1+ssNum*(k-1):end));
                            end
                        end

                        midBlocks = zeros(numScanVols,1);
                        for k = 1:numScanVols
                            if(ssNum*k<numel(midSlices))
                                midBlocks(k) = sum(midSlices(1+ssNum*(k-1):ssNum*k));
                            else
                                midBlocks(k) = sum(midSlices(1+ssNum*(k-1):end));
                            end
                        end

                        postBlocks = zeros(numScanVols,1);
                        for k = 1:numScanVols
                            if(ssNum*k<numel(postSlices))
                                postBlocks(k) = sum(postSlices(1+ssNum*(k-1):ssNum*k));
                            else
                                postBlocks(k) = sum(postSlices(1+ssNum*(k-1):end));
                            end
                        end
                        
 
                        
                        
                        
                        
                        % Normalize the response of each imaging trial
                        % section to each other section and trial by limiting 
                        % the number of frames averaged for each slice to 
                        % 2 or approximately 800 - 1000 ms for a 13 slice
                        % dataset
                        zCnt = 2*imgSz(3);
                        for k = numScanVols:-1:1
                            cur = preBlocks(k);
                            if(cur<=zCnt)
                                zCnt = zCnt - cur;
                                % dont trim cur count
                            else
                                cur = zCnt;
                                zCnt = 0;
                                % trim cur count
                            end
                            preBlocks(k) =  cur;
                        end

                        zCnt = 2*imgSz(3);
                        for k = 1:numScanVols
                            cur = midBlocks(k);
                            if(cur<=zCnt)
                                zCnt = zCnt - cur;
                                % dont trim cur count
                            else
                                cur = zCnt;
                                zCnt = 0;
                                % trim cur count
                            end
                            midBlocks(k) =  cur;
                        end

                        zCnt = 2*imgSz(3);
                        for k = 1:numScanVols
                            cur = postBlocks(k);
                            if(cur<=zCnt)
                                zCnt = zCnt - cur;
                                % dont trim cur count
                            else
                                cur = zCnt;
                                zCnt = 0;
                                % trim cur count
                            end
                            postBlocks(k) =  cur;
                        end


  
                        
                        % construct image volume matricies for current trial
                        preFrames4D = uint16(zeros(imgSz(1),imgSz(2),imgSz(3),numScanVols));
                        for k = 1:numScanVols
                            if(size(imageStacks{k},3)<size(preFrames4D,3))
                                temp = imageStacks{k};
                                temp(end,end,size(preFrames4D,3))=0;
                                preFrames4D(:,:,:,k) = temp;
                            else
                                preFrames4D(:,:,:,k) = imageStacks{k};
                            end
                        end 
                        midFrames4D = preFrames4D;
                        postFrames4D = preFrames4D;
                       
                        % Trim partial scan volumes - pre
                        initial = 1;
                        for k = 1:numScanVols
                            if(preBlocks(k)>0)
                                if(preBlocks(k)<ssNum)
                                    if(initial==1)
                                        % Trim front of volume
                                        curStack = preFrames4D(:,:,:,k);
                                        curStack(:,:,1:size(curStack,3)-preBlocks(k))=NaN;
                                        preFrames4D(:,:,:,k) = curStack;
                                        initial=0;
                                    else
                                        % Trim end of volume
                                        curStack = preFrames4D(:,:,:,k);
                                        curStack(:,:,preBlocks(k):end)=NaN;
                                        preFrames4D(:,:,:,k) = curStack;
                                    end
                                end
                            else
                                curStack = preFrames4D(:,:,:,k);
                                preFrames4D(:,:,:,k) = NaN(size(curStack));
                            end
                        end
                        

                        
                        % Trim partial scan volumes - mid 
                        initial = 1;
                        for k = 1:numScanVols
                            if(midBlocks(k)>0)
                                if(midBlocks(k)<ssNum)
                                    if(initial==1)
                                        % Trim front of volume
                                        curStack = midFrames4D(:,:,:,k);
                                        curStack(:,:,1:size(curStack,3)-midBlocks(k))=NaN;
                                        midFrames4D(:,:,:,k) = curStack;
                                        initial=0;
                                    else
                                        % Trim end of volume
                                        curStack = midFrames4D(:,:,:,k);
                                        curStack(:,:,midBlocks(k)+1:end)=NaN;
                                        midFrames4D(:,:,:,k) = curStack;
                                    end
                                end
                            else
                                curStack = midFrames4D(:,:,:,k);
                                midFrames4D(:,:,:,k) = NaN(size(curStack));
                            end
                        end


                        

                        % Trim partial scan volumes - post
                        initial=1;
                        for k = 1:numScanVols
                            if(postBlocks(k)>0)
                                if(postBlocks(k)<ssNum)
                                    if(initial==1)
                                        % Trim front of volume
                                        curStack = postFrames4D(:,:,:,k);
                                        curStack(:,:,1:size(curStack,3)-postBlocks(k))=NaN;
                                        postFrames4D(:,:,:,k) = curStack;
                                        initial=0;
                                    else
                                        % Trim end of volume
                                        curStack = postFrames4D(:,:,:,k);
                                        curStack(:,:,postBlocks(k)+1:end)=NaN;
                                        postFrames4D(:,:,:,k) = curStack;
                                    end
                                end
                            else
                                curStack = postFrames4D(:,:,:,k);
                                postFrames4D(:,:,:,k) = NaN(size(curStack));
                            end
                        end


                        
                        trialPreStimVol = nanmean(preFrames4D,4);
                        trialMidStimVol = nanmean(midFrames4D,4);
                        trialPostStimVol = nanmean(postFrames4D,4);
                        preStim = single(trialPreStimVol);
                        midStim = single(trialMidStimVol);
                        postStim = single(trialPostStimVol);

                        save(TrialKeys{trialCount},'preStim','midStim','postStim');
                        
                        intensAve(trialCount) = mean(preStim,'all','omitnan');
                        intensSTD(trialCount) = std(preStim,[],'all','omitnan');
                        n = n+1;
                        sumVol = sumVol + preStim;
                        sumSqVol = sumSqVol + (preStim .* preStim);
                            
                        % across Channel metrics
                        nPM = nPM+1;
                        sumVolPM = sumVolPM + preStim;
                        sumSqVolPM = sumSqVolPM + (preStim .* preStim);

                        nPM = nPM+1;
                        sumVolPM =  sumVolPM + midStim;
                        sumSqVolPM = sumSqVolPM + (midStim .* midStim);
                        
                        % Save z-projections of pre and mid stim volumes (for
                        % debugging)
                        curFig = figure();
                        imgMax = max([max((preStim),[],'all','omitnan'), max((midStim),[],'all','omitnan')]);
                        subplot(1,2,1)
                        zProjected = max((preStim),[],3,'omitnan');
                        imshow(zProjected/imgMax)
                        title('Prestim')
                        xlabel(strcat('Ave Intensity: ',num2str(mean(preStim,'all'))))
                        subplot(1,2,2)
                        zProjected = max((midStim),[],3,'omitnan');
                        imshow(zProjected/imgMax)
                        title('Midstim')
                        xlabel(strcat('Ave Intensity: ',num2str(mean(midStim,'all'))))
                        combineFileName = strcat(trialPath,'\trialVolumes',num2str(trialCount),'.tif');
                        saveas(curFig,combineFileName);
                        close(curFig)
                        
                        

                            
                                    

                        deltaAll = deltaAll + (midStim-preStim);
                        
                        
                       
                                                
                        trialCount = trialCount+1; % increment trial count
                    end
                end
            end
            fprintf('\n')
            

            
            %% iterate through each channel                        
            uChans = unique(TrialInfo(:,1)); % Channel #
            standardDevCh = cell(numel(uChans),1); % each channel must have it's own key, no repeats on a channel for a given day ( break into seperate 
            meanCh = cell(numel(uChans),1);
            stdMod = zeros(numel(uChans),1);
            for cInd = 1:numel(uChans)
                
                
                chan = uChans(cInd);
                validInds = TrialInfo(:,1)==chan;
                nCH = 0;
                sumVolCH = zeros(imgSz);
                sumSqVolCH = zeros(imgSz);
                chanTrialNum = 0;

                % construct variable for analyzing trial by trial metrics 
                preStimMean = zeros(sum(validInds),1);
                midStimMean = zeros(sum(validInds),1);
                preStimSTD = zeros(sum(validInds),1);
                midStimSTD = zeros(sum(validInds),1);
                deltaStimMean = zeros(sum(validInds),1);
                deltaStimSTD = zeros(sum(validInds),1);
                deltaStimMeanTop = zeros(sum(validInds),1);
                
                
                for curInd = 1:size(TrialInfo,1)
                    if(validInds(curInd)==1) % process trials belonging to current channel
                            
                        % load trial data... this'll take a while
                        trialDat = load(TrialKeys{curInd},'preStim','midStim');
                        preStim = trialDat.preStim;
                        midStim = trialDat.midStim;
                        
                        % Calculate standard deviations
                        % channel specific metrics
                        nCH = nCH +1;
                        sumVolCH = sumVolCH + preStim;
                        sumSqVolCH = sumSqVolCH + (preStim .* preStim);
 
                        
                        % calculate trial by trial metrics for each channel
                        chanTrialNum = chanTrialNum + 1;

                        preStimMean(chanTrialNum) = mean(preStim,'all','omitnan');
                        preStimSTD(chanTrialNum) = std(preStim,[],'all','omitnan');
                        midStimMean(chanTrialNum) = mean(midStim,'all','omitnan');
                        midStimSTD(chanTrialNum) = std(midStim,[],'all','omitnan');

                        deltaVol = midStim-preStim;
                        deltaStimMean(chanTrialNum) = mean(deltaVol,'all','omitnan');
                        deltaStimSTD(chanTrialNum) = std(deltaVol,[],'all','omitnan');

                        % calculate intensity of to 5% of delta voxels
                        deltaVolF = reshape(deltaVol,1,[]);
                        n2=round(.05*numel(deltaVolF));
                        [~,idx]=sort(deltaVolF);
                        top=deltaVolF(idx(end-n2+1:end));
                        deltaStimMeanTop(nCH) = mean(top,'all','omitnan');
                    end
                end
                
                
                % Calculate channel specific STD volumes. Will require
                % parsing data following loading all data blocks
                meanVol = sumVolCH/n;
                variance = (sumSqVolCH/n) - (meanVol.*meanVol);
                standardDevCh{cInd} = sqrt(variance);
                meanCh{cInd} = meanVol;

                curFig = figure();
                zProjected = max(standardDevCh{cInd},[],3,'omitnan');
                imshow(zProjected/max(zProjected,[],'all'))
                title(strcat('Ch',num2str(chan),'  Standard Dev : Ave ',num2str(mean(standardDevCh{cInd},'all','omitnan'))))
                combineFileName = strcat(trialPath,'\Ch',num2str(chan),'_STD.fig');
                saveas(curFig,combineFileName);  
                combineFileName = strcat(trialPath,'\Ch',num2str(chan),'_STD.tif');
                saveas(curFig,combineFileName);
                close(curFig)  
                
                
                %trim trial by trial metrics
                preStimMean = preStimMean(1:chanTrialNum);
                midStimMean = midStimMean(1:chanTrialNum);
                preStimSTD = preStimSTD(1:chanTrialNum);
                midStimSTD = midStimSTD(1:chanTrialNum);
                deltaStimMean = deltaStimMean(1:chanTrialNum);
                deltaStimSTD = deltaStimSTD(1:chanTrialNum);
                deltaStimMeanTop = deltaStimMeanTop(1:chanTrialNum);

                % For well segmenting trials the threshold is usually between
                % the delta mean and the mean of the top 5%
                topVal = mean(deltaStimMeanTop,'all','omitnan'); 
                midVal = mean(deltaStimMean,'all','omitnan');
                targetVal = mean([topVal, midVal]);
                stdBaseVal = mean(standardDevCh{cInd},'all');
                stdMod(cInd) = targetVal/stdBaseVal;
                if(threshOverride~=0)
                    stdMod(cInd) = threshOverride;
                end


                % Plot the trial by trial response  
                curFig = figure();
                alpha = 0.1;  
                h = errorbar(preStimMean,preStimSTD,'Color',[1, 0, 0, 0.2]);
                hold on
                set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
                h = errorbar(midStimMean,midStimSTD,'Color',[0, 1, 0, 0.2]);
                set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
                h = errorbar(deltaStimMean,deltaStimSTD,'Color',[0, 0, 1, 0.2]);
                set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', 'ColorData', [h.Line.ColorData(1:3); 255*alpha])
                plot(preStimMean,'Color',[1, 0, 0, 1],'LineWidth',2)
                plot(midStimMean,'Color',[0, 1, 0, 1],'LineWidth',2)
                plot(deltaStimMean,'Color',[0, 0, 1, 1],'LineWidth',2)
                plot(deltaStimMeanTop,'Color',[1 0 1],'LineWidth',2)
                stdBase = stdMod(cInd)*mean(standardDevCh{cInd},'all');
                yline(stdBase,'r--',strcat(num2str(stdMod(cInd)),'x STD Baseline'));
                xlabel('Trial')
                ylabel('Flourescence Intensity')
                legend('','','','PreStim','MidStim','Delta','Delta: Top 5%','','')
                curFig.WindowState = 'maximized';
                title('Individual Trial Intensity')
                combineFileName = strcat(trialPath,'\trialIntensityCh',num2str(chan),'.fig');
                saveas(curFig,combineFileName);  
                combineFileName = strcat(trialPath,'\trialIntensityCH',num2str(chan),'.tif');
                saveas(curFig,combineFileName);
                close(curFig)
            end
        


            
            % calculate mean and variance for entire session
            overallBaselineMean = sumVol/n; % Pre-stim values
            overallBaselineVariance = (sumSqVol/n) - (overallBaselineMean.*overallBaselineMean);
            combineFileName = strcat(trialPath,'\baselineSTD_Overall.mat');
            save(combineFileName,'overallBaselineVariance','overallBaselineMean');  
            standardDev = sqrt(overallBaselineVariance);
            
            meanVolPM = sumVolPM/nPM; % Pre and Mid-stim values
            varianceAllPM = (sumSqVolPM/nPM) - (meanVolPM.*meanVolPM);
            standardDevAllPM = sqrt(varianceAllPM);
            combineFileName = strcat(trialPath,'\baselineSTD_OverallPM.mat');
            save(combineFileName,'varianceAllPM','meanVolPM','standardDevAllPM');  
            varianceAll = varianceAllPM;
            standardDevAll = standardDevAllPM;      
            
            % Save trial keys and infomration datasets
            keypath = strcat(trialPath,'\KeyRing.mat'); % trial key filename
            save(keypath,'TrialKeys','TrialInfo');
            combineFileName = strcat(trialPath,'\trialPreStimIntensity.mat');
            save(combineFileName,'intensAve','intensSTD');   
            
            combineFileName = strcat(trialPath,'\channelSpecificSTD.mat');
            save(combineFileName,'standardDevCh','meanCh','stdMod','deltaAll');  

            

        else
            % If data was previously processed, load trial keys and infomration datasets
            keypath = strcat(trialPath,'\KeyRing.mat');
            temp = load(keypath,'TrialKeys','TrialInfo');
            TrialKeys = temp.TrialKeys;
            TrialInfo = temp.TrialInfo;
            
            % Load STD deviation data (overall and channel specific)
            destFile = strcat(trialPath,'\baselineSTD_Overall.mat');
            fprintf('Preprocessed Data Loading');
            temp = load(destFile,'overallBaselineVariance');
            variance = temp.overallBaselineVariance;
            standardDev = sqrt(variance);
            temp = load(destFile,'overallBaselineMean');
            overallBaselineMean = temp.overallBaselineMean;
            
            
            combineFileName = strcat(trialPath,'\baselineSTD_OverallPM.mat');
            temp = load(combineFileName,'varianceAllPM','meanVolPM','standardDevAllPM');  
            varianceAll = temp.varianceAllPM;
            standardDevAll = temp.standardDevAllPM;
            
            combineFileName = strcat(trialPath,'\channelSpecificSTD.mat');
            temp = load(combineFileName,'standardDevCh','meanCh','stdMod','deltaAll');  
            standardDevCh = temp.standardDevCh;
            meanCh = temp.meanCh;
            stdMod = temp.stdMod;
            deltaAll = temp.deltaAll;
            
            
            % Load pre-processing baseline data
            combineFileName = strcat(trialPath,'\trialPreStimIntensity.mat');
            intenDat = load(combineFileName,'intensAve','intensSTD'); 
            intensAve = intenDat.intensAve;
            intensSTD = intenDat.intensSTD; 
        end
        
        stdMod(stdMod==2)=1;
        
        %% Perform secondary processing aggrigating trials based on channel #, current, hit status, or combinations therein
        fprintf('\nSeconday processing initiated\n');  
        
        % Create output folder for storing results        
        selpathOut = strcat(pullPath,'\BehaveAnalysis_D4_',datestr(datetime('now'),'mm-dd-yyyy_HH-MM'));
        mkdir(selpathOut); %Create output folder

        % Save z-projection of baseline imageBlock standard deviation
        curFig = figure();
        zProjected = max(standardDev,[],3,'omitnan');
        imshow(zProjected/max(zProjected,[],'all'))
        combineFileName = strcat(selpathOut,'\baselineSTD_Overall.fig');
        saveas(curFig,combineFileName);  
        combineFileName = strcat(selpathOut,'\baselineSTD_Overall.tif');
        saveas(curFig,combineFileName);
        close(curFig)

        curFig = figure();
        zProjected = max(overallBaselineMean,[],3,'omitnan');
        imshow(zProjected/max(zProjected,[],'all'))
        combineFileName = strcat(selpathOut,'\baselineMean_Overall.fig');
        saveas(curFig,combineFileName);  
        combineFileName = strcat(selpathOut,'\baselineMean_Overall.tif');
        saveas(curFig,combineFileName);
        close(curFig)
        
        % save image of intensity 
        curFig = figure();
        errorbar(intensAve,intensSTD)
        title('Average Trial Intensity')
        combineFileName = strcat(selpathOut,'\trialPreStimIntensity.fig');
        saveas(curFig,combineFileName);  
        combineFileName = strcat(selpathOut,'\trialPreStimIntensity.tif');
        saveas(curFig,combineFileName);
        close(curFig)
        
        
        
        
        
        
        
        % examine the flourescence reponse of the pre stimulation period
        volTraces = zeros(size(overallBaselineMean,1),size(overallBaselineMean,2),size(overallBaselineMean,3),5*numel(TrialKeys),'uint16');
        for trialNum = 1:numel(TrialKeys)
            trialDat = load(TrialKeys{trialNum});
            preStim = single(trialDat.preStim);
            midStim = single(trialDat.midStim);
            postStim = single(trialDat.postStim);
            curInd = 5*(trialNum-1);
            volTraces(:,:,:,curInd+1) = (preStim);
%             volTraces(:,:,:,curInd+2) = 0.75*(preStim) + 0.25*(midStim);
            volTraces(:,:,:,curInd+2) = 0.5*preStim + 0.5*midStim;
%             volTraces(:,:,:,curInd+4) = 0.25*(preStim) + 0.75*(midStim);
            volTraces(:,:,:,curInd+3) = (midStim);
%             volTraces(:,:,:,curInd+6) = 0.75*(midStim) + 0.25*(postStim);
            volTraces(:,:,:,curInd+4) = 0.5*midStim + 0.5*postStim;
%             volTraces(:,:,:,curInd+7) = 0.25*(midStim) + 0.75*(postStim);
            volTraces(:,:,:,curInd+5) = (postStim);
        end
        
%         provide planar imaging data to Extract algorithm to segment
%         individual cell regions for 
        
        
        config=[];
        config = get_defaults(config);
        
        %Set some important settings
        config.use_gpu=0;
        config.avg_cell_radius=4;
        config.trace_output_option='nonneg'; 
        config.num_partitions_x=1;
        config.num_partitions_y=1;
        config.cellfind_min_snr=20;
        config.max_iter = 2;
        config.verbose=0;
        config.spatial_highpass_cutoff=inf; % no need for highpass filtering
        config.cellfind_filter_type ='none';
%         config.visualize_cellfinding=1;
        


        roiMask = false(size(overallBaselineMean));
        for layer = 1:size(overallBaselineMean,3)
            M = squeeze(volTraces(:,:,layer,:));
            output = extractor(M,config);
            sliceMask = max(output.spatial_weights,[],3);
            bwMask = false(size(sliceMask));
            if(size(sliceMask,3)>0)
                bwMask(sliceMask>0.5)=1;
                bwMask = bwareaopen(bwMask,10);
                roiMask(:,:,layer)=bwMask;
            end 
        end




        
        
        curFig = figure();
        zProjected = max(standardDevAll,[],3,'omitnan');
        imshow(zProjected/max(zProjected,[],'all'))
        combineFileName = strcat(selpathOut,'\STD_Segmented.fig');
        saveas(curFig,combineFileName);  
        combineFileName = strcat(selpathOut,'\STD_Segmented.tif');
        saveas(curFig,combineFileName);
        close(curFig)
        combineFileName = strcat(selpathOut,'\STD_Segmented.mat');
        save(combineFileName,'standardDevAll');    
        
        
        curFig = figure();
        zProjected = max(deltaAll,[],3,'omitnan');
        imshow(zProjected/max(zProjected,[],'all'))
        combineFileName = strcat(selpathOut,'\DeltaAll.fig');
        saveas(curFig,combineFileName);  
        combineFileName = strcat(selpathOut,'\DeltaAll.tif');
        saveas(curFig,combineFileName);
        close(curFig)
        
        
        


        
        
        %% Perform Segmentation of Regions across all channels
        MIJ.start('C:\Users\rjl6\Desktop\fiji-win64\Fiji.app\plugins','C:\Users\rjl6\Desktop\fiji-win64\Fiji.app\plugins',0);
        MIJ.createImage(uint8(squeeze(512*standardDevAll/max(standardDevAll,[],'all')))) % Pass variable to ImageJ
        MIJ.run("3D Iterative Thresholding", "min_vol_pix=5 max_vol_pix=1000 min_threshold=5 min_contrast=2 criteria_method=EDGES threshold_method=STEP segment_results=All value_method=1 starts");
        
        
        % Transfer data from ImageJ back into Matlab
        mijOut = MIJ.getListImages();
        chkStr = '';
        for k = 1:numel(mijOut)
            chkStr = strcat(chkStr,mijOut(k).toCharArray');
        end
        segROIs = false(size(standardDevAll));

        if(contains(chkStr,'draw')) % Check if ouput figure was generated
            MIJ.selectWindow("draw");
            MIJ.run('Measure Stack...', 'channels slices frames order=czt(default)');
            maxChan = max(MIJ.getColumn('Ch'));
            MIJ.run("Clear Results");
            if(isempty(maxChan))
                MIJ.selectWindow("draw");
                curChanVol = MIJ.getCurrentImage();
                segROIs(curChanVol>0) = 1;
            else
                for chanN = 1:maxChan % Process all channels in output segmentation
                    MIJ.run('Duplicate...','duplicate channels=1');
                    curChanVol = MIJ.getCurrentImage();
                    segROIs(curChanVol>0) = 1;
                end
            end
        end
        MIJ.run('Close All'); 
        try
            MIJ.closeAllWindows(); % Close MIJI figures
        catch
        end   
        
        
        
        MIJ.createImage(uint8(squeeze(512*abs(varianceAll)/max(abs(varianceAll),[],'all')))) % Pass variable to ImageJ
        MIJ.run("3D Iterative Thresholding", "min_vol_pix=10 max_vol_pix=100 min_threshold=10 min_contrast=5 criteria_method=EDGES threshold_method=STEP segment_results=All value_method=5 starts");

        % Transfer data from ImageJ back into Matlab
        mijOut = MIJ.getListImages();
        chkStr = '';
        for k = 1:numel(mijOut)
            chkStr = strcat(chkStr,mijOut(k).toCharArray');
        end
        segROIsV = false(size(standardDevAll));

        if(contains(chkStr,'draw')) % Check if ouput figure was generated
            MIJ.selectWindow("draw");
            MIJ.run('Measure Stack...', 'channels slices frames order=czt(default)');
            maxChan = max(MIJ.getColumn('Ch'));
            MIJ.run("Clear Results");
            if(isempty(maxChan))
                MIJ.selectWindow("draw");
                curChanVol = MIJ.getCurrentImage();
                segROIsV(curChanVol>0) = 1;
            else
                for chanN = 1:maxChan % Process all channels in output segmentation
                    MIJ.run('Duplicate...', strcat('duplicate channels=',num2str(chanN)));
                    curChanVol = MIJ.getCurrentImage();
                    segROIsV(curChanVol>0) = 1;
                end
            end
        end
        MIJ.run('Close All'); 
        try
            MIJ.closeAllWindows(); % Close MIJI figures
        catch
        end    
        
        MIJ.createImage(uint8(squeeze(512*abs(deltaAll)/max(abs(deltaAll),[],'all')))) % Pass variable to ImageJ
        MIJ.run("3D Iterative Thresholding", "min_vol_pix=10 max_vol_pix=100 min_threshold=10 min_contrast=5 criteria_method=EDGES threshold_method=STEP segment_results=All value_method=5 starts");

        % Transfer data from ImageJ back into Matlab
        mijOut = MIJ.getListImages();
        chkStr = '';
        for k = 1:numel(mijOut)
            chkStr = strcat(chkStr,mijOut(k).toCharArray');
        end
        segROIsD = false(size(standardDevAll));

        if(contains(chkStr,'draw')) % Check if ouput figure was generated
            MIJ.selectWindow("draw");
            MIJ.run('Measure Stack...', 'channels slices frames order=czt(default)');
            maxChan = max(MIJ.getColumn('Ch'));
            MIJ.run("Clear Results");
            if(isempty(maxChan))
                MIJ.selectWindow("draw");
                curChanVol = MIJ.getCurrentImage();
                segROIsD(curChanVol>0) = 1;
            else
                for chanN = 1:maxChan % Process all channels in output segmentation
                    MIJ.run('Duplicate...', strcat('duplicate channels=',num2str(chanN)));
                    curChanVol = MIJ.getCurrentImage();
                    segROIsD(curChanVol>0) = 1;
                end
            end
        end
        MIJ.run('Close All'); 
        try
            MIJ.closeAllWindows(); % Close MIJI figures
        catch
        end   
        % When finished processing source data close MIJI extension
        try
            MIJ.exit % close the ImageJ interface
        catch
            disp('MIJI Close Error');
        end
        
        
        % Collect all regions detected
        segROIs(segROIsV==1)=1; % merge two different segmentations
        segROIs(segROIsD==1)=1;
        SE = strel("disk",2);
        segROIs(imdilate(roiMask,SE)==1)=0; % mask out the initial segmentation approach, it takes precidence
        for k = 1:size(segROIs,3)
            segROIs(:,:,k) = bwareaopen(segROIs(:,:,k),10); % Remove small regions, likely noise.
        end
        roiMask(segROIs==1)=1;

        conn = false(3,3,3);
        conn(:,:,2)=1;
        CC = bwconncomp(roiMask,conn);
        regions=regionprops(CC,'PixelList');

        % Plot detected regions
        curFig = figure();
        zProjectedBW = max(roiMask,[],3,'omitnan');
        imshow(zProjectedBW/max(zProjectedBW,[],'all'))
        combineFileName = strcat(selpathOut,'\Overall_Seg.fig');
        saveas(curFig,combineFileName);  
        combineFileName = strcat(selpathOut,'\Overall_Seg.tif');
        saveas(curFig,combineFileName);
        close(curFig)
        combineFileName = strcat(selpathOut,'\Overall_Seg.mat');
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
        combineFileName = strcat(selpathOut,'\Overlay_Seg.fig');
        saveas(curFig,combineFileName);  
        combineFileName = strcat(selpathOut,'\Overlay_Seg.tif');
        saveas(curFig,combineFileName);
        close(curFig)
        
        
        
        
        
        
        
        chans = TrialInfo(:,1);
        uChans = unique(chans);
        combineFileName = strcat(selpathOut,'\channels.mat');
        save(combineFileName,'uChans'); 
        
        
        
        %% collect  baseline intensity values for each ROI
        baselineROItemp = cell(numel(regions),numel(TrialKeys));
        for i = 1:numel(TrialKeys) % load each trial dataset
            trialDat = load(TrialKeys{i}, 'preStim');
            curVol = trialDat.preStim;

            % for each ROI extract values from baseline
            for r = 1:numel(regions)
                curPixels = regions(r).PixelList;
                vals = zeros(size(curPixels,1),1);
                for m = 1:size(curPixels,1)
                    vals(m) = curVol(curPixels(m,2),curPixels(m,1),curPixels(m,3));
                end
                baselineROItemp{r,i} = vals;
            end
        end
        % condense all the values to individual 1-D arrays
        baselineROIvals = cell(numel(regions),1);
        for r = 1:numel(regions)
            ta = cell2mat(baselineROItemp(r,:));
            baselineROIvals{r} = reshape(ta,1,[]);
        end
        clear baselineROItemp
        
        

        %% Analyze the hits, misses, and currents for each channel
        for i = 1:numel(uChans)
            disp(strcat('  Processing channel: ',num2str(uChans(i))));
            
            curChTrialInfo = TrialInfo(TrialInfo(:,1)==uChans(i),:);
            uFreq = unique(curChTrialInfo(:,6));% trial frequency
            uPulse = unique(curChTrialInfo(:,7)); % trial pulse duration            
            for fi = 1:numel(uFreq)
                for pi = 1:numel(uPulse)

                    % Analyze each block
                    currents = unique(curChTrialInfo(:,3));

                    % Sweep through all currents in wide range of threshold
                    SingleCurpathOut = strcat(selpathOut,'\CurrentAnalysis');
                    mkdir(SingleCurpathOut); %Create output folder

                    rasterAllHit = false(numel(regions),numel(currents));
                    rasterAllMiss = false(numel(regions),numel(currents));

                    hitCounts = zeros(numel(currents),1);
                    MissCounts = zeros(numel(currents),1);
                    HitOnlyCounts = zeros(numel(currents),1);

                    currentTracesH = cell(numel(regions),numel(currents));
                    currentTracesM = cell(numel(regions),numel(currents));
                    hitTrials = zeros(numel(currents),1);
                    missTrials = zeros(numel(currents),1);

                    fprintf(strcat('    Processing currents for FREQ',num2str(uFreq(fi)),'PULSE',num2str(uPulse(pi))));
                    strCR = -1;
                    for curInd = 1:numel(currents)

                        % Display processing status
                        pVal = round(curInd/numel(currents)*100);     
                        percentageOut = [num2str(pVal) '%%'];
                        strOut = [percentageOut repmat(' ',1,10-length(percentageOut)-1)];
                        fprintf([strCR strOut]);
                        strCR = repmat('\b',1,length(strOut)-1);

                        k = currents(curInd); % select current value
                        curMask = zeros(size(TrialInfo,1),1);
                        curMask(TrialInfo(:,1) == uChans(i))=1; % only select trials for correct channel and current
                        curMask(TrialInfo(:,3) ~= k)=0; % stimulation current
                        curMask(TrialInfo(:,6) ~= uFreq(fi))=0; % stimulation frequency
                        curMask(TrialInfo(:,7) ~= uPulse(pi))=0; % stimulation pulse duration


                        if(sum(curMask)>0) % Only run on currents with data

                            curMaskHit = curMask;
                            curMaskHit(TrialInfo(:,4) == 0)=0; % Hit
                            fileNameAppHit = strcat('Ch',num2str(uChans(i)),'-CUR',num2str(k),'FREQ',num2str(uFreq(fi)),'PULSE',num2str(uPulse(pi)),'HIT');
                            curMaskMiss = curMask;
                            curMaskMiss(TrialInfo(:,4) == 1)=0; % Miss
                            fileNameAppMiss = strcat('Ch',num2str(uChans(i)),'-CUR',num2str(k),'FREQ',num2str(uFreq(fi)),'PULSE',num2str(uPulse(pi)),'MISS');

                            
                            % Process HITs
                            if(sum(curMaskHit)>0)
                                [rasterHit,curHROI,roiTracesHit,rasterHitInhibit,~,rasterPreStimH,~,postRasterH,postRasterInhibitH,~,~] = SecondaryProcess(TrialKeys,curMaskHit,fileNameAppHit,size(overallBaselineMean),standardDevCh{i},meanCh{i},selpathOut,regions,stdMod(i),baselineROIvals);
                            
                                rasterAllHit(rasterHit==1,curInd)=1;
                                currentTracesH(:,curInd) = roiTracesHit;
                                hitTrials(curInd) = sum(curMaskHit);
                                hitCounts(curInd) = sum(rasterHit);
                            else
                                rasterHit = false(numel(regions),1);
                                curHROI = false(size(overallBaselineMean));
                                roiTracesHit = false(numel(regions),1);
                                rasterHitInhibit = false(numel(regions),1);
                                rasterPreStimH = false(numel(regions),1);
                                postRasterH = false(numel(regions),1);
                                postRasterInhibitH = false(numel(regions),1);
                                
                                rasterAllHit(rasterHit==1,curInd)=1;
                                hitTrials(curInd) = sum(curMaskHit);
                                hitCounts(curInd) = sum(rasterHit);
                            end
                            
                            % Process MISSs
                            if(sum(curMaskMiss)>0)
                                [rasterMiss,curMROI,roiTracesMiss,rasterMissInhibit,~,rasterPreStimM,~,postRasterM,postRasterInhibitM,~,~] = SecondaryProcess(TrialKeys,curMaskMiss,fileNameAppMiss,size(overallBaselineMean),standardDevCh{i},meanCh{i},selpathOut,regions,stdMod(i),baselineROIvals);
                            
                                rasterAllMiss(rasterMiss==1,curInd)=1;
                                currentTracesM(:,curInd) = roiTracesMiss;
                                missTrials(curInd) = sum(curMaskMiss);
                                MissCounts(curInd) = sum(rasterMiss);
                            else
                                rasterMiss = false(numel(regions),1);
                                curMROI = false(size(overallBaselineMean));
                                roiTracesMiss = false(numel(regions),1);
                                rasterMissInhibit = false(numel(regions),1);
                                rasterPreStimM = false(numel(regions),1);
                                postRasterM = false(numel(regions),1);
                                postRasterInhibitM = false(numel(regions),1);
                                
                                rasterAllMiss(rasterMiss==1,curInd)=1;
                                missTrials(curInd) = sum(curMaskMiss);
                                MissCounts(curInd) = sum(rasterMiss);
                            end
                            
                            
                            
                            
                            
                            
                            
                            
                             % If both HITs and MISSs were present process
                             % comparision figures
                             if(sum(curMaskHit)>0 && sum(curMaskMiss)>0 )
                                hO = hitCounts;
                                hO(rasterMiss==1)=0;
                                HitOnlyCounts(curInd) = sum(hO);
                            

                            
                            
                            
%                             if(sum(curMaskHit)>0 && sum(curMaskMiss)>0 ) % only analyze where hits and missses are present for each current
%                                 [rasterHit,curHROI,roiTracesHit,rasterHitInhibit,~,rasterPreStimH,validROIPreStimH,postRasterH,postRasterInhibitH,validROIPostH,validROIPostInhibitH] = SecondaryProcess(TrialKeys,curMaskHit,fileNameAppHit,size(overallBaselineMean),standardDevCh{i},meanCh{i},selpathOut,regions,stdMod(i),baselineROIvals);
%                                 [rasterMiss,curMROI,roiTracesMiss,rasterMissInhibit,~,rasterPreStimM,validROIPreStimM,postRasterM,postRasterInhibitM,validROIPostM,validROIPostInhibitM] = SecondaryProcess(TrialKeys,curMaskMiss,fileNameAppMiss,size(overallBaselineMean),standardDevCh{i},meanCh{i},selpathOut,regions,stdMod(i),baselineROIvals);

%                                 rasterAllHit(rasterHit==1,curInd)=1;
%                                 rasterAllMiss(rasterMiss==1,curInd)=1;
% 
%                                 currentTracesH(:,curInd) = roiTracesHit;
%                                 currentTracesM(:,curInd) = roiTracesMiss;
%                                 hitTrials(curInd) = sum(curMaskHit);
%                                 missTrials(curInd) = sum(curMaskMiss);
% 
%                                 hitCounts(curInd) = sum(rasterHit);
%                                 MissCounts(curInd) = sum(rasterMiss);
%                                 hO = hitCounts;
%                                 hO(rasterMiss==1)=0;
%                                 HitOnlyCounts(curInd) = sum(hO);

                                % plot overlap of ROIs as lised in raster
                                fileNameApp = strcat(SingleCurpathOut,'\HITvMISS-OVERLAY-ROI-Ch',num2str(uChans(i)),'-CUR',num2str(k),'FREQ',num2str(uFreq(fi)),'PULSE',num2str(uPulse(pi)),'-roi');
                                hitVol = false(size(overallBaselineMean));
                                missVol = false(size(overallBaselineMean));
                                for k2 = 1:numel(regions)
                                    % generate current raster region
                                    curPixels = regions(k2).PixelList;
                                    if(rasterHit(k2)==1)
                                        for m = 1:size(curPixels,1)
                                            hitVol(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
                                        end 
                                    end
                                    if(rasterMiss(k2)==1)
                                        for m = 1:size(curPixels,1)
                                            missVol(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
                                        end   
                                    end
                                end
                                zProH = max(hitVol,[],3,'omitnan');
                                zProM = max(missVol,[],3,'omitnan');
                                curFig = figure();
                                image = zeros(size(curHROI,1),size(curHROI,2)); % Dark background
                                imshow(image)
                                zProHO = zProH;
                                zProHO(zProM==1)=0;
                                zProMO = zProM;
                                zProMO(zProH==1)=0;
                                zProB = zProH;
                                zProB(zProHO==1)=0;
                                A = imoverlay(image,max(zProMO,[],3),'red');
                                B = imoverlay(A,max(zProHO,[],3),'green');
                                C = imoverlay(B,max(zProB,[],3),'blue');
                                imshow(C)
                                titleStr1 = strcat('Current=',num2str(k),'\muA; Hit n=',num2str(sum(curMaskHit)),'; Miss n=',num2str(sum(curMaskMiss)),';  Events n=',num2str(sum(curMaskHit)+sum(curMaskMiss)));
                                title(titleStr1)
                                combineFileName = strcat(fileNameApp,'.fig');
                                saveas(curFig,combineFileName);  
                                combineFileName = strcat(fileNameApp,'.tif');
                                saveas(curFig,combineFileName);
                                close(curFig)
                                
                                
                                
                                
                                % Inhibited neural response
                                rasterAllHit(rasterHitInhibit==1,curInd)=1;
                                rasterAllMiss(rasterMissInhibit==1,curInd)=1;

                                currentTracesH(:,curInd) = roiTracesHit;
                                currentTracesM(:,curInd) = roiTracesMiss;
                                hitTrials(curInd) = sum(curMaskHit);
                                missTrials(curInd) = sum(curMaskMiss);

                                hitCounts(curInd) = sum(rasterHitInhibit);
                                MissCounts(curInd) = sum(rasterMissInhibit);
                                hO = hitCounts;
                                hO(rasterMissInhibit==1)=0;
                                HitOnlyCounts(curInd) = sum(hO);

                                % plot overlap of ROIs as lised in raster
                                fileNameApp = strcat(SingleCurpathOut,'\HITvMISS-OVERLAY-INHIBIT-ROI-Ch',num2str(uChans(i)),'-CUR',num2str(k),'FREQ',num2str(uFreq(fi)),'PULSE',num2str(uPulse(pi)),'-roi');
                                hitVol = false(size(overallBaselineMean));
                                missVol = false(size(overallBaselineMean));
                                for k2 = 1:numel(regions)
                                    % generate current raster region
                                    curPixels = regions(k2).PixelList;
                                    if(rasterHitInhibit(k2)==1)
                                        for m = 1:size(curPixels,1)
                                            hitVol(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
                                        end 
                                    end
                                    if(rasterMissInhibit(k2)==1)
                                        for m = 1:size(curPixels,1)
                                            missVol(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
                                        end   
                                    end
                                end
                                zProH = max(hitVol,[],3,'omitnan');
                                zProM = max(missVol,[],3,'omitnan');
                                curFig = figure();
                                image = zeros(size(curHROI,1),size(curHROI,2)); % Dark background
                                imshow(image)
                                zProHO = zProH;
                                zProHO(zProM==1)=0;
                                zProMO = zProM;
                                zProMO(zProH==1)=0;
                                zProB = zProH;
                                zProB(zProHO==1)=0;
                                A = imoverlay(image,max(zProMO,[],3),'red');
                                B = imoverlay(A,max(zProHO,[],3),'green');
                                C = imoverlay(B,max(zProB,[],3),'blue');
                                imshow(C)
                                titleStr1 = strcat('Current=',num2str(k),'\muA; Hit n=',num2str(sum(curMaskHit)),'; Miss n=',num2str(sum(curMaskMiss)),';  Events n=',num2str(sum(curMaskHit)+sum(curMaskMiss)));
                                title(titleStr1)
                                combineFileName = strcat(fileNameApp,'.fig');
                                saveas(curFig,combineFileName);  
                                combineFileName = strcat(fileNameApp,'.tif');
                                saveas(curFig,combineFileName);
                                close(curFig)
                                
                                
                                
                                
                                
                                
                                
                                
                                % Post stimulation
                                % plot overlap of ROIs as lised in raster
                                fileNameApp = strcat(SingleCurpathOut,'\PostStim-HITvMISS-OVERLAY-ROI-Ch',num2str(uChans(i)),'-CUR',num2str(k),'FREQ',num2str(uFreq(fi)),'PULSE',num2str(uPulse(pi)),'-roi');
                                hitVol = false(size(overallBaselineMean));
                                missVol = false(size(overallBaselineMean));
                                for k2 = 1:numel(regions)
                                    % generate current raster region
                                    curPixels = regions(k2).PixelList;
                                    if(postRasterH(k2)==1)
                                        for m = 1:size(curPixels,1)
                                            hitVol(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
                                        end 
                                    end
                                    if(postRasterM(k2)==1)
                                        for m = 1:size(curPixels,1)
                                            missVol(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
                                        end   
                                    end
                                end
                                zProH = max(hitVol,[],3,'omitnan');
                                zProM = max(missVol,[],3,'omitnan');
                                curFig = figure();
                                image = zeros(size(curHROI,1),size(curHROI,2)); % Dark background
                                imshow(image)
                                zProHO = zProH;
                                zProHO(zProM==1)=0;
                                zProMO = zProM;
                                zProMO(zProH==1)=0;
                                zProB = zProH;
                                zProB(zProHO==1)=0;
                                A = imoverlay(image,max(zProMO,[],3),'red');
                                B = imoverlay(A,max(zProHO,[],3),'green');
                                C = imoverlay(B,max(zProB,[],3),'blue');
                                imshow(C)
                                titleStr1 = strcat('Current=',num2str(k),'\muA; Hit n=',num2str(sum(curMaskHit)),'; Miss n=',num2str(sum(curMaskMiss)),';  Events n=',num2str(sum(curMaskHit)+sum(curMaskMiss)));
                                title(titleStr1)
                                combineFileName = strcat(fileNameApp,'.fig');
                                saveas(curFig,combineFileName);  
                                combineFileName = strcat(fileNameApp,'.tif');
                                saveas(curFig,combineFileName);
                                close(curFig)
                                
                                
                                
                                
                                % Inhibited neural response post
                                % stimulation

                                % plot overlap of ROIs as lised in raster
                                fileNameApp = strcat(SingleCurpathOut,'\PostStim-HITvMISS-OVERLAY-INHIBIT-ROI-Ch',num2str(uChans(i)),'-CUR',num2str(k),'FREQ',num2str(uFreq(fi)),'PULSE',num2str(uPulse(pi)),'-roi');
                                hitVol = false(size(overallBaselineMean));
                                missVol = false(size(overallBaselineMean));
                                for k2 = 1:numel(regions)
                                    % generate current raster region
                                    curPixels = regions(k2).PixelList;
                                    if(postRasterInhibitH(k2)==1)
                                        for m = 1:size(curPixels,1)
                                            hitVol(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
                                        end 
                                    end
                                    if(postRasterInhibitM(k2)==1)
                                        for m = 1:size(curPixels,1)
                                            missVol(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
                                        end   
                                    end
                                end
                                zProH = max(hitVol,[],3,'omitnan');
                                zProM = max(missVol,[],3,'omitnan');
                                curFig = figure();
                                image = zeros(size(curHROI,1),size(curHROI,2)); % Dark background
                                imshow(image)
                                zProHO = zProH;
                                zProHO(zProM==1)=0;
                                zProMO = zProM;
                                zProMO(zProH==1)=0;
                                zProB = zProH;
                                zProB(zProHO==1)=0;
                                A = imoverlay(image,max(zProMO,[],3),'red');
                                B = imoverlay(A,max(zProHO,[],3),'green');
                                C = imoverlay(B,max(zProB,[],3),'blue');
                                imshow(C)
                                titleStr1 = strcat('Current=',num2str(k),'\muA; Hit n=',num2str(sum(curMaskHit)),'; Miss n=',num2str(sum(curMaskMiss)),';  Events n=',num2str(sum(curMaskHit)+sum(curMaskMiss)));
                                title(titleStr1)
                                combineFileName = strcat(fileNameApp,'.fig');
                                saveas(curFig,combineFileName);  
                                combineFileName = strcat(fileNameApp,'.tif');
                                saveas(curFig,combineFileName);
                                close(curFig)
                                
                                
                                
                                
                                
                                
                                
                                % Perform comparision between mid stim and
                                % post stim (HITS)
                                
                                % Excitation
                                mid2Post = rasterHit;
                                mid2Post(postRasterH==0)=0;
                                pre2Mid = rasterHit;
                                pre2Mid(rasterPreStimH==0)=0;
                                
                                
                                
                                %plot total values, delta of values, and
                                %consistency of regions
                                curFig = figure();
                                X = categorical({'Pre','Mid','Post','Pre2Mid','Mid2Post'});
                                Y = [sum(rasterPreStimH),sum(rasterHit),sum(postRasterH),sum(pre2Mid),sum(mid2Post)];
                                bar(X,Y)
                                titleStr1 = strcat('Activation Consistency HIT Current=',num2str(k),'\muA');
                                title(titleStr1)
                                fileNameApp = strcat(SingleCurpathOut,'\ActivationConsistencyHIT-Ch',num2str(uChans(i)),'-CUR',num2str(k),'FREQ',num2str(uFreq(fi)),'PULSE',num2str(uPulse(pi)),'-roi');
                                combineFileName = strcat(fileNameApp,'.fig');
                                saveas(curFig,combineFileName);  
                                combineFileName = strcat(fileNameApp,'.tif');
                                saveas(curFig,combineFileName);
                                close(curFig)
                                
                                
                                
                                % Inhibition (can only look at consistency
                                % for mid to post, there isn't really a
                                % good way to quantify baseline inhibition
                                % by definition)
                                mid2PostInhib = rasterHitInhibit;
                                mid2PostInhib(postRasterInhibitH==0)=0;
                                
                                curFig = figure();
                                X = categorical({'Mid','Post','Mid2Post'});
                                Y = [sum(rasterHitInhibit), sum(postRasterInhibitH), sum(mid2PostInhib)];
                                bar(X,Y)
                                titleStr1 = strcat('Inhibition Consistency HIT Current=',num2str(k),'\muA');
                                title(titleStr1)
                                fileNameApp = strcat(SingleCurpathOut,'\InhibitionConsistencyHIT-Ch',num2str(uChans(i)),'-CUR',num2str(k),'FREQ',num2str(uFreq(fi)),'PULSE',num2str(uPulse(pi)),'-roi');
                                combineFileName = strcat(fileNameApp,'.fig');
                                saveas(curFig,combineFileName);  
                                combineFileName = strcat(fileNameApp,'.tif');
                                saveas(curFig,combineFileName);
                                close(curFig)
                                
                                
                                
                                
                                
                                % Perform comparision between mid stim and
                                % post stim (MISS)
                                
                                % Excitation
                                mid2Post = rasterMiss;
                                mid2Post(postRasterM==0)=0;
                                pre2Mid = rasterMiss;
                                pre2Mid(rasterPreStimM==0)=0;
                                curFig = figure();
                                X = categorical({'Pre','Mid','Post','Pre2Mid','Mid2Post'});
                                Y = [sum(rasterPreStimM),sum(rasterMiss),sum(postRasterM),sum(pre2Mid),sum(mid2Post)];
                                bar(X,Y)
                                titleStr1 = strcat('Activation Consistency MISS Current=',num2str(k),'\muA');
                                title(titleStr1)
                                fileNameApp = strcat(SingleCurpathOut,'\ActivationConsistencyMISS-Ch',num2str(uChans(i)),'-CUR',num2str(k),'FREQ',num2str(uFreq(fi)),'PULSE',num2str(uPulse(pi)),'-roi');
                                combineFileName = strcat(fileNameApp,'.fig');
                                saveas(curFig,combineFileName);  
                                combineFileName = strcat(fileNameApp,'.tif');
                                saveas(curFig,combineFileName);
                                close(curFig)
                                                               
                                
                                % Inhibition
                                mid2PostInhib = rasterMissInhibit;
                                mid2PostInhib(postRasterInhibitM==0)=0;
                                curFig = figure();
                                X = categorical({'Mid','Post','Mid2Post'});
                                Y = [sum(rasterMissInhibit),sum(postRasterInhibitM),sum(mid2PostInhib)];
                                bar(X,Y)
                                titleStr1 = strcat('Inhibition Consistency MISS Current=',num2str(k),'\muA');
                                title(titleStr1)
                                fileNameApp = strcat(SingleCurpathOut,'\InhibitionConsistencyMISS-Ch',num2str(uChans(i)),'-CUR',num2str(k),'FREQ',num2str(uFreq(fi)),'PULSE',num2str(uPulse(pi)),'-roi');
                                combineFileName = strcat(fileNameApp,'.fig');
                                saveas(curFig,combineFileName);  
                                combineFileName = strcat(fileNameApp,'.tif');
                                saveas(curFig,combineFileName);
                                close(curFig)


                                % Detected Regions
                                compairROI = curHROI;
                                compairROI(curMROI==1)=0;
                                fileNameApp = strcat(SingleCurpathOut,'\HITvMISS_Ch',num2str(uChans(i)),'_CUR',num2str(k),'FREQ',num2str(uFreq(fi)),'PULSE',num2str(uPulse(pi)),'_roi');
                                curFig = figure();
                                zProjected = max(compairROI,[],3,'omitnan');
                                imshow(zProjected)
                                titleStr1 = strcat('Current=',num2str(k),'\muA; Hit n=',num2str(sum(curMaskHit)),'; Miss n=',num2str(sum(curMaskMiss)),';  Events n=',num2str(sum(curMaskHit)+sum(curMaskMiss)));
                                title(titleStr1)
                                combineFileName = strcat(fileNameApp,'.fig');
                                saveas(curFig,combineFileName);  
                                combineFileName = strcat(fileNameApp,'.tif');
                                saveas(curFig,combineFileName);
                                close(curFig)
                             end

                            combineFileName = strcat(SingleCurpathOut,'\HITvMISS_METRICS_Ch',num2str(uChans(i)),'_CUR',num2str(k),'FREQ',num2str(uFreq(fi)),'PULSE',num2str(uPulse(pi)));
                            activeCurr = k;
                            activeChan = uChans(i);
                            save(combineFileName,'activeChan','activeCurr','curHROI','curMROI','elecPos'); 

%                             end
                        end
                    end
                    fprintf('\n')

% %                     hitOnly = rasterAllHit;
% %                     hitOnly(rasterAllMiss==1)=0;
% %                     hitSum = sum(hitOnly,2);
% % 
% %                     % Determine valid ROIs to plot (i.e. only those that were
% %                     % active at some current)
% %                     validInds = find(hitSum);
% %                     if(validInds>0)
% %         
% % 
% %                         % Plot traces from the current channel over the currents
% %                         % analyzed highlighting trials that were hits and misses
% %                         % Save z-projection of mean mid-stimulus imageBlock
% %                         curFig = figure();
% %                         vertStep = 20;
% % 
% % 
% %                         validInds2 = find(sum(rasterAllHit,2));
% %                         vertRange = vertStep*numel(validInds2)+10;
% %                         totLength = sum(hitTrials) + sum(missTrials);
% %                         plot([1,totLength],[1,2],'Color',[0, 0, 0, 0.01]) % plot initial trace to establish figure
% %                         hold on
% % 
% %                         priorTrials = 0;
% %                         for c = 1:numel(currents)% highlight individual current regions
% %                             HT = hitTrials(c);
% %                             MT = missTrials(c);
% %                             hs = priorTrials + 0.5;
% %                             he = priorTrials + HT + 0.5;
% %                             me = priorTrials + HT + MT + 0.5;
% % 
% %                             if(HT~=0)
% %                                 patch([hs he he hs], [0 0 vertRange vertRange], 'green', 'FaceAlpha', 0.5, 'EdgeColor', 'none')
% %                             end
% %                             if(MT~=0)
% %                                 patch([he me me he], [0 0 vertRange vertRange], 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'none')
% %                             end
% %                             xline(me,'--'); % Plot vertical line deliniating currents
% % 
% %                             priorTrials = priorTrials + HT + MT;
% %                         end
% %                         hold on
% % 
% %                         priorTrials = 0;
% %                         for c = 1:numel(currents)
% %                             HT = hitTrials(c);
% %                             MT = missTrials(c);
% %                             ts = priorTrials + 1;
% %                             te = priorTrials + HT + MT;
% %                             for r = 1:numel(regions) % plot traces
% %                                 if(rasterAllHit(r,c)==1)
% %                                     vs = find(validInds2==r);
% % 
% %                                     curHTrace = currentTracesH{r,c};
% %                                     curMTrace = currentTracesM{r,c}; 
% %                                     curTrace = [curHTrace;curMTrace];
% % 
% %                                     plot(ts:te,curTrace+vs*vertStep,'Color',[0, 0, 0, 1])
% %                                 end
% %                             end
% %                             priorTrials = priorTrials + HT + MT;
% %                         end
% %                         combineFileName = strcat(selpathOut,'\ROItracesV2_ch',num2str(uChans(i)),'.fig');
% %                         saveas(curFig,combineFileName);  
% %                         combineFileName = strcat(selpathOut,'\ROItracesV2_ch',num2str(uChans(i)),'.tif');
% %                         saveas(curFig,combineFileName);
% %                         close(curFig)
% %                     end
            
                end
            end
        end
        
        
        
        %% perform interblock analsysis of trails
        
        % generate trial block alignment variable
        trialTot = 10;
        TrialBlockKey = zeros(trialTot,1);
        trialCount = 1;
        for curBlock = 1:numel(keyListing) % Process each channel/staircase/Block dataset

           % Load timing key data
            dataCellArr = load(strcat(keypathO,'\',keyListing(keyInds(curBlock)).name));
            dataCellArr = dataCellArr.dataCellArr;
            eCell = zeros(numel(dataCellArr),1);
            for k = 1:numel(dataCellArr)
                if(isempty(dataCellArr{k}))
                    eCell(k)=1;
                end
            end
            dataCellArr(eCell==1)=[];

            for trialNum = 1:numel(dataCellArr) % process each trial

                % Valid trial, process data. Trials with error flags
                % will be ignored/discarded
                trial_i = dataCellArr{trialNum};
                validTrial = trial_i{7};

                if(validTrial)
                    TrialBlockKey(trialCount) = curBlock; % assign block ID
                    trialCount = trialCount+1; % increment trial count
                end
            end
        end
        uBlocks = unique(TrialBlockKey);

        

        %% seperate analysis by the hits, misses, and currents for each
        % channel and compare average volume flourescence
        
        hitMeans = zeros(numel(uBlocks),1);
        deltaMeans = zeros(numel(uBlocks),1);
        hitROInum = zeros(numel(uBlocks),1);
        comCnt = 0;
        imgSz = size(overallBaselineMean);
        
        for i = 1:numel(uChans)

            curChTrialInfo = TrialInfo(TrialInfo(:,1)==uChans(i),:);
            uFreq = unique(curChTrialInfo(:,6));% trial frequency
            uPulse = unique(curChTrialInfo(:,7)); % trial pulse duration    
            currents = unique(curChTrialInfo(:,3)); % currents

            for fi = 1:numel(uFreq)
                for pi = 1:numel(uPulse)
                    for curInd = 1:numel(currents) % Sweep through all currents
                        curr = currents(curInd);
                        
                        curMask = zeros(size(TrialInfo,1),1);
                        curMask(TrialInfo(:,1) == uChans(i))=1; % only select trials for correct channel and current
                        curMask(TrialInfo(:,3) ~= curr)=0; % stimulation current
                        curMask(TrialInfo(:,6) ~= uFreq(fi))=0; % stimulation frequency
                        curMask(TrialInfo(:,7) ~= uPulse(pi))=0; % stimulation pulse duration
                        
                        if(sum(curMask)>0) % Only run on currents with data
                            blockMean = zeros(numel(uBlocks),size(overallBaselineMean,1),size(overallBaselineMean,2),size(overallBaselineMean,3));
                            blockDelta = zeros(numel(uBlocks),size(overallBaselineMean,1),size(overallBaselineMean,2),size(overallBaselineMean,3));
                            invFlag = 0;
                            for b = 1:numel(uBlocks)

                                curMaskHit = curMask;
                                curMaskHit(TrialBlockKey ~= uBlocks(b))=0; % block check
                                curMaskHit(TrialInfo(:,4) == 0)=0; % Hit

                                if(sum(curMaskHit)>0)
                                    % Construct variables pulling data from across all included HIT trials
                                    preStim = cell(sum(curMaskHit),1);
                                    midStim = cell(sum(curMaskHit),1);
                                    trialInds = find(curMaskHit);
                                    for tInds = 1:numel(trialInds) % load each trial dataset
                                        trialDat = load(TrialKeys{trialInds(tInds)}, 'preStim','midStim');
                                        preStim{tInds} = trialDat.preStim;
                                        midStim{tInds} = trialDat.midStim;
                                    end    
                                    % Average flourescence response of all trials calculate regions with intesity delta over baseline
                                    meanStimDelta = zeros(imgSz(1),imgSz(2),imgSz(3),numel(trialInds));
                                    meanMidStimOrig = zeros(imgSz(1),imgSz(2),imgSz(3),numel(trialInds));
                                    for k = 1:numel(trialInds)
                                        meanStimDelta(:,:,:,k) = midStim{k}-preStim{k};
                                        meanMidStimOrig(:,:,:,k) = midStim{k};
                                    end
                                    blockMean(b,:,:,:) = squeeze(mean(meanStimDelta,4));
                                    blockDelta(b,:,:,:) = squeeze(mean(meanMidStimOrig,4));
                                else
                                    invFlag = 1;
                                end
                            end
                            
                            if(invFlag==0)
                                comCnt = comCnt+1;

                                
                                
                                for b = 1:numel(uBlocks)
                                    curMaskHit = curMask;
                                    curMaskHit(TrialBlockKey ~= uBlocks(b))=0; % block check
                                    curMaskHit(TrialInfo(:,4) == 0)=0; % Hit
                                    [rasterHit,~,~,~,~,~,~,~,~,~,~] = SecondaryProcess(TrialKeys,curMaskHit,"blockanalysis",size(overallBaselineMean),standardDevCh{i},meanCh{i},selpathOut,regions,stdMod(i),baselineROIvals);

                                    
                                    hitROInum(b,comCnt) = sum(rasterHit); % The number of ROIs active/detected
                                    hitMeans(b,comCnt) = mean(blockMean(b,:,:,:) ,'all','omitnan');
                                    deltaMeans(b,comCnt) = mean(blockDelta(b,:,:,:) ./ blockDelta(1,:,:,:),'all','omitnan');
                                end
                            end
                        end                        
                    end
                end
            end
        end
        
        
       %% % plot block to block variation for matched trials
        meanVals = [];
        deltaVals = [];
        ROIvals = [];
        g = [];
        
        for b = 1:numel(uBlocks)
            ROIvals = [ROIvals; hitROInum(b,:)'];
            meanVals = [meanVals; hitMeans(b,:)'];
            deltaVals = [deltaVals; deltaMeans(b,:)'];
            g = [g; b*ones(size(deltaMeans,2),1)];

        end
        
        % Mean Mid-stimulation flourescence
        curFig = figure();
        boxplot(meanVals,g);
        ylabel('Trial Intensity')
        xlabel('Imaging block')
        title('Midstim HIT Intensity')
        
        combineFileName = strcat(selpathOut,'BlockComparisionMidStim.fig');
        saveas(curFig,combineFileName);  
        combineFileName = strcat(selpathOut,'BlockComparisionMidStim.tif');
        saveas(curFig,combineFileName);  
        close(curFig)
        
        % Normalized Delta of flouresecence
        curFig = figure();
        boxplot(deltaVals,g);
        ylabel('Normalized Trial Intensity')
        xlabel('Imaging block')
        title('Delta HIT Intensity')
        
        combineFileName = strcat(selpathOut,'BlockComparisionDelta.fig');
        saveas(curFig,combineFileName);  
        combineFileName = strcat(selpathOut,'BlockComparisionDelta.tif');
        saveas(curFig,combineFileName);  
        close(curFig)
        
        % Number of detected ROIs
        curFig = figure();
        boxplot(ROIvals,g);
        ylabel('Number of detected ROIs')
        xlabel('Imaging block')
        title('Detected ROI variability')
        
        combineFileName = strcat(selpathOut,'BlockComparisionROInum.fig');
        saveas(curFig,combineFileName);  
        combineFileName = strcat(selpathOut,'BlockComparisionROInum.tif');
        saveas(curFig,combineFileName);  
        close(curFig)
        
        
        
        combineFileName = strcat(selpathOut,'BlockComparisionData.mat');
        save(combineFileName,'hitMeans','deltaMeans','ROIvals','uBlocks')
        
        % Number of HIT vs MISS per block
        
        
        
    end
    fprintf('Finished Processing\n\n');
end




















function [midRaster,validROI,roiTraces,midRasterInhibit,validROIInhibit,rasterPreStim,validROIPreStim,postRaster,postRasterInhibit,validROIPost,validROIPostInhibit] = SecondaryProcess(TrialKeys,trialMask,fileNameApp,imgSz,standardDev,overallBaselineMean,selpathOut,regions,stdMod,baselineROIvals)
% SecondaryProcess
%   Function to perform secondary processing of trials which are included
%   in the mask of trials.
%     stdMod=2; % variable determining detection threshold 
  
    curTraceDir = strcat(selpathOut,'\',fileNameApp);
    mkdir(curTraceDir)
    
    % Construct variables pulling data from across all included trials
    preStim = cell(sum(trialMask),1);
    midStim = cell(sum(trialMask),1);
    postStim = cell(sum(trialMask),1);

    trialInds = find(trialMask);
    midRaster = zeros(numel(regions),1);%zeros(numel(regions),numel(trialInds));
    rasterPreStim = zeros(numel(regions),1);
    midRasterInhibit = zeros(numel(regions),1);
    postRaster = zeros(numel(regions),1);
    postRasterInhibit = zeros(numel(regions),1);
    
    
    for i = 1:numel(trialInds) % load each trial dataset
        trialDat = load(TrialKeys{trialInds(i)}, 'preStim','midStim','postStim');
        preStim{i} = trialDat.preStim;
        midStim{i} = trialDat.midStim;
        postStim{i} = trialDat.postStim;
    end    
    
    

    

    %% Average flourescence response of all trials calculate regions with intesity delta over baseline
    meanPreStimOrig = zeros(imgSz(1),imgSz(2),imgSz(3),numel(trialInds));
    meanMidStimOrig = zeros(imgSz(1),imgSz(2),imgSz(3),numel(trialInds));
    meanPostStimOrig = zeros(imgSz(1),imgSz(2),imgSz(3),numel(trialInds));
    for k = 1:numel(trialInds)
        meanPreStimOrig(:,:,:,k) = preStim{k};
        meanMidStimOrig(:,:,:,k) = midStim{k};
        meanPostStimOrig(:,:,:,k) = postStim{k};
    end
    
    meanPreStim = mean(meanPreStimOrig,4,'omitnan'); % Baseline neural data
    meanMidStim = mean(meanMidStimOrig,4,'omitnan');
    meanPostStim = mean(meanPostStimOrig,4,'omitnan');
    
%     % plot individual images of pre and mid-stim volumes (debugging)
%     maxPreStim = max(meanPreStimOrig,[],'all','omitnan'); % Baseline neural data
%     maxMidStim = max(meanMidStimOrig,[],'all','omitnan');
%     for k = 1:numel(trialInds)
%         curFig = figure();
%         zProjected = max(preStim{k},[],3,'omitnan');
%         imshow(zProjected/maxPreStim)
%         title(strcat('PreStim-',num2str(k),'-AveInt:',num2str(mean(preStim{k},'all','omitnan'))))
%         combineFileName = strcat(curTraceDir,'\RawPreStim',num2str(k),'.tif');
%         saveas(curFig,combineFileName);
%         close(curFig)
% 
%         curFig = figure();
%         zProjected = max(midStim{k},[],3,'omitnan');
%         imshow(zProjected/maxMidStim)
%         title(strcat('MidStim-',num2str(k),'-AveInt:',num2str(mean(midStim{k},'all','omitnan'))))
%         combineFileName = strcat(curTraceDir,'\RawMidStim',num2str(k),'.tif');
%         saveas(curFig,combineFileName);
%         close(curFig)
%     end
    
    
    
    midDelta = meanMidStim - meanPreStim;
    postDelta = meanPostStim - meanPreStim;
    midActive = false(size(midDelta));
    midActive(midDelta>(stdMod*standardDev))=1;
    postActive = false(size(postDelta));
    postActive(postDelta>(stdMod*standardDev))=1;
    midInactive = false(size(midDelta));
    midInactive(midDelta<(-stdMod*standardDev))=1;
    
    % Save mean pre-stimulus imageBlock
    combineFileName = strcat(curTraceDir,'\PreStim.mat');
    save(combineFileName,'meanPreStim');  
    
    % Save z-projection of mean mid-stimulus imageBlock
    curFig = figure();
    zProjected = max(meanMidStim,[],3,'omitnan');
    imshow(zProjected/max(zProjected,[],'all'))
    combineFileName = strcat(curTraceDir,'\MidStim.fig');
    saveas(curFig,combineFileName);  
    combineFileName = strcat(curTraceDir,'\MidStim.tif');
    saveas(curFig,combineFileName);
    close(curFig)
    combineFileName = strcat(curTraceDir,'\MidStim.mat');
    save(combineFileName,'meanMidStim');  


    
    
   
   rasterPostInhibit = zeros(numel(regions),numel(trialInds));
   rasterExcite = zeros(numel(regions),numel(trialInds));
   rasterInhibit = zeros(numel(regions),numel(trialInds));
   rasterPre = zeros(numel(regions),numel(trialInds));
   rasterExcite2 = zeros(numel(regions),numel(trialInds));
   rasterInhibit2 = zeros(numel(regions),numel(trialInds));
   
	for i = 1:numel(regions)
        curPixels = regions(i).PixelList;
        
        % Quantify baseline activity for current ROI over all trials
        baselineVals = zeros(numel(trialInds),size(curPixels,1));
        postStimVals = zeros(numel(trialInds),size(curPixels,1));
        midStimVals = zeros(numel(trialInds),size(curPixels,1));
        for k = 1:numel(trialInds)
            preVol = preStim{k};
            midVol = midStim{k};
            postVol = postStim{k};
            for m = 1:size(curPixels,1)
                baselineVals(k,m) = preVol(curPixels(m,2),curPixels(m,1),curPixels(m,3));
                midStimVals(k,m) = midVol(curPixels(m,2),curPixels(m,1),curPixels(m,3));
                postStimVals(k,m) = postVol(curPixels(m,2),curPixels(m,1),curPixels(m,3));
            end  
        end
%         preStimFlat = reshape(baselineVals,1,[]);
        baselineFlat = baselineROIvals{i};
       
        % Check to see if mid stim values are significantly greater or
        % lower than baseline activity
        for k = 1:numel(trialInds)
            if(ttest2(midStimVals(k,:),baselineFlat))
                if(mean(midStimVals(k,:),'all')>mean(baselineFlat,'all'))
                    % Excitation
                    rasterExcite(i,k)=1;
                else
                    % Inhibition
                    rasterInhibit(i,k)=1;
                end
            end
            
            if(ttest2(postStimVals(k,:),baselineFlat))
                if(mean(postStimVals(k,:),'all')<mean(baselineFlat,'all'))
                    % Inhibition
                    rasterPostInhibit(i,k)=1;
                end
            end
            
            
            if(ttest2(midStimVals(k,:),baselineVals(k,:)))
                if(mean(midStimVals(k,:),'all')>mean(baselineVals(k,:),'all'))
                    % Excitation
                    rasterExcite2(i,k)=1;
                else
                    % Inhibition
                    rasterInhibit2(i,k)=1;
                end
            end
            
            % Check if the pre-stimulation values are significantly greater
            % than average baseline
            if(ttest2(baselineVals(k,:),baselineFlat))
                if(mean(baselineVals(k,:),'all')>mean(baselineFlat,'all'))
                    % Excitation
                    rasterPre(i,k)=1;
                end
            end
        end
	end
   
	midRasterInhibit(mean(rasterInhibit,2)>0.5)=1;
    postRasterInhibit(mean(rasterPostInhibit,2)>0.5)=1;
    rasterPreStim(mean(rasterPre,2)>0.5)=1;
    validROIInhibit = false(size(overallBaselineMean));
    validROIPreStim = false(size(overallBaselineMean));
    validROIPostInhibit = false(size(overallBaselineMean));
    
    % Activated pre-stimulation regions  
    for i = 1:numel(regions)
        if(rasterPreStim(i)==1)
            curPixels = regions(i).PixelList;
            for m = 1:size(curPixels,1)
                validROIPreStim(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
            end
        end
    end
    
    % Inhibited regions
    for i = 1:numel(regions)
        if(midRasterInhibit(i)==1)
            curPixels = regions(i).PixelList;
            for m = 1:size(curPixels,1)
                validROIInhibit(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
            end
        end
    end
    
    % Inhibited post stim regions
    for i = 1:numel(regions)
        if(postRasterInhibit(i)==1)
            curPixels = regions(i).PixelList;
            for m = 1:size(curPixels,1)
                validROIPostInhibit(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
            end
        end
    end
   
    %% Process all regions 
    % Check each region to see if 50% of voxels are detected as over intensity threshold
    validROIPost = false(size(overallBaselineMean));
    validROI = false(size(overallBaselineMean));
    possibleROI = false(size(overallBaselineMean));
    for i = 1:numel(regions)
        % Pull sum of activations data from overall dataset to 
        curPixels = regions(i).PixelList;
        curMidVals = zeros(size(curPixels,1),1);
        curPostVals = zeros(size(curPixels,1),1);
        curVol = false(size(overallBaselineMean));
        for m = 1:size(curPixels,1)
            curMidVals(m) = midActive(curPixels(m,2),curPixels(m,1),curPixels(m,3));
            curPostVals(m) = postActive(curPixels(m,2),curPixels(m,1),curPixels(m,3));
            curVol(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
        end   
        if(mean(curMidVals)>=0.5)
            midRaster(i,:)=1;
            validROI(curVol)=1;
        end
        if(mean(curPostVals)>=0.5)
            postRaster(i,:)=1;
            validROIPost(curVol)=1;
        end
        if(sum(curMidVals)>=1)
            possibleROI(curVol)=1;
        end
    end
        
    curFig = figure();
    zProjected = max(midActive,[],3,'omitnan');
    imshow(zProjected/max(zProjected,[],'all'))
    title(strcat('Mid Stimulation activation : ',num2str(numel(trialInds)),' Trials'))
    combineFileName = strcat(curTraceDir,'\MidActive.fig');
    saveas(curFig,combineFileName);  
    combineFileName = strcat(curTraceDir,'\MidActive.tif');
    saveas(curFig,combineFileName);
    close(curFig)
    
    
    curFig = figure();
    zProjected = max(midInactive,[],3,'omitnan');
    imshow(zProjected/max(zProjected,[],'all'))
    title(strcat('Mid Stimulation Inactivation : ',num2str(numel(trialInds)),' Trials'))
    combineFileName = strcat(curTraceDir,'\MidInactive.fig');
    saveas(curFig,combineFileName);  
    combineFileName = strcat(curTraceDir,'\MidInactive.tif');
    saveas(curFig,combineFileName);
    close(curFig)
    
    
    % Save mid-stimulation valid ROIs
    curFig = figure();
    zProjectedBW = max(validROI,[],3,'omitnan');
    imshow(zProjectedBW/max(zProjectedBW,[],'all'))
    title('Z projected ROI segmentation')
    combineFileName = strcat(curTraceDir,'\ROI.fig');
    saveas(curFig,combineFileName);  
    combineFileName = strcat(curTraceDir,'\ROI.tif');
    saveas(curFig,combineFileName);
    close(curFig)
    combineFileName = strcat(curTraceDir,'\ROI.mat');
    save(combineFileName,'validROI','midRaster','regions','midActive','validROIInhibit','midInactive','rasterPreStim','validROIPreStim','postRaster','postRasterInhibit','validROIPost','validROIPostInhibit'); 
    
    curFig = figure();
    zProjectedBW = max(validROIInhibit,[],3,'omitnan');
    imshow(zProjectedBW/max(zProjectedBW,[],'all'))
    title('Z projected Inactive ROI segmentation')
    combineFileName = strcat(curTraceDir,'\ROI-Inhibit.fig');
    saveas(curFig,combineFileName);  
    combineFileName = strcat(curTraceDir,'\ROI-Inhibit.tif');
    saveas(curFig,combineFileName);
    close(curFig)
    
    curFig = figure();
    zProjectedBW = max(validROIPost,[],3,'omitnan');
    imshow(zProjectedBW/max(zProjectedBW,[],'all'))
    title('Z projected Post Stim ROI segmentation')
    combineFileName = strcat(curTraceDir,'\ROI-Post.fig');
    saveas(curFig,combineFileName);  
    combineFileName = strcat(curTraceDir,'\ROI-Post.tif');
    saveas(curFig,combineFileName);
    close(curFig)
        
    curFig = figure();
    zProjectedBW = max(validROIPostInhibit,[],3,'omitnan');
    imshow(zProjectedBW/max(zProjectedBW,[],'all'))
    title('Z projected Post Stim Inactive ROI segmentation')
    combineFileName = strcat(curTraceDir,'\ROI-Post-Inhibit.fig');
    saveas(curFig,combineFileName);  
    combineFileName = strcat(curTraceDir,'\ROI-Post-Inhibit.tif');
    saveas(curFig,combineFileName);
    close(curFig)
    
    curFig = figure();
    zProjectedBW = max(validROIPreStim,[],3,'omitnan');
    imshow(zProjectedBW/max(zProjectedBW,[],'all'))
    title('Z projected PreStim ROI segmentation')
    combineFileName = strcat(curTraceDir,'\ROI-PreStimActive.fig');
    saveas(curFig,combineFileName);  
    combineFileName = strcat(curTraceDir,'\ROI-PreStimActive.tif');
    saveas(curFig,combineFileName);
    close(curFig)

    curFig = figure();
    zProjectedBW = max(possibleROI,[],3,'omitnan');
    imshow(zProjectedBW/max(zProjectedBW,[],'all'))
    title('Z projected possible ROI segmentation')
    combineFileName = strcat(curTraceDir,'\ROIPoss.fig');
    saveas(curFig,combineFileName);  
    combineFileName = strcat(curTraceDir,'\ROIPoss.tif');
    saveas(curFig,combineFileName);
    close(curFig)
    
    
    % Save delta volume
    curFig = figure();
    zProjected = max(midDelta,[],3,'omitnan');
    imshow(zProjected/max(zProjected,[],'all'))
    title(strcat('Delta of stimulated Volume : ',num2str(numel(trialInds)),' Trials'))
    combineFileName = strcat(curTraceDir,'\DeltaVol.fig');
    saveas(curFig,combineFileName);  
    combineFileName = strcat(curTraceDir,'\DeltaVol.tif');
    saveas(curFig,combineFileName);
    close(curFig)
    combineFileName = strcat(curTraceDir,'\DeltaVol.mat');
    save(combineFileName,'midDelta'); 
    
    
    % Save delta volume vs baseline STD
    sharedMax = 1.1*max(stdMod*standardDev,[],'all');
    curFig = figure();
    subplot(2,1,1)
    zProjected = max(midDelta,[],3,'omitnan');
    imshow(zProjected/sharedMax)
    title(strcat('Delta of stimulated Volume : ',num2str(numel(trialInds)),' Trials'))
    xlabel(num2str(mean(midDelta,'all','omitnan')))
    subplot(2,1,2)
    zProjected = max(stdMod*standardDev,[],3,'omitnan');
    imshow(zProjected/sharedMax)
    title('Standard Deviation')
    xlabel(num2str(mean(stdMod*standardDev,'all','omitnan')))
    combineFileName = strcat(curTraceDir,'\DeltaVolSTD.fig');
    saveas(curFig,combineFileName);  
    combineFileName = strcat(curTraceDir,'\DeltaVolSTD.tif');
    saveas(curFig,combineFileName);
    close(curFig)
   
    
    
    % Plot outlines of detected regions
    zProjectedBW = max(validROI,[],3,'omitnan');
    boundaries = bwboundaries(zProjectedBW);
    curFig = figure();
    title('Overlay of ')
    zProjected = max(midDelta,[],3,'omitnan');
    imshow(zProjected/max(zProjected,[],'all'))
    hold on
    for k=1:numel(boundaries)
       b = boundaries{k};
       plot(b(:,2),b(:,1),'r','LineWidth',1);
    end
    combineFileName = strcat(curTraceDir,'\Overlay_ROI.fig');
    saveas(curFig,combineFileName);  
    combineFileName = strcat(curTraceDir,'\Overlay_ROI.tif');
    saveas(curFig,combineFileName);
    close(curFig)
    
    
    
    % Plot z projections of volumes processed
    curFig = figure();
    stdMax = 1.2 * max((stdMod*standardDev),[],'all','omitnan');
    subplot(2,2,1)
    ylabel('Pre-Stimulation')
    zProjected = max(meanPreStim,[],3,'omitnan');
    imshow(zProjected/stdMax)
    xlabel(strcat('Ave Intensity: ',num2str(mean(meanPreStim,'all'))))
    subplot(2,2,2)
    ylabel('Mid-Stimulation')
    zProjected = max(meanMidStim,[],3,'omitnan');
    imshow(zProjected/stdMax)
    xlabel(strcat('Ave Intensity: ',num2str(mean(meanMidStim,'all'))))
    subplot(2,2,3)
    ylabel('Stimulation Delta')
    zProjected = max(midDelta,[],3,'omitnan');
    imshow(zProjected/stdMax)
    xlabel(strcat('Ave Intensity: ',num2str(mean(midDelta,'all'))))
    subplot(2,2,4)
    ylabel('Baseline STD')
    zProjected = max((stdMod*standardDev),[],3,'omitnan');
    imshow(zProjected/stdMax)
    xlabel(strcat('Ave Intensity: ',num2str(mean(stdMod*standardDev,'all'))))
    combineFileName = strcat(curTraceDir,'\VolumeProcessing.fig');
    saveas(curFig,combineFileName);  
    combineFileName = strcat(curTraceDir,'\VolumeProcessing.tif');
    saveas(curFig,combineFileName);
    close(curFig)
    
    

    
    
    
    % Plot traces for ROIs that were detected
    roiTraces = cell(numel(regions),1);
    for r = 1:numel(regions)
        curPixels = regions(r).PixelList;
        curVol = false(size(overallBaselineMean));
%         curThresh = 1;
        trialPreMean = zeros(numel(trialInds),1);
        trialPreSTD = zeros(numel(trialInds),1);
        trialMidMean = zeros(numel(trialInds),1);
        trialMidSTD = zeros(numel(trialInds),1);
        trialDeltaMean = zeros(numel(trialInds),1);
        trialDeltaSTD = zeros(numel(trialInds),1);

        for t = 1:numel(trialInds) % iterate through all trials and collect values in current ROI
            curMidVals = zeros(size(curPixels,1),1);
            curPreVals = zeros(size(curPixels,1),1);
            curSTDVals = zeros(size(curPixels,1),1);
            curPre = preStim{t};
            curMid = midStim{t};
%             curDelta = midStim{t}-preStim{t};

            for m = 1:size(curPixels,1)
                curPreVals(m) = curPre(curPixels(m,2),curPixels(m,1),curPixels(m,3));
                curMidVals(m) = curMid(curPixels(m,2),curPixels(m,1),curPixels(m,3));
                curVol(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
                curSTDVals(m) = standardDev(curPixels(m,2),curPixels(m,1),curPixels(m,3));
            end   
%             curThresh = stdMod*mean(curSTDVals);
            trialPreMean(t) = mean(curPreVals);
            trialMidMean(t) =mean(curMidVals);
            trialPreSTD(t) =std(curPreVals);
            trialMidSTD(t) =std(curMidVals);

            trialDeltaMean(t) =mean(curMidVals-curPreVals);
            trialDeltaSTD(t) =std(curMidVals-curPreVals);
        end
        roiTraces{r} = trialDeltaMean;

%         if(midRaster(r)==1) % only analyze regions that neurons were detected in
% 
%             % Generate figure with subplots examining intensity of ROI in
%             % trials analyzed at current current (he he)
%             curFig = figure();
%             subplot(2,2,1)
%             errorbar(trialPreMean,trialPreSTD)
%             hold on
%             errorbar(trialMidMean,trialMidSTD)
%             xlim([0,numel(trialMidMean)+1])
%             ylabel('Stimulation Flourescence')
%             xlabel('Trial')
%             legend({'Pre-Stim','Mid-Stim'})
%             title('Pre Vs Mid Intensity')
% 
% 
%             subplot(2,2,2)
%             ylabel('BaselineSTD')
%             stdMax = 1.2 * max((stdMod*standardDev),[],'all','omitnan');
%             zProjected = max((stdMod*standardDev),[],3,'omitnan');
%             imshow(zProjected/stdMax)
%             xlabel(strcat('Ave Threshold: ',num2str(mean(stdMod*standardDev,'all'))))
% 
%             subplot(2,2,3)
%             errorbar(trialDeltaMean,trialDeltaSTD)
%             hold on
%             yline(curThresh,'r--',strcat(num2str(stdMod),'x STD'));
%             yline(mean(trialDeltaMean),'g');
%             xlim([0,numel(trialDeltaMean)+1])
%             ylabel('Stimulation Flourescence Delta')
%             xlabel('Trial')
%             title('Stimulation Delta Intensity')
% 
%             subplot(2,2,4)
%             zProjectedBW = max(curVol,[],3,'omitnan');
%             boundaries = bwboundaries(zProjectedBW);
%             zProjected = max(midDelta,[],3,'omitnan');
%             imshow(zProjected/stdMax)
%             hold on
%             for k=1:numel(boundaries)
%                b = boundaries{k};
%                plot(b(:,2),b(:,1),'r','LineWidth',1);
%             end
%             title('ROI Analyzed')
% 
%             combineFileName = strcat(curTraceDir,'\ROI',num2str(r),'_Intensity.fig');
%             saveas(curFig,combineFileName);  
%             combineFileName = strcat(curTraceDir,'\ROI',num2str(r),'_Intensity.tif');
%             saveas(curFig,combineFileName);
%             close(curFig)
% 
%         end
    end
end










