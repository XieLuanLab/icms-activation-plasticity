% function thresholdMerging()
function thresholdMergingSept2025(daysThresholds83,daysThresholds92,daysThresholds93,daysThresholds98,daysThresholds100,daysThresholds101,contactPos83,contactPos92,contactPos93,contactPos98,contactPos100,contactPos101,sourceFolders83, sourceFolders92, sourceFolders93, sourceFolders98, sourceFolders100, sourceFolders101)
    
    sourceFolders = {sourceFolders83, sourceFolders92, sourceFolders93, sourceFolders98, sourceFolders100, sourceFolders101};
    % vertStep = 25; % electrode spacing information
    % contactPos = {contactPos83,contactPos92,contactPos93,contactPos98,contactPos100,contactPos101};
% contactsOrdered = [9,31,11,29,13,27,15,25,7,17,5,19,3,21,1,23,10,32,12,30,14,28,16,26,8,18,6,20,4,22,2,24];

%% pull variables from saved files
% [10;29]
% [9;11;12]
% [9;11;12]
% [4;5;6]
% [3;5;6]
% [10;11;12]

    % pathR{1} = 'Z:\Roy\2P Image Segmentation_12-20-24\2P Image Segmentation\raster83.mat';
    % pathR{2} = 'Z:\Roy\2P Image Segmentation_12-20-24\2P Image Segmentation\raster92.mat';
    % pathR{3} = 'Z:\Roy\2P Image Segmentation_12-20-24\2P Image Segmentation\raster93.mat';
    % pathR{4} = 'Z:\Roy\2P Image Segmentation_12-20-24\2P Image Segmentation\raster98.mat';
    % pathR{5} = 'Z:\Roy\2P Image Segmentation_12-20-24\2P Image Segmentation\raster100.mat';
    % pathR{6} = 'Z:\Roy\2P Image Segmentation_12-20-24\2P Image Segmentation\raster101.mat';
    % 
    % 
    % pathF{1} = 'Z:\Roy\AlternativeConcistencyDataMiss_ver6MISS83.mat';
    % pathF{2} = 'Z:\Roy\AlternativeConcistencyDataMiss_ver6MISS92.mat';
    % pathF{3} = 'Z:\Roy\AlternativeConcistencyDataMiss_ver6MISS93.mat';
    % pathF{4} = 'Z:\Roy\AlternativeConcistencyDataMiss_ver6MISS98.mat';
    % pathF{5} = 'Z:\Roy\AlternativeConcistencyDataMiss_ver6MISS100.mat';
    % pathF{6} = 'Z:\Roy\AlternativeConcistencyDataMiss_ver6MISS101.mat';
    % 
    % % 
    % % % % % Previously used n 4/16/24
    % pathC{1} = 'Z:\Roy\AlternativeConcistencyDataMiss_ver6MISS83.mat';
    % pathC{2} = 'Z:\Roy\AlternativeConcistencyDataMiss_ver6MISS92.mat';
    % pathC{3} = 'Z:\Roy\AlternativeConcistencyDataMiss_ver6MISS93.mat';
    % pathC{4} = 'Z:\Roy\AlternativeConcistencyDataMiss_ver6MISS98.mat';
    % pathC{5} = 'Z:\Roy\AlternativeConcistencyDataMiss_ver6MISS100.mat';
    % pathC{6} = 'Z:\Roy\AlternativeConcistencyDataMiss_ver6MISS101.mat';

    
    

    pathR{1} = 'Z:\xl_stimulation\Roy\2P Image Segmentation_12-20-24\2P Image Segmentation\raster83.mat';
    pathR{2} = 'Z:\xl_stimulation\Roy\2P Image Segmentation_12-20-24\2P Image Segmentation\raster92.mat';
    pathR{3} = 'Z:\xl_stimulation\Roy\2P Image Segmentation_12-20-24\2P Image Segmentation\raster93.mat';
    pathR{4} = 'Z:\xl_stimulation\Roy\2P Image Segmentation_12-20-24\2P Image Segmentation\raster98.mat';
    pathR{5} = 'Z:\xl_stimulation\Roy\2P Image Segmentation_12-20-24\2P Image Segmentation\raster100.mat';
    pathR{6} = 'Z:\xl_stimulation\Roy\2P Image Segmentation_12-20-24\2P Image Segmentation\raster101.mat';


    pathF{1} = 'Z:\xl_stimulation\Roy\AlternativeConcistencyDataMiss_ver6MISS83.mat';
    pathF{2} = 'Z:\xl_stimulation\Roy\AlternativeConcistencyDataMiss_ver6MISS92.mat';
    pathF{3} = 'Z:\xl_stimulation\Roy\AlternativeConcistencyDataMiss_ver6MISS93.mat';
    pathF{4} = 'Z:\xl_stimulation\Roy\AlternativeConcistencyDataMiss_ver6MISS98.mat';
    pathF{5} = 'Z:\xl_stimulation\Roy\AlternativeConcistencyDataMiss_ver6MISS100.mat';
    pathF{6} = 'Z:\xl_stimulation\Roy\AlternativeConcistencyDataMiss_ver6MISS101.mat';

    % 
    % % % % Previously used n 4/16/24
    pathC{1} = 'Z:\xl_stimulation\Roy\AlternativeConcistencyDataMiss_ver6MISS83.mat';
    pathC{2} = 'Z:\xl_stimulation\Roy\AlternativeConcistencyDataMiss_ver6MISS92.mat';
    pathC{3} = 'Z:\xl_stimulation\Roy\AlternativeConcistencyDataMiss_ver6MISS93.mat';
    pathC{4} = 'Z:\xl_stimulation\Roy\AlternativeConcistencyDataMiss_ver6MISS98.mat';
    pathC{5} = 'Z:\xl_stimulation\Roy\AlternativeConcistencyDataMiss_ver6MISS100.mat';
    pathC{6} = 'Z:\xl_stimulation\Roy\AlternativeConcistencyDataMiss_ver6MISS101.mat';


    anID = [83,92,93,98,100,101];
    uCurrAll = cell(numel(pathF),1);
    uCurrAll2 = cell(numel(pathF),1);
    uChanAll = cell(numel(pathF),1);
    daysTrainedAll = cell(numel(pathF),1);
    CurrChansAll = cell(numel(pathF),1);
    daysThresholdsAll = cell(numel(pathF),1);
    daysTrainedAll2 = cell(numel(pathF),1);
    CurrChansAll2 = cell(numel(pathF),1);
    
    allPopulationAll = cell(numel(pathF),1);
    allDensityAll = cell(numel(pathF),1);
    enclosingRadiusAll = cell(numel(pathF),1);
    populationDistanceAll = cell(numel(pathF),1);
    populationHistogramAll = cell(numel(pathF),1);
    groupedCurrChansAll = cell(numel(pathF),1);
    groupedCurrChansThreshAll = cell(numel(pathF),1);
    allRadiiThreshAll = cell(numel(pathF),1);
    allDensityThreshAll = cell(numel(pathF),1);
    populationDistancesAll = cell(numel(pathF),1);
    uChanCurrsRef = cell(numel(pathF),1);
    popDatesAll = cell(numel(pathF),1);
    groupedCurrChansBaseMeanAll = cell(numel(pathF),1);
    groupedCurrChansBaseSTDAll = cell(numel(pathF),1);
    groupedCurrChansBaseMeanUpper10All = cell(numel(pathF),1);
    groupedCurrChansThreshBaseMeanAll = cell(numel(pathF),1);
    groupedCurrChansThreshBaseSTDAll = cell(numel(pathF),1);
    groupedCurrChansThreshBaseMeanUpper10All = cell(numel(pathF),1);
    neighborComparisions = cell(numel(pathF),1);
    correlationDates = cell(numel(pathF),1);
    correlationMatricies = cell(numel(pathF),1);
    correlationMatricies3D = cell(numel(pathF),1);
    correlationMatriciesBase = cell(numel(pathF),1);
    correlationMatricies3DBase = cell(numel(pathF),1);
    diffMatricies = cell(numel(pathF),1);
    diffMatriciesSTD = cell(numel(pathF),1);
    diffMatriciesVals = cell(numel(pathF),1);
    diffMatricies3D = cell(numel(pathF),1);
    correlationMatriciesROI = cell(numel(pathF),1);
    correlationMatriciesROI3D = cell(numel(pathF),1);
    correlationMatriciesNeuronCount = cell(numel(pathF),1);
    correlationMatriciesNeuronCountNorm = cell(numel(pathF),1);
    uChanCurrsAll = cell(numel(pathF),1);
    sharedROIMappingAll = cell(numel(pathF),1);
    groupedCurrChansConsistencyAll = cell(numel(pathF),1);
    allRasters = cell(numel(pathF),1);
    allRegions = cell(numel(pathF),1);
    allVols =  cell(numel(pathF),1);
    for i = 1:numel(pathF)
        tempData = load(pathC{i});

%         tempData = load(pathF{i});
        uChanCurrsAll{i} = tempData.uChanCurrs;
        % uCurrAll{i} = tempData.uCurr;
        temp = tempData.uChanCurrs;
        uCurrAll{i} = unique(temp(:,2));
        uChanAll{i} = unique(tempData.uChanCurrs(:,1));
        if(i==1)
            temp = uChanAll{i};
            temp(2) = [];
            uChanAll{i} = temp;
        end
        CurrChansAll{i} = tempData.groupedCurrChans;
        daysTrainedAll{i} = tempData.daysTrained;
        allPopulationAll{i} = tempData.allPopulation;
        allDensityAll{i} = tempData.allDensity;
        enclosingRadiusAll{i} = tempData.enclosingRadius;
        populationDistanceAll{i} = tempData.populationDistance;
        populationHistogramAll{i} = tempData.populationHistogram;
        populationDistancesAll{i} = tempData.populationDistances;
        groupedCurrChansAll{i} = tempData.groupedCurrChans;
        groupedCurrChansThreshAll{i} = tempData.groupedCurrChansThresh;
        allRadiiThreshAll{i} = tempData.allRadiiThresh;
        allDensityThreshAll{i} = tempData.allDensityThresh;
        daysThresholdsAll{i} = tempData.daysThresholds;
        uChanCurrsRef{i} = tempData.uChanCurrs;
        popDatesAll{i} = tempData.popDates;
        
        groupedCurrChansBaseMeanAll{i} = tempData.groupedCurrChansBaseMean;
        groupedCurrChansBaseSTDAll{i} = tempData.groupedCurrChansBaseSTD;
        groupedCurrChansBaseMeanUpper10All{i} = tempData.groupedCurrChansBaseMeanUpper10;
        groupedCurrChansThreshBaseMeanAll{i} = tempData.groupedCurrChansThreshBaseMean;
        groupedCurrChansThreshBaseSTDAll{i} = tempData.groupedCurrChansThreshBaseSTD;
        groupedCurrChansThreshBaseMeanUpper10All{i} = tempData.groupedCurrChansThreshBaseMeanUpper10;
        
        
        
        tempData = load(pathC{i});
        uCurrAll2{i} = tempData.uCurr;
        CurrChansAll2{i} = tempData.groupedCurrChans;
        daysTrainedAll2{i} = tempData.daysTrained;
        neighborComparisions{i} = tempData.neighborComps;
        sharedROIMappingAll{i} = tempData.sharedROIMapping;
        groupedCurrChansConsistencyAll{i} = tempData.groupedCurrChansConsistency;
        
        correlationDates{i} = tempData.correlationDates;
        correlationMatricies{i} = tempData.correlationMatricies;
        correlationMatricies3D{i} = tempData.correlationMatricies3D;
        correlationMatriciesBase{i} = tempData.correlationMatriciesBase;
        correlationMatricies3DBase{i} = tempData.correlationMatricies3DBase;
        diffMatricies{i} = tempData.diffMatricies;
        diffMatriciesSTD{i} = tempData.diffMatriciesSTD;
        diffMatriciesVals{i} = tempData.diffMatriciesVals;
        diffMatricies3D{i} = tempData.diffMatricies3D;
        correlationMatriciesROI{i} = tempData.correlationMatriciesROI;
        correlationMatriciesROI3D{i} = tempData.correlationMatriciesROI3D;
        correlationMatriciesNeuronCount{i} = tempData.correlationMatriciesNeuronCount;
        correlationMatriciesNeuronCountNorm{i} = tempData.correlationMatriciesNeuronCountNorm;
        
        
        tempData = load(pathR{i});
        allRasters{i} = tempData.allRasters; % raster data not even used.....
        
        % pull region metrics from original datasets
        curSrc = sourceFolders{i};
        anmlVols = cell(numel(curSrc),1);
        anmlRegions = cell(numel(curSrc),1);
        for cs = 1:numel(curSrc)
            
            % find analysis folders
            selpath = curSrc{cs};
            idcs = strfind(selpath,'\');
            pullPath = selpath(1:idcs(end)-1);
            segListing = dir(pullPath);
            for s = 1:numel(segListing)
                segName = segListing(s).name;
                if(contains(segName, 'A1'))
                    trueSourceFolder = strcat(pullPath,'\',segName);
                end
            end
            analysisListing = dir(trueSourceFolder);
            analysisListing(1:2) = [];
            
            fldText=[];
            for fi = 1:numel(analysisListing)
                fldText = strcat(analysisListing(fi).folder,'\',analysisListing(fi).name);
                if(isfolder(fldText))
                    break
                end
            end
            
            roiData = load(strcat(fldText,'\ROI.mat')); % load data
            regions = roiData.regions;
            anmlRegions{cs} = regions;
            anmlVols{cs} = roiData.midActive;
        end
        allRegions{i} = anmlRegions;
        allVols{i} = anmlVols;
    end
    
    daysThresholdsAll{1} = daysThresholds83; % Mannually loaded day thresholds
    daysThresholdsAll{2} = daysThresholds92;
    daysThresholdsAll{3} = daysThresholds93;
    daysThresholdsAll{4} = daysThresholds98;
    daysThresholdsAll{5} = daysThresholds100;
    daysThresholdsAll{6} = daysThresholds101;
    

    dDiv = 7;


    
    currColors = {[1 0 0],[1 0 1],[0 0 1],[0 1 0],[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0 0 0],[1 1 0]}; % Hopefully these sets of colors and symbols will be enough
    currColorsAlpha = {[1 0 0 0.3],[1 0 1 0.3],[0 0 1 0.3],[0 1 0 0.3],[0.8500 0.3250 0.0980 0.3],[0.4660 0.6740 0.1880 0.3],[0 0 0 0.3],[1 1 0 0.3]}; % Hopefully these sets of colors and symbols will be enough

    
    
  
    %%
    uCurrOverallExaust = unique(cell2mat(uCurrAll));
    densityAll = cell(numel(anID),1); % Updated Density calculation, units adjusted for readability
    densityAllAgr = cell(numel(anID),1); % Updated Density calculation, units adjusted for readability

    % 
    % strelVol = cell(10,1);
    % for i = 1:10
    %     strelVol{i} = strel3d((95*(i))-0);
    % end


    % Process each animal
    for aInd = 1:numel(anID)
        uChan = uChanAll{aInd};
        uCurr = uCurrAll{aInd};
        daysTrained = daysTrainedAll{aInd};
        uChanCurrs = uChanCurrsRef{aInd};
        popDates = popDatesAll{aInd};
        populationDistances = populationDistancesAll{aInd};

        anVols = allVols{aInd};
        densityCur = densityAll{aInd};
        densityAgr = densityAllAgr{aInd};

        if(isempty(densityCur))
            % initialize variable if needed
            densityCur = NaN(numel(uCurrOverallExaust),3,numel(daysTrained),10);
            densityAgr = NaN(numel(uCurrOverallExaust),3,numel(daysTrained));
        end

        for cInd = 1:numel(uCurrOverallExaust) % process each current over weeks
            curInd = find(uCurr==uCurrOverallExaust(cInd));
            if(curInd>0) % if current current is present for given animal, continue analysis
                for chInd = 1:numel(uChan) % process each channel
                    if(aInd==5 && chInd==3)
                        continue;
                        % Skip ICMS 100 channel 6 (error was found in this
                        % dataset)
                    end
                    curCurr = uCurr(curInd);
                    for p = 1:size(uChanCurrs,1)
                        if(uChanCurrs(p,2)==curCurr) % if valid entry exists
                            if(uChanCurrs(p,1)==uChan(chInd)) % if valid entry exists
    
                                curDistances = populationDistances{p};
                                nDates = popDates{p};
                                curDays = daysTrained(nDates);  
                                curDays = ceil(curDays/dDiv);
                                for wk = 1:numel(curDays)
                                    d = curDays(wk);
                                    activeDists = curDistances{wk};
                                    activeDists(activeDists==0)=[];
                                    mergedDensity = zeros(10,1);
                                    if(sum(activeDists<1000)>1)
                                        for m = 1:10
                                            curcnt = activeDists>=(m-1)*100 & activeDists<m*100;
                                            mergedDensity(m) = sum(curcnt)/(((4/3)*pi()*(0.1*m)^3)-((4/3)*pi()*(0.1*(m-1))^3));  % currently calculated in mm rather than um
                                        end
                                        densityCur(cInd,chInd,wk,:) = mergedDensity;

                                        % curElecPos   (chInd)
                                        % curSegData =  false(size(anVols{wk},1),size(anVols{wk},2),size(anVols{wk},3)*25);
                                        % for k = 1:10
                                        %     % Create initial volume, must be large enough to encompase
                                        %     measureVol = false(size(curSegData);
                                        %     se = strelVol{k}; % Create masking structure and apply to volume
                                        % 
                                        %     seMask = getnhood(se);
                                        %     seMidX = floor(size(seMask,1)/2);
                                        %     seMidY = floor(size(seMask,2)/2);
                                        %     seMidZ = floor(size(seMask,3)/2);
                                        % 
                                        %     minXVol = max([1 round(curElecPos(1))-seMidX])+1;
                                        %     minYVol = max([1 round(curElecPos(2))-seMidY])+1;
                                        %     minZVol = max([1 round(curElecPos(3))-seMidZ])+1;
                                        %     maxXVol = min([size(curSegData,1) round(curElecPos(1))+seMidX]);
                                        %     maxYVol = min([size(curSegData,2) round(curElecPos(2))+seMidY]);
                                        %     maxZVol = min([size(curSegData,3) round(curElecPos(3))+seMidZ]);            
                                        % 
                                        % 
                                        %     xMinRange = round(curElecPos(1))- minXVol;
                                        %     xMaxRange = maxXVol - round(curElecPos(1));
                                        %     yMinRange = round(curElecPos(2))- minYVol;
                                        %     yMaxRange = maxYVol - round(curElecPos(2));
                                        %     zMinRange = round(curElecPos(3))- minZVol;
                                        %     zMaxRange = maxZVol - round(curElecPos(3));
                                        % 
                                        %     minXSE = seMidX-xMinRange;
                                        %     minYSE = seMidY-yMinRange;
                                        %     minZSE = seMidZ-zMinRange;
                                        %     maxXSE = seMidX+xMaxRange;
                                        %     maxYSE = seMidY+yMaxRange;
                                        %     maxZSE = seMidZ+zMaxRange;
                                        %     measureVol(minXVol:maxXVol, minYVol:maxYVol, minZVol:maxZVol) = seMask(minXSE:maxXSE, minYSE:maxYSE, minZSE:maxZSE);
                                        %     contactPos = contactPos{aInd};
                                        %     curElecPos = contactPos(uChan(chInd))
                                        %     if(k>1)
                                        %         % Subtract inner portion of volume to create shell
                                        %         % Create initial volume, must be large enough to encompase
                                        %         subtractVol = false(size(curSegData));
                                        %         se = strelVol{k-1}; % Create masking structure and apply to volume
                                        %         seMask = getnhood(se);
                                        %         seMidX = floor(size(seMask,1)/2);
                                        %         seMidY = floor(size(seMask,2)/2);
                                        %         seMidZ = floor(size(seMask,3)/2);
                                        % 
                                        %         minXVol = max([1 round(curElecPos(1))-seMidX])+1;
                                        %         minYVol = max([1 round(curElecPos(2))-seMidY])+1;
                                        %         minZVol = max([1 round(curElecPos(3))-seMidZ])+1;
                                        %         maxXVol = min([size(curSegData,1) round(curElecPos(1))+seMidX]);
                                        %         maxYVol = min([size(curSegData,2) round(curElecPos(2))+seMidY]);
                                        %         maxZVol = min([size(curSegData,3) round(curElecPos(3))+seMidZ]);            
                                        % 
                                        % 
                                        %         xMinRange = round(curElecPos(1))- minXVol;
                                        %         xMaxRange = maxXVol - round(curElecPos(1));
                                        %         yMinRange = round(curElecPos(2))- minYVol;
                                        %         yMaxRange = maxYVol - round(curElecPos(2));
                                        %         zMinRange = round(curElecPos(3))- minZVol;
                                        %         zMaxRange = maxZVol - round(curElecPos(3));
                                        % 
                                        %         minXSE = seMidX-xMinRange;
                                        %         minYSE = seMidY-yMinRange;
                                        %         minZSE = seMidZ-zMinRange;
                                        %         maxXSE = seMidX+xMaxRange;
                                        %         maxYSE = seMidY+yMaxRange;
                                        %         maxZSE = seMidZ+zMaxRange;
                                        %         subtractVol(minXVol:maxXVol, minYVol:maxYVol, minZVol:maxZVol) = seMask(minXSE:maxXSE, minYSE:maxYSE, minZSE:maxZSE);
                                        %         measureVol(subtractVol)=0; 
                                        %     end
                                        %     trimmedVol = sum(measureVol,'all') * (0.002187 * 0.002187 * 0.002); % multiply number of voxels by voxel volume
                                        % 
                                        %     % calculate 3D density for area surveyed
                                        %     if(trimmedVol==0)
                                        %         density3DFullTemp_SITE(siteCNT,k)=0;
                                        %     else
                                        %         density3DFullTemp_SITE(siteCNT,k) = density3D(k)/trimmedVol;
                                        %     end
                                        %     density3DFullTempCOM_SITE(siteCNT,k) = density3DCOM(k)/(((4/3)*pi()*(0.1*k)^3)-((4/3)*pi()*(0.1*(k-1))^3));
                                        % end
                                        % 
                                        % 



                                        curcnt = activeDists>=0 & activeDists<500;
                                        densityAgr(cInd,chInd,wk) = sum(curcnt)/((4/3)*pi()*0.5^3);  % currently calculated in mm rather than um
                                    end
                                end
                            end
                        end
                    end
                end 
            end
        end  
        densityAllAgr{aInd} = densityAgr;
        densityAll{aInd} = densityCur;

    end
    
    
    
 
    
    
    %% ACTIVATED POPULATION vs WEEKS ----------------------------------------------
    figure()
    hold on

    uCurrOverall = unique(cell2mat(uCurrAll));
    uCurrOverall(uCurrOverall>6)=[];
    uChanOverall = unique(cell2mat(uChanAll));
    daysTrainedOverall = unique(cell2mat(daysTrainedAll'));
    weeksTrainedOverall = unique(ceil(daysTrainedOverall/dDiv));
    for cInd = 1:numel(uCurrOverall) % process each current

        % Structures for averaging current specific response
        mergingTraces = NaN(numel(weeksTrainedOverall),numel(anID)*numel(uChanOverall));
        traceNum=0;
        
            
        % Process each animal
        for aInd = 1:numel(anID)
            uChan = uChanAll{aInd};
            uCurr = uCurrAll{aInd};
            groupedCurrChans = CurrChansAll{aInd};
            daysTrained = ceil(daysTrainedAll{aInd}/dDiv);

            curInd = find(uCurr==uCurrOverall(cInd));
            
            if(curInd>0) % if current current is present for given animal, continue analysis
                for chInd = 1:numel(uChan) % process each channel
                    if(aInd==5 && chInd==3)
                        continue;
                        % Skip ICMS 100 channel 6
                    end
                    neuCnt = groupedCurrChans(:,chInd,curInd);
                    curDays = daysTrained(~isnan(neuCnt));
                    neuCnt = neuCnt(~isnan(neuCnt));

                    % Only plot if there is data
                    if(~isempty(neuCnt))
                        scatter(curDays,neuCnt,[],currColors{cInd},'filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.8);
                        
                        % Add trace to current average
                        traceNum = traceNum+1;
                        alignedDays = zeros(numel(curDays),1);
                        for d = 1:numel(curDays)
                            alignedDays(d) = find(weeksTrainedOverall==curDays(d)); % Align days
                        end
                        mergingTraces(alignedDays,traceNum) = neuCnt;
                    end 
                end 
            end
        end
        
       
        
        % Average response for current and plot
        curAve = mean(mergingTraces,2,'omitnan');
        curDays = weeksTrainedOverall;
        curDays = curDays(~isnan(curAve));
        curAve = curAve(~isnan(curAve));

        p1 = plot(curDays',curAve);
        tcolor = currColors{cInd};
        p1.Color = tcolor(1:3);
        p1.LineWidth=2;

    end
    title('Neural Activation over Time')
    ylabel('Number of  Neurons')
    xlabel('Weeks of Training')

    lgd = legend('','Location', 'eastoutside');
    col_names = num2str(uCurrOverall);
    for j =1:length(col_names)
        plot([NaN NaN], [NaN NaN], 'Color', currColors{j}, 'DisplayName', strcat(col_names(j),' \muA'))
    end
    title(lgd,'Current')
  
    
    
    
    
    
    
    
    % ACTIVATED POPULATION vs WEEKS showing all source traces ----------------------------------------------
    figure()
    hold on

    uCurrOverall = unique(cell2mat(uCurrAll));
    uCurrOverall(uCurrOverall>6)=[];
    uChanOverall = unique(cell2mat(uChanAll));
    daysTrainedOverall = unique(cell2mat(daysTrainedAll'));
    weeksTrainedOverall = unique(ceil(daysTrainedOverall/dDiv));
    for cInd = 1:numel(uCurrOverall) % process each current

        % Structures for averaging current specific response
        mergingTraces = NaN(numel(weeksTrainedOverall),numel(anID)*numel(uChanOverall));
        traceNum=0;
        
            
        % Process each animal
        for aInd = 1:numel(anID)
            uChan = uChanAll{aInd};
            uCurr = uCurrAll{aInd};
            groupedCurrChans = CurrChansAll{aInd};
            daysTrained = ceil(daysTrainedAll{aInd}/dDiv);

            curInd = find(uCurr==uCurrOverall(cInd));
            
            if(curInd>0) % if current current is present for given animal, continue analysis
                for chInd = 1:numel(uChan) % process each channel
                    if(aInd==5 && chInd==3)
                        continue;
                        % Skip ICMS 100 channel 6
                    end
                    neuCnt = groupedCurrChans(:,chInd,curInd);
                    curDays = daysTrained(~isnan(neuCnt));
                    neuCnt = neuCnt(~isnan(neuCnt));

                    % Only plot if there is data
                    if(~isempty(neuCnt))
                        % merge datapoints from the same days
                        curDaysP = curDays;
                        neuCntP = neuCnt;
                        curDaysU = unique(curDaysP);
                        neuCntU = zeros(numel(curDaysU),1);
                        for u = 1:numel(curDaysU)
                            temp = neuCntP(curDaysP==curDaysU(u));
                            neuCntU(u) = mean(temp);
                        end
                        
                        plot(curDaysU,neuCntU,'-','Color',currColorsAlpha{cInd});
                        
                        % Add trace to current average
                        traceNum = traceNum+1;
                        alignedDays = zeros(numel(curDays),1);
                        for d = 1:numel(curDays)
                            alignedDays(d) = find(weeksTrainedOverall==curDays(d)); % Align days
                        end
                        mergingTraces(alignedDays,traceNum) = neuCnt;
                    end 
                end 
            end
        end
        
       
        
        % Average response for current and plot
        curAve = mean(mergingTraces,2,'omitnan');
        curDays = weeksTrainedOverall;
        curDays = curDays(~isnan(curAve));
        curAve = curAve(~isnan(curAve));

        p1 = plot(curDays',curAve);
        tcolor = currColors{cInd};
        p1.Color = tcolor(1:3);
        p1.LineWidth=2;

    end
    title('Neural Activation over Time')
    ylabel('Number of  Neurons')
    xlabel('Weeks of Training')

    lgd = legend('','Location', 'eastoutside');
    col_names = num2str(uCurrOverall);
    for j =1:length(col_names)
        plot([NaN NaN], [NaN NaN], 'Color', currColors{j}, 'DisplayName', strcat(col_names(j),' \muA'))
    end
    title(lgd,'Current')
  
    
    
    
    
    
    
    
    
    
     %% ACTIVATED POPULATION vs WEEKS ----------------------------------------------
    figure()
    hold on

    uCurrOverall = unique(cell2mat(uCurrAll));
    uCurrOverall(uCurrOverall>6)=[];
    uChanOverall = unique(cell2mat(uChanAll));
    daysTrainedOverall = unique(cell2mat(daysTrainedAll'));
    weeksTrainedOverall = unique(ceil(daysTrainedOverall/dDiv));
    for cInd = 1:numel(uCurrOverall) % process each current

        % Structures for averaging current specific response
        mergingTraces = NaN(numel(weeksTrainedOverall),numel(anID)*numel(uChanOverall));
        traceNum=0;
        
            
        % Process each animal
        for aInd = 1:numel(anID)
            uChan = uChanAll{aInd};
            uCurr = uCurrAll{aInd};
            groupedCurrChans = CurrChansAll{aInd};
            daysTrained = ceil(daysTrainedAll{aInd}/dDiv);

            curInd = find(uCurr==uCurrOverall(cInd));
            
            if(curInd>0) % if current current is present for given animal, continue analysis
                for chInd = 1:numel(uChan) % process each channel
                    if(aInd==5 && chInd==3)
                        continue;
                        % Skip ICMS 100 channel 6
                    end
                    neuCnt = groupedCurrChans(:,chInd,curInd);
                    curDays = daysTrained(~isnan(neuCnt));
                    neuCnt = neuCnt(~isnan(neuCnt));

                    % Only plot if there is data
                    if(~isempty(neuCnt))
                        
                        % Add trace to current average
                        traceNum = traceNum+1;
                        alignedDays = zeros(numel(curDays),1);
                        for d = 1:numel(curDays)
                            alignedDays(d) = find(weeksTrainedOverall==curDays(d)); % Align days
                        end
                        mergingTraces(alignedDays,traceNum) = neuCnt;
                    end 
                end 
            end
        end
        
       
        
        % Average response for current and plot
        curAve = mean(mergingTraces,2,'omitnan');
        curSTD = std(mergingTraces,[],2,'omitnan');
        curDays = weeksTrainedOverall;
        curDays = curDays(~isnan(curAve));
        curSTD = curSTD(~isnan(curAve));
        curAve = curAve(~isnan(curAve));

        
        p1 = errorbar(curDays',curAve,curSTD);
        tcolor = currColors{cInd};
        p1.Color = tcolor(1:3);
        p1.LineWidth=2;

    end
    title('Neural Activation over Time')
    ylabel('Number of  Neurons')
    xlabel('Weeks of Training')

    lgd = legend('','Location', 'eastoutside');
    col_names = num2str(uCurrOverall);
    for j =1:length(col_names)
        plot([NaN NaN], [NaN NaN], 'Color', currColors{j}, 'DisplayName', strcat(col_names(j),' \muA'))
    end
    title(lgd,'Current')
    
    
    
    
    
    
    
    
    
    %% ACTIVATED POPULATION vs WEEKS statistical analysis----------------------------------------------

    uCurrOverall = unique(cell2mat(uCurrAll));
    uCurrOverall(uCurrOverall>6)=[];
    for cInd = 1:numel(uCurrOverall) % process each current

        % Structures for averaging current specific response
        dpDays = [];
        dpCounts = [];
        % Process each animal
        for aInd = 1:numel(anID)
            uChan = uChanAll{aInd};
            uCurr = uCurrAll{aInd};
            groupedCurrChans = CurrChansAll{aInd};
            daysTrained = ceil(daysTrainedAll{aInd}/dDiv);

            curInd = find(uCurr==uCurrOverall(cInd));
            
            if(curInd>0) % if current current is present for given animal, continue analysis
                for chInd = 1:numel(uChan) % process each channel
                    if(aInd==5 && chInd==3)
                        continue;
                        % Skip ICMS 100 channel 6
                    end
                    neuCnt = groupedCurrChans(:,chInd,curInd);
                    curDays = daysTrained(~isnan(neuCnt));
                    neuCnt = neuCnt(~isnan(neuCnt));
                    neuCnt(curDays>4) = [];
                    curDays(curDays>4) = [];

                    % Only plot if there is data
                    if(~isempty(neuCnt))
                        dpDays = [dpDays; curDays'];
                        dpCounts = [dpCounts; neuCnt];
                        
                    end 
                end 
            end
        end
        
%         [p,tbl,stats] = anova1(dpCounts,dpDays);
%         results = multcompare(stats);

    end
 
    
    
    
    
    
    
    
    
    
    
    
    % ACTIVATED POPULATION vs WEEKS (Linear Regression) ----------------------------------------------
    figure()
    hold on

    uCurrOverall = unique(cell2mat(uCurrAll));
    uCurrOverall(uCurrOverall>6)=[];
    daysTrainedOverall = unique(cell2mat(daysTrainedAll'));
    for cInd = 1:numel(uCurrOverall) % process each current

        % Structures for averaging current specific response
        dpDays = [];
        dpCounts = [];
        
        % Process each animal
        for aInd = 1:numel(anID)
            uChan = uChanAll{aInd};
            uCurr = uCurrAll{aInd};
            groupedCurrChans = CurrChansAll{aInd};
            daysTrained = ceil(daysTrainedAll{aInd}/dDiv);

            curInd = find(uCurr==uCurrOverall(cInd));
            
            if(curInd>0) % if current current is present for given animal, continue analysis
                for chInd = 1:numel(uChan) % process each channel
                    if(aInd==5 && chInd==3)
                        continue;
                        % Skip ICMS 100 channel 6
                    end
                    neuCnt = groupedCurrChans(:,chInd,curInd);
                    curDays = daysTrained(~isnan(neuCnt));
                    neuCnt = neuCnt(~isnan(neuCnt));
                    
                    neuCnt(curDays>4) = [];
                    curDays(curDays>4) = [];

                    % Only plot if there is data
                    if(~isempty(neuCnt))
%                         scatter(curDays,neuCnt,[],currColors{cInd},'filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.8);
                        
                        dpDays = [dpDays; curDays'];
                        dpCounts = [dpCounts; neuCnt];
                    end 
                end 
            end
        end
        
        
        % Average response for current and plot
        x = unique(dpDays);
        p = polyfit(dpDays,dpCounts,1);
        f = polyval(p,x); 
        p1 = plot(x,f,'-');
        tcolor = currColors{cInd};
        p1.Color = tcolor(1:3);
        p1.LineWidth=2;

    end
    title('Neural Activation over Time')
    ylabel('Number of  Neurons')
    xlabel('Weeks of Training')

    lgd = legend('','Location', 'eastoutside');
    col_names = num2str(uCurrOverall);
    for j =1:length(col_names)
        plot([NaN NaN], [NaN NaN], 'Color', currColors{j}, 'DisplayName', strcat(col_names(j),' \muA'))
    end
    title(lgd,'Current')
    
    
 
    % ACTIVATED POPULATION vs WEEKS (Linear Regression - Statistical Analysis) ----------------------------------------------
    uCurrOverall = unique(cell2mat(uCurrAll));
    uCurrOverall(uCurrOverall>6)=[];
    for cInd = 1:numel(uCurrOverall) % process each current

        % Structures for averaging current specific response
        dpDays = [];
        dpCounts = [];
        
        % Process each animal
        for aInd = 1:numel(anID)
            uChan = uChanAll{aInd};
            uCurr = uCurrAll{aInd};
            groupedCurrChans = CurrChansAll{aInd};
            daysTrained = ceil(daysTrainedAll{aInd}/dDiv);

            curInd = find(uCurr==uCurrOverall(cInd));
            
            if(curInd>0) % if current current is present for given animal, continue analysis
                for chInd = 1:numel(uChan) % process each channel
                    if(aInd==5 && chInd==3)
                        continue;
                        % Skip ICMS 100 channel 6
                    end
                    neuCnt = groupedCurrChans(:,chInd,curInd);
                    curDays = daysTrained(~isnan(neuCnt));
                    neuCnt = neuCnt(~isnan(neuCnt));
                    

                    % Only plot if there is data
                    if(~isempty(neuCnt))
                        
                        dpDays = [dpDays; curDays'];
                        dpCounts = [dpCounts; neuCnt];
                    end 
                end 
            end
        end
%         figure()
%         
%         mdl = fitlm(dpDays,dpCounts);
%         anova(mdl,'summary')
%         plot(mdl)
%         title(strcat('Neural Activation over Time Current ',num2str(uCurrOverall(cInd))))
    end



    
    
    
    
%     
%     
%    % ACTIVATED POPULATION vs WEEKS Quantify temporal stats---------------------------------------------
%     uCurrOverall = unique(cell2mat(uCurrAll));
%     uCurrOverall(uCurrOverall>6)=[];
%     for cInd = 1:numel(uCurrOverall) % process each current
% 
%         allCounts = [];
%         g = [];
%             
%         % Process each animal
%         for aInd = 1:numel(anID)
%             uChan = uChanAll{aInd};
%             uCurr = uCurrAll{aInd};
%             groupedCurrChans = CurrChansAll{aInd};
%             daysTrained = ceil(daysTrainedAll{aInd}/dDiv);
%             curInd = find(uCurr==uCurrOverall(cInd));
%             if(curInd>0) % if current current is present for given animal, continue analysis
%                 for chInd = 1:numel(uChan) % process each channel
%                     if(aInd==5 && chInd==3)
%                         continue;
%                         % Skip ICMS 100 channel 6
%                     end
%                     neuCnt = groupedCurrChans(:,chInd,curInd);
%                     curDays = daysTrained(~isnan(neuCnt));
%                     neuCnt = neuCnt(~isnan(neuCnt));
% 
%                     % Only plot if there is data
%                     if(~isempty(neuCnt))
%                         scatter(curDays,neuCnt,[],currColors{cInd},'filled');
%                         g = [g; curDays'];
%                         allCounts = [allCounts;neuCnt];
%                     end 
%                 end 
%             end
%         end
%         
%         
%         % perform statistical analysis
%         figure()
%         [~,~,stats] = anova1(allCounts,g);
%         results = multcompare(stats);
%         title(strcat('Current:',num2str(uCurrOverall(cInd))))
%     end
%     
%     
%     % ACTIVATED POPULATION vs WEEKS Quantify current stats---------------------------------------------
%     uCurrOverall = unique(cell2mat(uCurrAll));
%     uCurrOverall(uCurrOverall>6)=[];
%     allCounts = [];
%     g = [];
%     g2 = [];
%     for cInd = 1:numel(uCurrOverall) % process each current
% 
%         % Process each animal
%         for aInd = 1:numel(anID)
%             uChan = uChanAll{aInd};
%             uCurr = uCurrAll{aInd};
%             groupedCurrChans = CurrChansAll{aInd};
%             daysTrained = ceil(daysTrainedAll{aInd}/dDiv);
%             curInd = find(uCurr==uCurrOverall(cInd));
%             if(curInd>0) % if current current is present for given animal, continue analysis
%                 for chInd = 1:numel(uChan) % process each channel
%                     if(aInd==5 && chInd==3)
%                         continue;
%                         % Skip ICMS 100 channel 6
%                     end
%                     neuCnt = groupedCurrChans(:,chInd,curInd);
%                     curDays = daysTrained(~isnan(neuCnt));
%                     neuCnt = neuCnt(~isnan(neuCnt));
% 
%                     % Only plot if there is data
%                     if(~isempty(neuCnt))
% %                         scatter(curDays,neuCnt,[],currColors{cInd},'filled');
%                         curG = cInd*ones(numel(neuCnt),1);
%                         g = [g; curDays'];
%                         g2 = [g2; curG];
%                         allCounts = [allCounts;neuCnt];
%                     end 
%                 end 
%             end
%         end
%     end   
%     % perform statistical analysis
% %     figure()
% %     [~,~,stats] = anova1(allCounts,g);
% %     results = multcompare(stats);
%   
%     figure()
%     [~,~,stats] = anovan(allCounts,{g,g2});
%     results = multcompare(stats);
%     
%     
    
    
    
    
    
    
    
    
    % ACTIVATED POPULATION Vs. Threshold -------------------------------
    figure()
    hold on
    uCurrOverall = unique(cell2mat(uCurrAll));
    uCurrOverall(uCurrOverall>6)=[];
    uChanOverall = unique(cell2mat(uChanAll));
    
    allDTA = [];
    for a = 1:numel(daysThresholdsAll)
        temp = daysThresholdsAll{a};
        temp = reshape(temp,1,[]);
        allDTA = [allDTA, temp];
    end
    daysThresholdOverall = unique(round(allDTA));
%     daysThresholdOverall = unique(unique(round(cell2mat(daysThresholdsAll'))));
    daysThresholdOverall(daysThresholdOverall==0)=[];
    daysThresholdOverall(daysThresholdOverall>7)=[];
    threshDensityData = NaN(numel(uCurrOverall),numel(daysThresholdOverall));
    threshDensityDataCell = cell(numel(uCurrOverall),numel(daysThresholdOverall));

    for cInd = 1:numel(uCurrOverall) % process each current

        % Structures for averaging current specific response
        mergingTraces = NaN(numel(daysThresholdOverall),numel(anID)*numel(uChanOverall));
        traceNum=0;

        % Process each animal
        for aInd = 1:numel(anID)
            uChan = uChanAll{aInd};
            uCurr = uCurrAll{aInd};
%             groupedCurrChans = CurrChansAll{aInd};
            daysThresholds = round(daysThresholdsAll{aInd});
            uThresh = unique(reshape(daysThresholds,1,[]));
            uThresh(uThresh==0)=[]; % remove thresholds at 0    
            groupedCurrChansThresh = groupedCurrChansThreshAll{aInd};            

            curInd = find(uCurr==uCurrOverall(cInd));
            if(curInd>0) % if current current is present for given animal, continue analysis
                for chInd = 1:numel(uChan) % process each channel
                    if(aInd==5 && chInd==3)
                        continue;
                        % Skip ICMS 100 channel 6
                    end
                    neuCnt = groupedCurrChansThresh(:,chInd,curInd);
%                     curDayThresh = daysThresholds(chInd,:);
%                     curDays = curDayThresh(~isnan(neuCnt));
                    curDays = uThresh(~isnan(neuCnt));
                    neuCnt = neuCnt(~isnan(neuCnt));
                    
                    % Only plot if there is data
                    if(~isempty(neuCnt))
                        neuCnt(curDays>7)=[];
                        curDays(curDays>7)=[];
                        
                        scatter(curDays,neuCnt,[],currColors{cInd},'filled');

                        % Add trace to current average
                        traceNum = traceNum+1;
                        alignedDays = zeros(numel(curDays),1);
                        for d = 1:numel(curDays)
                            alignedDays(d) = find(daysThresholdOverall==curDays(d)); % Align thresholds
                        end
                        mergingTraces(alignedDays,traceNum) = neuCnt;
                    end 
                end 
            end
        end
        % Average response for current and plot
        curAve = mean(mergingTraces,2,'omitnan');
        threshDensityData(cInd,:) = curAve;
        curDays = daysThresholdOverall;
        curDays = curDays(~isnan(curAve));
        curAve = curAve(~isnan(curAve));

        p1 = plot(curDays',curAve);
        tcolor = currColors{cInd};
        p1.Color = tcolor(1:3);
        p1.LineWidth=2;
    end
    title('Neural Activation over Time')
    ylabel('Number of  Neurons')
    xlabel('Detection Threshold (\muA)')

    lgd = legend('','Location', 'eastoutside');
    col_names = num2str(uCurrOverall);
    for j =1:length(col_names)
        plot([NaN NaN], [NaN NaN], 'Color', currColors{j}, 'DisplayName', strcat(col_names(j),' \muA'))
    end
    title(lgd,'Current')
  
    
    
    
    
    
    
   
    
    
 % ACTIVATED POPULATION Vs. Threshold Trendline-------------------------------
    figure()
    hold on
    uCurrOverall = unique(cell2mat(uCurrAll));
    uCurrOverall(uCurrOverall>6)=[];
    uChanOverall = unique(cell2mat(uChanAll));
    
    allDTA = [];
    for a = 1:numel(daysThresholdsAll)
        temp = daysThresholdsAll{a};
        temp = reshape(temp,1,[]);
        allDTA = [allDTA, temp];
    end
    daysThresholdOverall = unique(round(allDTA));
%     daysThresholdOverall = unique(unique(round(cell2mat(daysThresholdsAll'))));
    daysThresholdOverall(daysThresholdOverall==0)=[];
    daysThresholdOverall(daysThresholdOverall>7)=[];
%     threshDensityData = NaN(numel(uCurrOverall),numel(daysThresholdOverall));

    for cInd = 1:numel(uCurrOverall) % process each current

        % Structures for averaging current specific response
        dpDays = [];
        dpCounts = [];

        % Process each animal
        for aInd = 1:numel(anID)
            uChan = uChanAll{aInd};
            uCurr = uCurrAll{aInd};
%             groupedCurrChans = CurrChansAll{aInd};
            daysThresholds = round(daysThresholdsAll{aInd});
            uThresh = unique(reshape(daysThresholds,1,[]));
            uThresh(uThresh==0)=[]; % remove thresholds at 0    
            groupedCurrChansThresh = groupedCurrChansThreshAll{aInd};            

            curInd = find(uCurr==uCurrOverall(cInd));
            if(curInd>0) % if current current is present for given animal, continue analysis
                for chInd = 1:numel(uChan) % process each channel
                    if(aInd==5 && chInd==3)
                        continue;
                        % Skip ICMS 100 channel 6
                    end
                    neuCnt = groupedCurrChansThresh(:,chInd,curInd);
%                     curDayThresh = daysThresholds(chInd,:);
%                     curDays = curDayThresh(~isnan(neuCnt));
                    curDays = uThresh(~isnan(neuCnt));
                    neuCnt = neuCnt(~isnan(neuCnt));
                    
                    % Only plot if there is data
                    if(~isempty(neuCnt))
                        neuCnt(curDays>7)=[];
                        curDays(curDays>7)=[];
                        
                        scatter(curDays,neuCnt,[],currColors{cInd},'filled');
                        
                        dpDays = [dpDays; curDays'];
                        dpCounts = [dpCounts;neuCnt];
                    end 
                end 
            end
        end
        
        % Average response for current and plot
        x = unique(dpDays);
        p = polyfit(dpDays,dpCounts,1);
        f = polyval(p,x); 
        p1 = plot(x,f,'-');
        tcolor = currColors{cInd};
        p1.Color = tcolor(1:3);
        p1.LineWidth=2;
    end
    title('Neural Activation over Time')
    ylabel('Number of  Neurons')
    xlabel('Detection Threshold (\muA)')

    lgd = legend('','Location', 'eastoutside');
    col_names = num2str(uCurrOverall);
    for j =1:length(col_names)
        plot([NaN NaN], [NaN NaN], 'Color', currColors{j}, 'DisplayName', strcat(col_names(j),' \muA'))
    end
    title(lgd,'Current')
    
    
    
    
    
%     
%     
%     
%  % ACTIVATED POPULATION Vs. Threshold Statistics-------------------------------
% 
%     uCurrOverall = unique(cell2mat(uCurrAll));
%     uCurrOverall(uCurrOverall>6)=[];
%     
%     allDTA = [];
%     for a = 1:numel(daysThresholdsAll)
%         temp = daysThresholdsAll{a};
%         temp = reshape(temp,1,[]);
%         allDTA = [allDTA, temp];
%     end
%     daysThresholdOverall = unique(round(allDTA));
% %     daysThresholdOverall = unique(unique(round(cell2mat(daysThresholdsAll'))));
%     daysThresholdOverall(daysThresholdOverall==0)=[];
%     daysThresholdOverall(daysThresholdOverall>7)=[];
% 
%     for cInd = 1:numel(uCurrOverall) % process each current
% 
%         % Structures for averaging current specific response
%         dpDays = [];
%         dpCounts = [];
% 
%         % Process each animal
%         for aInd = 1:numel(anID)
%             uChan = uChanAll{aInd};
%             uCurr = uCurrAll{aInd};
% %             groupedCurrChans = CurrChansAll{aInd};
%             daysThresholds = round(daysThresholdsAll{aInd});
%             uThresh = unique(reshape(daysThresholds,1,[]));
%             uThresh(uThresh==0)=[]; % remove thresholds at 0    
%             groupedCurrChansThresh = groupedCurrChansThreshAll{aInd};            
% 
%             curInd = find(uCurr==uCurrOverall(cInd));
%             if(curInd>0) % if current current is present for given animal, continue analysis
%                 for chInd = 1:numel(uChan) % process each channel
%                     if(aInd==5 && chInd==3)
%                         continue;
%                         % Skip ICMS 100 channel 6
%                     end
%                     neuCnt = groupedCurrChansThresh(:,chInd,curInd);
% %                     curDayThresh = daysThresholds(chInd,:);
% %                     curDays = curDayThresh(~isnan(neuCnt));
%                     curDays = uThresh(~isnan(neuCnt));
%                     neuCnt = neuCnt(~isnan(neuCnt));
%                     
%                     % Only plot if there is data
%                     if(~isempty(neuCnt))
%                         neuCnt(curDays>7)=[];
%                         curDays(curDays>7)=[];
%                         
%                         
%                         dpDays = [dpDays; curDays'];
%                         dpCounts = [dpCounts;neuCnt];
%                     end 
%                 end 
%             end
%         end
%         
% %         mdl = fitlm(dpDays,dpCounts);
% %         anova(mdl,'summary')
% %         plot(mdl)
% %         title(strcat('Neural Activation over Threshold - Current ',num2str(uCurrOverall(cInd))))
% %     
%         
%         
%     end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
      
    % Plot activation density with respect to threshold (traces)
    currColorsT  = jet(numel(daysThresholdOverall));
    figure
    hold on
    for t = 1:numel(daysThresholdOverall)
        threshAves = squeeze(threshDensityData(:,t));
        threshAves(uCurrOverall<3)=NaN;
        p1 = plot(uCurrOverall,threshAves,'-o');
        tcolor = currColorsT(t,:);
        p1.Color = tcolor(1:3);
        p1.LineWidth=2;
    end
    ylabel('Number of Neurons ')
    xlabel('Stimulating Current (\muA)')
    title('Neural Activation Across Currents & Thresholds')


    lgd = legend('','Location', 'eastoutside');
    col_names = num2str(daysThresholdOverall');
    for j =1:length(col_names)
        plot([NaN NaN], [NaN NaN], 'Color', currColorsT(j,:), 'DisplayName', col_names(j))
    end
    title(lgd,'Detection Threshold (\muA)')
    
    
    
    
    
    
    
    
    
    
    
    % Neighboring Contact Co-activation vs. Time ----------------------------------------------
    figure()
    hold on
    uCurrOverall = unique(cell2mat(uCurrAll));
    uCurrOverall(uCurrOverall>6)=[];
    uCurrOverall(uCurrOverall<2)=[];
    uChanOverall = unique(cell2mat(uChanAll));
    daysTrainedOverall = unique(ceil(cell2mat(daysTrainedAll2')/dDiv));
    for cInd = 2:numel(uCurrOverall) % process each current

        % Structures for averaging current specific response
        mergingTraces = NaN(numel(daysTrainedOverall),numel(anID)*numel(uChanOverall));
        traceNum=0;

        % Process each animal
        for aInd = 1:numel(anID)
            uCurr = uCurrAll{aInd};
            curNeighData = neighborComparisions{aInd};
            weeksTrained = ceil(daysTrainedAll2{aInd}/dDiv);

            curInd = find(uCurr==uCurrOverall(cInd));
            if(curInd>0) % if current current is present for given animal, continue analysis
                
                for cc = 1:3 % for each pair of compared channels report results
                    curOverlap = 100*curNeighData(:,cc,3,curInd); % Percentage of overlapping regions
                    curDays = weeksTrained;
                    curDays(isnan(curOverlap))=[];
                    curOverlap(isnan(curOverlap))=[];
                    if(~isempty(curOverlap))
                        scatter(curDays,curOverlap,[],currColors{cInd},'filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);

                        % Add trace to current average
                        traceNum = traceNum+1;
                        alignedDays = zeros(numel(curDays),1);
                        for d = 1:numel(curDays)
                            alignedDays(d) = find(daysTrainedOverall==curDays(d)); % Align days
                        end
                        mergingTraces(alignedDays,traceNum) = curOverlap;
                    end
                end
            end
        end
        % Average response for current and plot
        curAve = mean(mergingTraces,2,'omitnan');
        curDays = daysTrainedOverall;
        curDays = curDays(~isnan(curAve));
        curAve = curAve(~isnan(curAve));

        curDays(curDays>5)=NaN;
        
        p1 = plot(curDays',curAve);
        tcolor = currColors{cInd};
        p1.Color = tcolor(1:3);
        p1.LineWidth=2;
    end
    title('Overlapping Regions Over Time')
    ylabel('Percentage of Overlapping Active Neurons')
    xlabel('Weeks of Training')
    
    lgd = legend('','Location', 'eastoutside');
    col_names = num2str(uCurrOverall);
    for j =1:length(col_names)
        plot([NaN NaN], [NaN NaN], 'Color', currColors{j}, 'DisplayName', strcat(col_names(j),' \muA'))
    end
    title(lgd,'Current')
    ylim([0 100])
    


% % 
% % 
% %     % Co-activation needs to be calculated between contact at their
% %     % threshold on given dates not between fixed currents... have to update
% % 
% %     % Calculate the number of neurons relative between neighboring
% %     % contact sites - Comparisons are only made on the same day between
% %     % currents that are closest to threshold
% %     neighborThreshCompsAll = cell(numel(anID),1);
% %     for aInd = 1:numel(anID)
% %         daysThresholds = daysThresholdsAll{aInd};
% %         if(aInd==1)
% %             daysThresholds = daysThresholds([1,3],:);
% %         end
% %         daysTrained = ceil(daysTrainedAll2{aInd}/dDiv);
% %         uChan = uChanAll{aInd};
% %         uCurr = uCurrAll2{aInd};
% %         curRasters = allRasters{aInd};
% %         curChanCurrs = uChanCurrsRef{aInd};
% %         curPopDates = popDatesAll{aInd};
% % %         groupedCurrChans = groupedCurrChansAll{aInd};
% % 
% %         neighborThreshComps = NaN(numel(daysTrained),3);
% %         for s = 1:numel(daysTrained) % iterate through each session
% %             for c1 = 1:numel(uChan) % 
% %                 % Load channel specific data
% %                 c1Thresh = daysThresholds(c1,s);
% %                 cInd1 = find(uCurr==round(c1Thresh));
% % 
% %                 if(isempty(cInd1)) % no data, go to next iteration
% %                     continue
% %                 end
% % 
% %                 % Find associated raster dataset
% %                 m1 = (curChanCurrs(:,1) == uChan(c1));
% %                 m2 = (curChanCurrs(:,2) == uCurr(cInd1));
% %                 mInd = m1 & m2;
% %                 curDates = curPopDates{mInd};
% %                 dateInd = find(curDates==s);
% %                 rasterSets = curRasters{mInd};
% % 
% %                 if(isempty(dateInd)) % no data, go to next iteration
% %                     continue
% %                 end
% %                 c1Raster = rasterSets{dateInd};
% % 
% % 
% % 
% %                 % compare to other channels
% %                 for c2 = c1+1:numel(uChan)
% % 
% %                     % Load channel specific data
% %                     c2Thresh = daysThresholds(c2,s);
% %                     cInd2 = find(uCurr==round(c2Thresh));
% % 
% %                     if(isempty(cInd2)) % no data, go to next iteration
% %                         continue
% %                     end
% % 
% %                     % Find associated raster dataset
% %                     m1 = (curChanCurrs(:,1) == uChan(c2));
% %                     m2 = (curChanCurrs(:,2) == uCurr(cInd2));
% %                     mInd2 = m1 & m2;
% %                     curDates = curPopDates{mInd2};
% %                     dateInd = find(curDates==s);
% %                     rasterSets = curRasters{mInd2};
% % 
% %                     if(isempty(dateInd)) % no data, go to next iteration
% %                         continue
% %                     end
% %                     c2Raster = rasterSets{dateInd};
% % 
% %                     % perform comparison
% %                     rasterOverlap = c1Raster;
% %                     rasterOverlap(c2Raster==0)=0;
% %                     sharedRaster = c1Raster | c2Raster;
% %                     n = c1 + c2 - 2; % Silly but should work for small case of 3 channels max
% %                     neighborThreshComps(s,n) = sum(rasterOverlap)/sum(sharedRaster); % percentage of shared regions
% %                 end
% %             end
% %         end
% %         neighborThreshCompsAll{aInd} = neighborThreshComps;
% %     end
% % 
% % 
    % % 
    % % 
    % % % Neighboring Contact Co-activation at Threshold over time ----------------------------------------------
    % % % Structures for averaging current specific response
    % % mergingTraces = NaN(numel(weeksTrainedOverall),numel(anID)*3);
    % % mt=0;
    % % figure()
    % % hold on
    % % 
    % % % Process each animal
    % % for aInd = 1:numel(anID)
    % %     for c = 1:3
    % %         curNeighData = neighborThreshCompsAll{aInd};
    % %         weeksTrained = ceil(daysTrainedAll2{aInd}./dDiv);
    % % 
    % %         curOverlap = squeeze(100*curNeighData(:,c)); % Percentage of overlapping regions
    % % 
    % %         % Plot points for activation at threshold for current animal/channel
    % %         s = scatter(weeksTrained,curOverlap,[],'b','filled');
    % %         s.MarkerFaceAlpha = 0.5;
    % % 
    % %         % Save data to merging structure to examine effect over all
    % %         % weeks, animals, and channels
    % %         for d = 1:numel(weeksTrained)
    % %             mt = mt+1;
    % %             mergingTraces(weeksTrainedOverall==weeksTrained(d),mt) = curOverlap(d);
    % %         end
    % %     end
    % % end
    % % 
    % % % Average response for each week and plot neural population
    % % mergingTraces(mergingTraces==0)=NaN;
    % % curAve = mean(mergingTraces,2,'omitnan');
    % % curWeeks = weeksTrainedOverall;
    % % curWeeks(curWeeks>4)=NaN;
    % % plot(curWeeks',curAve,'b','LineWidth',2);
    % % 
    % % title('Overlapping Neural Activation at Threshold')
    % % ylabel('Percentage of Overlapping Active Neurons')
    % % xlabel('Weeks of Training')
    % % 
    % % 

    
    
    
    
    
    
    
    
    
    
    
%     % NEURAL ACTIVATION AT THREHOLD OVER TIME ---------------------------------------------- 
%     figure()
%     hold on
%     
%     % Structures for averaging current specific response
%     daysTrainedOverall = unique(cell2mat(daysTrainedAll'));
%     weeksTrainedOverall = unique(ceil(daysTrainedOverall/dDiv));
%     mergingTraces = NaN(numel(weeksTrainedOverall),numel(anID)*3);
%     mergingThresholds = NaN(numel(weeksTrainedOverall),numel(anID)*3);
%     mt=0;
% 
%     % Process each animal
%     for aInd = 1:numel(anID)
%         uChan = uChanAll{aInd};
%         uCurr = uCurrAll{aInd};
%             curNeighData = neighborComparisions{aInd};
%         weeksTrained = ceil(daysTrainedAll{aInd}/dDiv);
%         daysThresholds = daysThresholdsAll{aInd};
%         daysThresholds(daysThresholds==0)=NaN;
%         for chInd = 1:numel(uChan) % Process each channel
%             if(aInd==5 && chInd==3)
%                 continue;
%                 % Skip ICMS 100 channel 6
%             end
%             
%             
%             % find theresponse at each day's detection threshold  and
%             % plot, nearest neighbor
%             threshNeuCnts = NaN(size(chanThresh));
%             for d = 1:numel(weeksTrained) % for each session
%                 ccInd = find(uCurr == round(chanThresh(d)));
%                 if(~isempty(ccInd))
%                     threshNeuCnts(d) = chanNeuCnt(d,ccInd);
%                 end
%                 
%                 % comparisions of activation need to be re-performed at
%                 % each channel's threshold currents for each day
%                 
%                 % Aquire raster for each channel on the given day for the
%                 % threshold current
%                 
%             end
%             
%             for cc = 1:3 % for each pair of compared channels report results
%                 curOverlap = 100*curNeighData(:,cc,3,curInd); % Percentage of overlapping regions
%                 curDays = weeksTrained;
%                 curDays(isnan(curOverlap))=[];
%                 curOverlap(isnan(curOverlap))=[];
%                 if(~isempty(curOverlap))
%                     scatter(curDays,curOverlap,[],currColors{cInd},'filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
% 
%                     % Add trace to current average
%                     traceNum = traceNum+1;
%                     alignedDays = zeros(numel(curDays),1);
%                     for d = 1:numel(curDays)
%                         alignedDays(d) = find(daysTrainedOverall==curDays(d)); % Align days
%                     end
%                     mergingTraces(alignedDays,traceNum) = curOverlap;
%                 end
%             end
% 
%             
%             
%             curNeighData
%             chanNeuCnt = squeeze(groupedCurrChans(:,chInd,:));
%             chanThresh = daysThresholds(chInd,:);
%             threshNeuCnts = NaN(size(chanThresh));
% 
%             % calculate response at each day's detection threshold  and
%             % plot, nearest neighbor
%             for d = 1:numel(weeksTrained)
%                 ccInd = find(uCurr == round(chanThresh(d)));
%                 if(~isempty(ccInd))
%                     threshNeuCnts(d) = chanNeuCnt(d,ccInd);
%                 end
%             end
%              
%             % Plot points for activation at threshold
%             yyaxis left
% %             plot(weeksTrained,threshNeuCnts)
% 
%             s = scatter(weeksTrained,threshNeuCnts,[],'b','filled');
%             s.MarkerFaceAlpha = 0.5;
%             
%             
%             % Plot points for threshold
%             yyaxis right
%             chanThresh(isnan(threshNeuCnts))=NaN;
%             s = scatter(weeksTrained+0.1,chanThresh,[],[0.8500 0.3250 0.0980],'filled');
%             s.MarkerFaceAlpha = 0.5;
%             
%             % Save data to merging structure to examine effect over all
%             % weeks, animals, and channels
%             for d = 1:numel(weeksTrained)
%                 mt = mt+1;
%                 mergingTraces(weeksTrainedOverall==weeksTrained(d),mt) = threshNeuCnts(d);
%                 mergingThresholds(weeksTrainedOverall==weeksTrained(d),mt) = chanThresh(d);
%             end
%         end
%     end
%     
%     % Average response for each week and plot neural population
%     yyaxis left
%     mergingTraces(mergingTraces==0)=NaN;
%     curAve = mean(mergingTraces,2,'omitnan');
%     curWeeks = weeksTrainedOverall;
%     curWeeks(curWeeks>8)=NaN;
% %     curWeeks(curWeeks<1)=NaN;
% 
%     
%     
%     plot(curWeeks',curAve,'b','LineWidth',2);
% 
%     title('Neural Activation at Threshold - NN')
%     ylabel('Number of  Neurons')
%     xlabel('Weeks of Training')
%     
%     % Change axis and plot average threshold
%     yyaxis right
%     mergingThresholds(mergingThresholds==0)=NaN;
%     curThreshAve = mean(mergingThresholds,2,'omitnan');
%     plot(curWeeks',curThreshAve,'Color',[0.8500 0.3250 0.0980],'LineWidth',2);
%     ylabel('Detection Threshold (\muA)')
%     
%     
%     
%     
%     
%     
%     
    
    
    
    
    
    
    

% 
% 
% % Neighboring Contact Co-activation vs. Time vs intersite distance----------------------------------------------
%     for d = 60:60:180
%         figure()
%         hold on
%         uCurrOverall = unique(cell2mat(uCurrAll));
%         uCurrOverall(uCurrOverall>6)=[];
%         uCurrOverall(uCurrOverall<2)=[];
%         uChanOverall = unique(cell2mat(uChanAll));
%         daysTrainedOverall = unique(ceil(cell2mat(daysTrainedAll')/dDiv));
%         for cInd = 2:numel(uCurrOverall) % process each current
% 
%             % Structures for averaging current specific response
%             mergingTraces = NaN(numel(daysTrainedOverall),numel(anID)*numel(uChanOverall));
%             traceNum=0;
% 
%             % Process each animal
%             for aInd = 1:numel(anID)
%                 uCurr = uCurrAll{aInd};
%                 curNeighData = neighborComparisions{aInd};
%                 weeksTrained = ceil(daysTrainedAll{aInd}/dDiv);
% 
%                 curInd = find(uCurr==uCurrOverall(cInd));
%                 if(curInd>0) % if current current is present for given animal, continue analysis
% 
%                     for cc = 1:3 % for each pair of compared channels report results
%                         currDistance = curNeighData(:,cc,1,curInd);
%                         curOverlap = 100*curNeighData(:,cc,3,curInd); % Percentage of overlapping regions
%                         curDays = weeksTrained;
%                         curOverlap(currDistance~=d)=NaN;
%                         curDays(isnan(curOverlap))=[];
%                         curOverlap(isnan(curOverlap))=[];
%                         if(~isempty(curOverlap))
%                             scatter(curDays,curOverlap,[],currColors{cInd},'filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
% 
%                             % Add trace to current average
%                             traceNum = traceNum+1;
%                             alignedDays = zeros(numel(curDays),1);
%                             for d2 = 1:numel(curDays)
%                                 alignedDays(d2) = find(daysTrainedOverall==curDays(d2)); % Align days
%                             end
%                             mergingTraces(alignedDays,traceNum) = curOverlap;
%                         end
%                     end
%                 end
%             end
%                         
%             
%             % Average response for current and plot
%             curAve = mean(mergingTraces,2,'omitnan');
%             curDays = daysTrainedOverall;
%             curDays = curDays(~isnan(curAve));
%             curAve = curAve(~isnan(curAve));
% 
%             if(~isempty(curAve))
%                 p1 = plot(curDays',curAve);
%                 tcolor = currColors{cInd};
%                 p1.Color = tcolor(1:3);
%                 p1.LineWidth=2;
%             end
%         end
%         title(strcat('Overlapping Regions Over Time at',{' '},num2str(d),' \mum'))
%         ylabel('Percentage of Overlapping Active Neurons')
%         xlabel('Weeks of Training')
% 
%         lgd = legend('','Location', 'eastoutside');
%         col_names = num2str(uCurrOverall);
%         for j =1:length(col_names)
%             plot([NaN NaN], [NaN NaN], 'Color', currColors{j}, 'DisplayName', strcat(col_names(j),' \muA'))
%         end
%         title(lgd,'Current')
%         ylim([0 100])
%     end
%     
%     
%     
    
    
    
    
    
    %% ENCLOSING RADIUS vs. Time, seperated ----------------------------------------------
    figure()
    hold on
    uCurrOverall = unique(cell2mat(uCurrAll));
    uCurrOverall(uCurrOverall>6)=[];
    uCurrOverall(uCurrOverall<2)=[];
    uChanOverall = unique(cell2mat(uChanAll));
    daysTrainedOverall = unique(ceil(cell2mat(daysTrainedAll')/dDiv));
    
    for cInd = 3%2:numel(uCurrOverall) % process each current

        % Structures for averaging current specific response
        mergingTraces = NaN(numel(daysTrainedOverall),numel(anID)*numel(uChanOverall));
        traceNum=0;

        % Process each animal
        for aInd = 1:numel(anID)
            uChan = uChanAll{aInd};
            uCurr = uCurrAll{aInd};
            groupedRadii = enclosingRadiusAll{aInd};
            daysTrained = ceil(daysTrainedAll{aInd}/dDiv);
            uChanCurrs = uChanCurrsRef{aInd};
            popDates = popDatesAll{aInd};

            curInd = find(uCurr==uCurrOverall(cInd));
            if(curInd>0) % if current current is present for given animal, continue analysis
                for chInd = 1:numel(uChan) % process each channel
                    if(aInd==5 && chInd==3)
                        continue;
                        % Skip ICMS 100 channel 6
                    end
                    curCh = uChan(chInd);
                    curCurr = uCurr(curInd);
                    for p = 1:size(uChanCurrs,1)
                        if(uChanCurrs(p,2)==curCurr && uChanCurrs(p,1)==curCh) % if valid entry exists
                            curRadii = groupedRadii{p};
                            nDates = popDates{p};
                            curDays = daysTrained(nDates);  
                            curDays(curRadii==0)=[];
                            curRadii(curRadii==0)=[];
                            
                            scatter(curDays,curRadii,[],currColors{cInd},'filled');

                            % Add trace to current average
                            traceNum = traceNum+1;
                            alignedDays = zeros(numel(curDays),1);
                            for d = 1:numel(curDays)
                                alignedDays(d) = find(daysTrainedOverall==curDays(d)); % Align days
                            end
                            mergingTraces(alignedDays,traceNum) = curRadii;
                        end
                    end
                end 
            end
        end
        % Average response for current and plot
        curAve = mean(mergingTraces,2,'omitnan');
        curDays = daysTrainedOverall;
        curDays = curDays(~isnan(curAve));
        curAve = curAve(~isnan(curAve));

    % % % % % % % %     f=fit(curDays',curAve,'exp1');
    % % % % % % % %     p1 = plot(f);
        % p1 = plot(curDays',curAve);
        p1 = plot(curDays(1:5)',curAve(1:5));

        tcolor = currColors{cInd};
        p1.Color = tcolor(1:3);
        p1.LineWidth=2;
    end
    title('Radius of Neural Activation over Time')
    ylabel('Enclosing Radius (\mum)')
    xlabel('Weeks of Training')
    
    lgd = legend('','Location', 'eastoutside');
    col_names = num2str(uCurrOverall);
    for j =1:length(col_names)
        plot([NaN NaN], [NaN NaN], 'Color', currColors{j}, 'DisplayName', strcat(col_names(j),' \muA'))
    end
    title(lgd,'Current')


    bob=1;
 %%   
%     
%     % ENCLOSING RADIUS vs. Threshold ----------------------------------------------
%     figure()
%     hold on
%     uCurrOverall = unique(cell2mat(uCurrAll));
%     uCurrOverall(uCurrOverall>6)=[];
%     uChanOverall = unique(cell2mat(uChanAll));
%     daysTrainedOverall = unique(ceil(cell2mat(daysTrainedAll')/dDiv));
%     for cInd = 2:numel(uCurrOverall) % process each current
% 
%         % Structures for averaging current specific response
%         mergingTraces = NaN(numel(daysTrainedOverall),numel(anID)*numel(uChanOverall));
%         traceNum=0;
% 
%         % Process each animal
%         for aInd = 1:numel(anID)
%             uChan = uChanAll{aInd};
%             uCurr = uCurrAll{aInd};
%             groupedRadii = enclosingRadiusAll{aInd};
%             daysTrained = ceil(daysTrainedAll{aInd}/dDiv);
%             uChanCurrs = uChanCurrsRef{aInd};
%             popDates = popDatesAll{aInd};
%             
%             
%             
%             daysThresholds = round(daysThresholdsAll{aInd});
%             uThresh = unique(reshape(daysThresholds,1,[]));
%             uThresh(uThresh==0)=[]; % remove thresholds at 0    
%             groupedCurrChansThresh = groupedCurrChansThreshAll{aInd};            
% 
%             curInd = find(uCurr==uCurrOverall(cInd));
%             if(curInd>0) % if current current is present for given animal, continue analysis
%                 for chInd = 1:numel(uChan) % process each channel
%                     if(aInd==5 && chInd==3)
%                         continue;
%                         % Skip ICMS 100 channel 6
%                     end
%                     neuCnt = groupedCurrChansThresh(:,chInd,curInd);
%                     curDays = uThresh(~isnan(neuCnt));
% 
% 
%           
%                     curCh = uChan(chInd);
%                     curCurr = uCurr(curInd);
%                     for p = 1:size(uChanCurrs,1)
%         
%                             curRadii = groupedRadii{p};
%                             nDates = popDates{p};
%                             curDays = daysTrained(nDates);  
%                             
%                             
%                             scatter(curDays,curRadii,[],currColors{cInd},'filled');
% 
%                             % Add trace to current average
%                             traceNum = traceNum+1;
%                             alignedDays = zeros(numel(curDays),1);
%                             for d = 1:numel(curDays)
%                                 alignedDays(d) = find(daysThresholdOverall==curDays(d)); % Align days
%                             end
%                             mergingTraces(alignedDays,traceNum) = curRadii;
%                     end
%                 end 
%             end
%         end
%         % Average response for current and plot
%         curAve = mean(mergingTraces,2,'omitnan');
%         curDays = daysTrainedOverall;
%         curDays = curDays(~isnan(curAve));
%         curAve = curAve(~isnan(curAve));
% 
%     % % % % % % % %     f=fit(curDays',curAve,'exp1');
%     % % % % % % % %     p1 = plot(f);
%         p1 = plot(curDays',curAve);
%         tcolor = currColors{cInd};
%         p1.Color = tcolor(1:3);
%         p1.LineWidth=2;
%     end
%     title('Radius of Neural Activation vs Threshold')
%     ylabel('Enclosing Radius (\mum)')
%     xlabel('Days of Training')
%     
%     lgd = legend('','Location', 'eastoutside');
%     col_names = num2str(uCurrOverall);
%     for j =1:length(col_names)
%         plot([NaN NaN], [NaN NaN], 'Color', currColors{j}, 'DisplayName', strcat(col_names(j),' \muA'))
%     end
%     title(lgd,'Current')
%     
    
    
    
    
    
    % ALL POPULATION DISTANCES ----------------------------------------------
%     for cInd = 1:numel(uCurrOverall) % process each current over days
%         figure()
%         hold on
%         g = [];
%         allDistances = [];
% 
%         % Process each animal
%         for aInd = 1:numel(anID)
%             uChan = uChanAll{aInd};
%             uCurr = uCurrAll{aInd};
%             daysTrained = daysTrainedAll{aInd};
%             uChanCurrs = uChanCurrsRef{aInd};
%             popDates = popDatesAll{aInd};
%             populationDistances = populationDistancesAll{aInd};
% 
%             curInd = find(uCurr==uCurrOverall(cInd));
%             if(curInd>0) % if current current is present for given animal, continue analysis
%                 for chInd = 1:numel(uChan) % process each channel
%                     curCurr = uCurr(curInd);
%                     for p = 1:size(uChanCurrs,1)
%                         if(uChanCurrs(p,2)==curCurr) % if valid entry exists
%                             curDistances = populationDistances{p};
%                             nDates = popDates{p};
%                             curDays = daysTrained(nDates);  
%                             for di = 1:numel(curDays)
%                                 d = curDays(di);
%                                 g = [g; d*ones(numel(curDistances{di}),1)];
%                                 allDistances = [allDistances; curDistances{di}];
%                             end
%                         end
%                     end
%                 end 
%             end
%         end
%         boxplot(allDistances,g)
%         aDays = unique(g);
%         title(strcat('Neural Activation Distance Vs Time : Curr',{' '},num2str(curCurr),' \muA'))
%         ylabel('Activation Distance (\mum)')
%         xlabel('Days of Training')
%         xticklabels(num2str(aDays))
%     end
    
    %While calculating distance also calculate an updated density
    %reflecting all 
    
    uCurrOverall = unique(cell2mat(uCurrAll));

    densityAll2 = cell(numel(anID),1);
    

    for cInd = 1:numel(uCurrOverall) % process each current over weeks
        figure()
        hold on
        g = [];
        allDistances = [];

        % Process each animal
        for aInd = 1:numel(anID)
            uChan = uChanAll{aInd};
            uCurr = uCurrAll{aInd};
            daysTrained = daysTrainedAll{aInd};
            uChanCurrs = uChanCurrsRef{aInd};
            popDates = popDatesAll{aInd};
            populationDistances = populationDistancesAll{aInd};

            densityCur = densityAll2{aInd};
            if(isempty(densityCur))
                % initialize variable if needed
                densityCur = NaN(numel(uCurrOverall),3,numel(daysTrained));
            end

            curInd = find(uCurr==uCurrOverall(cInd));
            if(curInd>0) % if current current is present for given animal, continue analysis
                for chInd = 1:numel(uChan) % process each channel
                    if(aInd==5 && chInd==3)
                        continue;
                        % Skip ICMS 100 channel 6
                    end
                    curCurr = uCurr(curInd);
                    for p = 1:size(uChanCurrs,1)
                        if(uChanCurrs(p,2)==curCurr) % if valid entry exists
                            curDistances = populationDistances{p};
                            nDates = popDates{p};
                            curDays = daysTrained(nDates);  
                            curDays = ceil(curDays/dDiv);
                            for di = 1:numel(curDays)
                                d = curDays(di);
                                g = [g; d*ones(numel(curDistances{di}),1)];
                                allDistances = [allDistances; curDistances{di}];
%                                 g = [g; d];
%                                 allDistances = [allDistances; mean(curDistances{di})];
                                mergedDensty = zeros(5,1);
                                for m = 1:5
                                    curcnt = allDistances>(m-1)*100 & allDistances<m*100;
                                    mergedDensity(m) = sum(curcnt)/(((4/3)*pi()*(100*m)^2)-((4/3)*pi()*(100*(m-1))^2));
                                end
                                
                                densityCur(cInd,chInd,di) = mean(mergedDensity);

                            end
                        end
                    end
                end 
            end
            
            densityAll2{aInd} = densityCur;
        end
%         boxplot(allDistances,g)
% %         violinplot(allDistances, g,'ViolinColor',[0 0 1]);
% 
%         aDays = unique(g);
%         title(strcat('Neural Activation Distance Vs Time : Curr',{' '},num2str(curCurr),' \muA'))
%         ylabel('Activation Distance (\mum)')
%         xlabel('Weeks of Training')
%         xticklabels(num2str(aDays))
        
        
        scatter(g,allDistances)

        aDays = unique(g);
        title(strcat('Neural Activation Distance Vs Time : Curr',{' '},num2str(curCurr),' \muA'))
        ylabel('Activation Distance (\mum)')
        xlabel('Weeks of Training')
        xticklabels(num2str(aDays))


        % [~,~,stats] = anova1(allDistances,g);
%         results = multcompare(stats);
        
    end
    
    



 uCurrSel = [3:6];
    for cInd = 1:numel(uCurrSel) % process each current over weeks
        figure()
        hold on
        g = [];
        allDistances = [];

        % Process each animal
        for aInd = 1:numel(anID)
            uChan = uChanAll{aInd};
            uCurr = uCurrAll{aInd};
            daysTrained = daysTrainedAll{aInd};
            uChanCurrs = uChanCurrsRef{aInd};
            popDates = popDatesAll{aInd};
            populationDistances = populationDistancesAll{aInd};

            curInd = find(uCurr==uCurrSel(cInd));
            if(curInd>0) % if current current is present for given animal, continue analysis
                for chInd = 1:numel(uChan) % process each channel
                    % if(aInd==5 && chInd==3)
                    %     continue;
                    %     % Skip ICMS 100 channel 6
                    % end
                    curCurr = uCurr(curInd);
                    for p = 1:size(uChanCurrs,1)
                        if(uChanCurrs(p,2)==curCurr) % if valid entry exists
                            if(uChanCurrs(p,1)==uChan(chInd)) % if valid entry exists
                                curDistances = populationDistances{p};
                                nDates = popDates{p};
                                curDays = daysTrained(nDates);  
                                curDays = ceil(curDays/dDiv);
                                for di = 1:numel(curDays)
                                    d = curDays(di);
                                    curDave =  prctile(curDistances{di},5);
                                    % curDave =  mean(curDistances{di});
                                    g = [g; d];
                                    allDistances = [allDistances; curDave];
    %                                 g = [g; d];
    %                                 allDistances = [allDistances; mean(curDistances{di})];
    
                                end
                            end
                        end
                    end
                end 
            end
        end
%         boxplot(allDistances,g)
% %         violinplot(allDistances, g,'ViolinColor',[0 0 1]);
% 
%         aDays = unique(g);
%         title(strcat('Neural Activation Distance Vs Time : Curr',{' '},num2str(curCurr),' \muA'))
%         ylabel('Activation Distance (\mum)')
%         xlabel('Weeks of Training')
%         xticklabels(num2str(aDays))
     
        g(allDistances==0)=[];
        allDistances(allDistances==0)=[];
        allDistances(g>4)=[];
        g(g>4)=[];
        g2 = double(g>1);
        scatter(g,allDistances)

        pAves = zeros(5,1);
        pSTDs = zeros(5,1);
        for h= 1:5
            pAves(h) = median(allDistances(g==(h-1)),'omitnan');
            pSTDs(h) = std(allDistances(g==(h-1)),'omitnan');
        end
        plot(0:4,pAves)

        % aDays = unique(g);
        title(strcat('Neural Activation Distance Vs Time : Curr',{' '},num2str(curCurr),' \muA'))
        % ylabel('Activation Distance (\mum)')
        % xlabel('Weeks of Training')
        % xticklabels(num2str(aDays))

        ranksum(allDistances(g<2),allDistances(g>1))
        % % % % % % [~,~,stats] = anovan(allDistances,g2);
%         results = multcompare(stats);
        
    end
    
%% at threshold analysis updated


    figure()
    hold on
    g = [];
    allDistances = [];

    % Process each animal
    for aInd = 1:numel(anID)
        uChan = uChanAll{aInd};
        % uCurr = uCurrAll{aInd};
        daysTrained = daysTrainedAll{aInd};
        weeksTrained = ceil(daysTrained./7);
        uChanCurrs = uChanCurrsRef{aInd};
        popDates = popDatesAll{aInd};
        populationDistances = populationDistancesAll{aInd};
        daysThresholds = daysThresholdsAll{aInd};

        for chInd = 1:numel(uChan) % process each channel
            for cDay = 1:numel(daysTrained) % process each day
                curThresh = round(daysThresholds(chInd,cDay));

                for p = 1:size(uChanCurrs,1)
                    if(uChanCurrs(p,2)==curThresh) % if valid entry exists
                        if(uChanCurrs(p,1)==uChan(chInd)) % if valid entry exists
                            curDistances = populationDistances{p};
                            nDates = popDates{p};
    
                            dInd = find(nDates==cDay); % see if target current was present for target day
                            if(dInd>0)
                                curDave =  prctile(curDistances{dInd},5);
                                % curDave =  mean(curDistances{di});
                                g = [g; weeksTrained(cDay)];
                                allDistances = [allDistances; curDave];
                            end
                        end
                    end
                end
            end
        end
    end

    g(allDistances==0)=[];
    allDistances(allDistances==0)=[];
    allDistances(g>4)=[];
    g(g>4)=[];
    g2 = double(g>1);
    scatter(g,allDistances,'filled')

    pAves = zeros(5,1);
    pSTDs = zeros(5,1);
    for h= 1:5
        pAves(h) = mean(allDistances(g==(h-1)),'omitnan');
        pSTDs(h) = std(allDistances(g==(h-1)),'omitnan');
    end
    plot(0:4,pAves)
    figure
    errorbar(0:4,pAves,pSTDs)

    aDays = unique(g);
    title(strcat('Neural Activation Distance At threshold : Curr',{' '},num2str(curCurr),' \muA'))
    ylabel('Activation Distance (\mum)')
    xlabel('Weeks of Training')
    xticklabels(num2str(aDays))


    % % % % [~,~,stats] = anovan(allDistances,g2);
%         results = multcompare(stats);
        
    


%% At threshold analysis
    for cInd = 1:numel(uCurrOverall) % process each current over Thresholds
        figure()
        hold on
        g = [];
        allDistances = [];

        % Process each animal
        for aInd = 1:numel(anID)
            uChan = uChanAll{aInd};
            uCurr = uCurrAll{aInd};
            uChanCurrs = uChanCurrsRef{aInd};
            popDates = popDatesAll{aInd};
            populationDistances = populationDistancesAll{aInd};
            daysThresholds = daysThresholdsAll{aInd};

            curInd = find(uCurr==uCurrOverall(cInd));
            if(curInd>0) % if current current is present for given animal, continue analysis
                for chInd = 1:numel(uChan) % process each channel
                    if(aInd==5 && chInd==3)
                        continue;
                        % Skip ICMS 100 channel 6
                    end
                    curCurr = uCurr(curInd);

                    for p = 1:size(uChanCurrs,1)
                        if(uChanCurrs(p,2)==curCurr) % if valid entry exists
                            curDistances = populationDistances{p};
                            nDates = popDates{p};
                            curThresh = round(daysThresholds(chInd,nDates));
                            for di = 1:numel(curThresh)
                                d = curThresh(di);
                                g = [g; d*ones(numel(curDistances{di}),1)];
                                allDistances = [allDistances; curDistances{di}];
                            end
                        end
                    end
                end 
            end
        end
        boxplot(allDistances,g)
%         violinplot(allDistances, g,'ViolinColor',[0 0 1]);
        aDays = unique(g);
        title(strcat('Neural Activation Distance Vs Threshold : Curr',{' '},num2str(curCurr),' \muA'))
        ylabel('Activation Distance (\mum)')
        xlabel('Detection Threshold (\muA)')
        xticklabels(num2str(aDays))
    end

%     
    
    


%    for cInd = 1:numel(uCurrOverall) % process each current over Thresholds
%         figure()
%         hold on
%         g = [];
%         allDistances = [];
% 
%         % Process each animal
%         for aInd = 1:numel(anID)
%             uChan = uChanAll{aInd};
%             uCurr = uCurrAll{aInd};
%             uChanCurrs = uChanCurrsRef{aInd};
%             popDates = popDatesAll{aInd};
%             populationDistances = populationDistancesAll{aInd};
%             daysThresholds = daysThresholdsAll{aInd};
% 
%             curInd = find(uCurr==uCurrOverall(cInd));
%             if(curInd>0) % if current current is present for given animal, continue analysis
%                 for chInd = 1:numel(uChan) % process each channel
%                     if(aInd==5 && chInd==3)
%                         continue;
%                         % Skip ICMS 100 channel 6
%                     end
%                     curCurr = uCurr(curInd);
% 
%                     for p = 1:size(uChanCurrs,1)
%                         if(uChanCurrs(p,2)==curCurr) % if valid entry exists
%                             curDistances = populationDistances{p};
%                             nDates = popDates{p};
%                             curThresh = round(daysThresholds(chInd,nDates));
%                             for di = 1:numel(curThresh)
%                                 d = curThresh(di);
%                                 g = [g; d*ones(numel(curDistances{di}),1)];
%                                 allDistances = [allDistances; curDistances{di}];
%                             end
%                         end
%                     end
%                 end 
%             end
%         end
%         boxplot(allDistances,g)
% %         violinplot(allDistances, g,'ViolinColor',[0 0 1]);
%         aDays = unique(g);
%         title(strcat('Neural Activation Distance Vs Threshold : Curr',{' '},num2str(curCurr),' \muA'))
%         ylabel('Activation Distance (\mum)')
%         xlabel('Detection Threshold (\muA)')
%         xticklabels(num2str(aDays))
%     end
% 
% 

    
    
    
    
    
    
    
%     ACTIVATION DENSITY vs Time ----------------------------------------------
    currColorsTurbo = jet(10);
    uCurrOverall = unique(cell2mat(uCurrAll));
    uCurrOverall(uCurrOverall>6)=[];
    uChanOverall = unique(cell2mat(uChanAll));
    daysTrainedOverall = unique(cell2mat(daysTrainedAll'));
    for cInd = 1:numel(uCurrOverall) % process each current
        
        for di = 1:10 % process each distance from electrode
            % Structures for averaging current specific response
            mergingTraces = NaN(numel(daysTrainedOverall),numel(anID)*numel(uChanOverall));
            traceNum=0;

            % Process each animal
            for aInd = 1:numel(anID)
                uCurr = uCurrAll{aInd};
                daysTrained = daysTrainedAll{aInd};
                allDensity = allDensityAll{aInd};
            
                curInd = find(uCurr==uCurrOverall(cInd));
                if(curInd>0) % if current current is present for given animal, continue analysis
                    curDat = allDensity(curInd,:,di);

                    toRmv = zeros(numel(curDat),1);
                    curDays = daysTrained;
                    curMeans = zeros(numel(curDays),1);
                    for cc = 1:numel(curDays)
                        if(isempty(curDat{cc}))
                            toRmv(cc)=1;
                        else
                            curMeans(cc) = mean(curDat{cc});
                        end
                    end
                    curMeans(toRmv==1)=[];
                    curDays(toRmv==1)=[];

                    % Only plot if there is data
                    if(~isempty(curMeans))
                        scatter(curDays,curMeans,[],currColorsTurbo(di,:),'filled');


                        % Add trace to current average
                        traceNum = traceNum+1;
                        alignedDays = zeros(numel(curDays),1);
                        for d = 1:numel(curDays)
                            alignedDays(d) = find(daysTrainedOverall==curDays(d)); % Align days
                        end
                        mergingTraces(alignedDays,traceNum) = curMeans;
                    end 
                end
            end
            
            % Average response for current and plot
            curAve = mean(mergingTraces,2,'omitnan');
            curDays = daysTrainedOverall;
            curDays = curDays(~isnan(curAve));
            curAve = curAve(~isnan(curAve));
% 
%             p1 = plot(curDays',curAve);
%             tcolor = currColorsTurbo(di,:);
%             p1.Color = tcolor(1:3);
%             p1.LineWidth=2;
        end

%         curCurr = uCurrOverall(cInd);
%         title(strcat('Overall Neural Activation Density Vs Time  :  Curr',num2str(curCurr)))
%         ylabel('Neural Activation Density (neurons \mum^3)')
%         xlabel('Days of Training') 
%         
%         lgd = legend('','Location', 'eastoutside');
%         col_names = {'0-100','100-200','200-300','300-400','400-500','500-600','600-700','700-800','800-900','900-1000'};
%         for j =1:length(col_names)
%             plot([NaN NaN], [NaN NaN], 'Color', currColorsTurbo(j,:), 'DisplayName', col_names{j})
%         end
%         title(lgd,'Distance Bin (\mum)')
    end



    
    % Continue temporal analysis but over weeks
    currColorsTurbo = jet(10);
    uCurrOverall = unique(cell2mat(uCurrAll));
    uCurrOverall(uCurrOverall>6)=[];
    uChanOverall = unique(cell2mat(uChanAll));
    daysTrainedOverall = unique(ceil(unique(cell2mat(daysTrainedAll'))/dDiv));
    
    weekDensityData = NaN(numel(uCurrOverall),numel(daysTrainedOverall),10);
    weekDensityDataSTD = NaN(numel(uCurrOverall),numel(daysTrainedOverall),10);
    
    
    for cInd = 1:numel(uCurrOverall) % process each current
        figure
        hold on
        currRawPnts = cell(numel(daysTrainedOverall),10);
        for di = 1:10 % process each distance from electrode
            % Structures for averaging current specific response
            mergingTraces = NaN(numel(daysTrainedOverall),numel(anID)*numel(uChanOverall));
            traceNum=0;

            % Process each animal
            for aInd = 1:numel(anID)
                uCurr = uCurrAll{aInd};
                daysTrained = daysTrainedAll{aInd};
                allDensity = allDensityAll{aInd};
            
                curInd = find(uCurr==uCurrOverall(cInd));
                if(curInd>0) % if current current is present for given animal, continue analysis
                    curDat = allDensity(curInd,:,di);

                    toRmv = zeros(numel(curDat),1);
                    curDays = ceil(daysTrained/dDiv);
                    curMeans = zeros(numel(curDays),1);
                    for cc = 1:numel(curDays)
                        if(isempty(curDat{cc}))
                            toRmv(cc)=1;
                        else
                            curMeans(cc) = mean(curDat{cc});
                        end
                    end
                    curMeans(toRmv==1)=[];
                    curDays(toRmv==1)=[];

                    % Only plot if there is data
                    if(~isempty(curMeans))
                        scatter(curDays,curMeans,[],currColorsTurbo(di,:),'filled');


                        % Add trace to current average
                        traceNum = traceNum+1;
                        alignedDays = zeros(numel(curDays),1);
                        for d = 1:numel(curDays)
                            alignedDays(d) = find(daysTrainedOverall==curDays(d)); % Align days
                            temp = currRawPnts{alignedDays(d),di};
                            temp = [temp,curMeans(d)];
                            currRawPnts{alignedDays(d),di} = temp;
                        end
                        mergingTraces(alignedDays,traceNum) = curMeans;
                    end 
                end
            end
            
            % Average response for current and plot
            curAve = mean(mergingTraces,2,'omitnan');
            curSTD = std(mergingTraces,[],2,'omitnan');
            curDays = daysTrainedOverall;
            curDays = curDays(~isnan(curAve));
            curAve = curAve(~isnan(curAve));
            curSTD = curSTD(~isnan(curSTD));

            
            weekInds = zeros(numel(curDays),1);
            for k = 1:numel(curDays)
                weekInds(k) = find(daysTrainedOverall==curDays(k));
            end
            
            weekDensityData(cInd,weekInds,di) = curAve;
            weekDensityDataSTD(cInd,weekInds,di) = curAve;

            p1 = plot(curDays',curAve);
            tcolor = currColorsTurbo(di,:);
            p1.Color = tcolor(1:3);
            p1.LineWidth=2;
        end

        curCurr = uCurrOverall(cInd);
        title(strcat('Activation Density Vs Time at',{' '},num2str(curCurr),'\muA'))
        ylabel('Neural Activation Density (neurons \mum^3)')
        xlabel('Weeks of Training') 

        lgd = legend('','Location', 'eastoutside');
        col_names = {'0-100','100-200','200-300','300-400','400-500','500-600','600-700','700-800','800-900','900-1000'};
        for j =1:length(col_names)
            plot([NaN NaN], [NaN NaN], 'Color', currColorsTurbo(j,:), 'DisplayName', col_names{j})
        end
        title(lgd,'Distance Bin (\mum)')

        grpA = [];
        grpB = [];
        for j = 1:5
            for k = 1:2
                grpA = [grpA, currRawPnts{k,j}];
            end
            for k = 3:5
                grpB = [grpB, currRawPnts{k,j}];
            end
        end
        if(~isempty(grpA))
            if(~isempty(grpB))
                ranksum(grpA,grpB)
            end
        end
    end

    
    
    
    %% calculate average density for each stimulation current over weeks of training
    weekDensityDataAll = cell(6,5,10);
    weekDensityDataAll500 = cell(6,5);
    weekDensityDataAll500Agr = cell(6,5);
    weekDensityDataAll1000 = cell(6,5);
    weekDensityData = NaN(6,5,10); % currents 1-6, weeks 0-4, distances 0 - 1000
    weekDensityDataSTD = NaN(6,5,10);
    weekDensityData500 = NaN(6,5); % currents 1-6, weeks 0-4, distances 0 - 1000
    weekDensityData500Agr = NaN(6,5); % currents 1-6, weeks 0-4, distances 0 - 1000
    weekDensityData1000 = NaN(6,5); % currents 1-6, weeks 0-4, distances 0 - 1000
    weekDensityDataSTD500 = NaN(6,5);
    weekDensityDataSTD500Agr = NaN(6,5);
    weekDensityDataSTD1000 = NaN(6,5);
    
    for cInd = 2:7 % process each current over weeks (should correspond to 1 through 6 uA)

        % Process each animal
        for aInd = 1:numel(anID)
            densityCur = densityAll{aInd};
            densityAgr = densityAllAgr{aInd};
            % densityCur(densityCur==0)=NaN;
            daysTrained = daysTrainedAll{aInd};
            weeksTrained = ceil(daysTrained./7);
            for wk = 0:4

                % Check for data during target week
                if(sum(weeksTrained==wk)>0)
                    temp = squeeze(densityCur(cInd,:,weeksTrained==wk,:));

                    td = ndims(temp);
                    if(td>2)
                        for k = 1:size(temp,2)
                            for di = 1:10
                                weekDensityDataAll{cInd-1,wk+1,di} = [weekDensityDataAll{cInd-1,wk+1,di}; squeeze(temp(:,k,di))'];
                            end
                            weekDensityDataAll500{cInd-1,wk+1} = [weekDensityDataAll500{cInd-1,wk+1}; mean(squeeze(temp(:,k,1:5)),2)'];
                            % squeeze(temp(:,k,1:5))
                            weekDensityDataAll1000{cInd-1,wk+1} = [weekDensityDataAll1000{cInd-1,wk+1}; squeeze(temp(:,k,:))'];
                        end
                    else
                        for di = 1:10
                            weekDensityDataAll{cInd-1,wk+1,di} = [weekDensityDataAll{cInd-1,wk+1,di}; temp(:,di)'];
                        end
                        weekDensityDataAll500{cInd-1,wk+1} = [weekDensityDataAll500{cInd-1,wk+1}; mean(temp(:,1:5),2)'];
                        weekDensityDataAll1000{cInd-1,wk+1} = [weekDensityDataAll1000{cInd-1,wk+1}; temp(:,:)'];
                    end

                    tempAgr = squeeze(densityAgr(cInd,:,weeksTrained==wk));
                    if(size(tempAgr,1)>1)
                        weekDensityDataAll500Agr{cInd-1,wk+1} = [weekDensityDataAll500Agr{cInd-1,wk+1}; tempAgr'];
                    else
                        weekDensityDataAll500Agr{cInd-1,wk+1} = [weekDensityDataAll500Agr{cInd-1,wk+1}; tempAgr];
                    end

                end
            end
        end
    end
   
    for cInd = 1:6
        for wk = 1:5
            for di = 1:10
                if(~isempty(weekDensityDataAll{cInd,wk,di}))
                    weekDensityData(cInd,wk,di) = mean(weekDensityDataAll{cInd,wk,di},'all','omitnan');
                    weekDensityDataSTD(cInd,wk,di) = std(weekDensityDataAll{cInd,wk,di},[],'all','omitnan');
                end
            end
            if(~isempty(weekDensityDataAll500{cInd,wk}))
                weekDensityData500(cInd,wk) = mean(weekDensityDataAll500{cInd,wk},'all','omitnan');
                weekDensityDataSTD500(cInd,wk) = std(weekDensityDataAll500{cInd,wk},[],'all','omitnan');


                weekDensityData1000(cInd,wk) = mean(weekDensityDataAll1000{cInd,wk},'all','omitnan');
                weekDensityDataSTD1000(cInd,wk) = std(weekDensityDataAll1000{cInd,wk},[],'all','omitnan');
                weekDensityData500Agr(cInd,wk) = mean(weekDensityDataAll500Agr{cInd,wk},'all','omitnan');
                weekDensityDataSTD500Agr(cInd,wk) = std(weekDensityDataAll500Agr{cInd,wk},[],'all','omitnan');
            end
        end
    end

    %% Plot updated density trends
    for cInd = 4:5
        figure()
        hold on
        g = [];
        allDense = [];
        for wk=0:4
           temp = weekDensityDataAll500{cInd,wk+1};
           temp(isnan(temp))=[];
           allDense = [allDense, reshape(temp,1,[])];
           g = [g, ones(1,numel(temp))*wk];
        end
        
        scatter(g,allDense,[],currColorsTurbo(cInd,:),'filled');
        plot(0:4,weekDensityData500(cInd,:))
        % errorbar(0:4,weekDensityData500(cInd,:),weekDensityDataSTD500(cInd,:));
        title(strcat('Activation Density Vs Time at',{' '},num2str(cInd),'\muA'))
        ylabel('Neural Activation Density (neurons mm^3)')
        xlabel('Weeks of Training') 
        disp(strcat("Stat for Density at ",{' '},num2str(cInd)))
        ranksum(allDense(g<2),allDense(g>1))
    
    end
    



    %% Plot updated density trends
    for cInd = 4:5
        figure()
        hold on
        g = [];
        allDense = [];
        for wk=0:4
           temp = weekDensityDataAll500Agr{cInd,wk+1};
           allDense = [allDense, reshape(temp,1,[])];
           g = [g, ones(1,numel(temp))*wk];
        end

        % scatter(g,allDense,[],currColorsTurbo(cInd,:),'filled');
        errorbar(0:4,weekDensityData500Agr(cInd,:),weekDensityDataSTD500Agr(cInd,:));
        title(strcat('Agrigated Activation Density Vs Time at',{' '},num2str(cInd),'\muA'))
        ylabel('Neural Activation Density (neurons mm^3)')
        xlabel('Weeks of Training') 
        disp(strcat("Stat for Density at ",{' '},num2str(cInd)))
        ranksum(allDense(g<2),allDense(g>1))

    end

%%
    for cInd = 3:6
        figure()
        hold on
        % g = [];
        % allDense = [];
        % for wk=0:4
        %    temp = weekDensityDataAll1000{cInd,wk+1};
        %    allDense = [allDense, reshape(temp,1,[])];
        %    g = [g, ones(1,numel(temp))*wk];
        % end

        for di = 1:10
            curMeans = weekDensityData(cInd,:,di);
            curSTDs = weekDensityDataSTD(cInd,:,di);
    
                p1 = errorbar(0:4,curMeans,curSTDs);
                tcolor = currColorsTurbo(di,:);
                p1.Color = tcolor(1:3);
                p1.LineWidth=2;
            
            % scatter(g,allDense,[],currColorsTurbo(cInd,:),'filled');
            
        end
        title(strcat('Activation Density Sweep Vs Time at',{' '},num2str(cInd),'\muA'))
        ylabel('Neural Activation Density (neurons mm^3)')
        xlabel('Weeks of Training') 

        lgd = legend('','Location', 'eastoutside');
        col_names = {'0-100','100-200','200-300','300-400','400-500','500-600','600-700','700-800','800-900','900-1000'};
        for j =1:length(col_names)
            plot([NaN NaN], [NaN NaN], 'Color', currColorsTurbo(j,:), 'DisplayName', col_names{j})
        end
        title(lgd,'Distance Bin (\mum)')
    end
    
 
    %% Ring Plots for each current level over the course of weeks
    cmap = jet(256);
    for c = 4:5
        maxDen3DCOM = max(weekDensityData(c,:,:),[],'all');
        for t = 1:5 % For each week of training

            threshAves = squeeze(weekDensityData(c,t,:));
            threshAves(threshAves==0) = 1e-15;
            if(~isnan(threshAves(1)))
                figName = strcat('Current: ',{' '},num2str(c),' uA - Week:',{' '},num2str(daysTrainedOverall(t)));
                figure('Name',figName{1},'NumberTitle','off');
                hold on
                for k = 10:-1:1
                    x1 = 0.5*abs(k-11);
                    pos = [x1 x1 k k]; 
                    if(threshAves(k)==0)
                        cmapS = [1,1,1];
                    else
                        cmapS = cmap(ceil(256*threshAves(k)/maxDen3DCOM),:);
                    end
                    rectangle('Position',pos,'Curvature',[1 1],'FaceColor',cmapS)
                    axis equal 
                end
                caxis([0 maxDen3DCOM])
                colormap(jet)
                cB = colorbar;
                set(gca,'fontsize', 28)
                axis off

            end
        end
    end
    
    
    
    
%% Ring PLots and 2D plot of activation at fixed currents



    cmap = jet(256);
    mDen = NaN(6,1);
    sDen = NaN(6,1);
    datAll = [];
    g = [];
    maxDen3DCOM = max(weekDensityData(3:6,:,:),[],'all');

    for c = 3:6
        threshAves = squeeze(weekDensityData(c,:,:));
       
        curDvals = reshape(threshAves(:,1:5),1,[]);
        mDen(c) = mean(curDvals);
        sDen(c) = std(curDvals);
        datAll = [datAll, curDvals];
        g = [g; ones(numel(curDvals),1)*c];

        threshAves(threshAves==0) = 1e-15;
        if(~isnan(threshAves(1)))
            figName = strcat('Current: ',{' '},num2str(uCurrOverall(c)),' uA');
            figure('Name',figName{1},'NumberTitle','off');
            hold on
            for k = 10:-1:1
                x1 = 0.5*abs(k-11);
                pos = [x1 x1 k k]; 
                if(threshAves(k)==0)
                    cmapS = [1,1,1];
                else
                    cmapS = cmap(ceil(256*threshAves(k)/maxDen3DCOM),:);
                end
                rectangle('Position',pos,'Curvature',[1 1],'FaceColor',cmapS)
%                     text(5.5,x1+0.25,ringText{k},'FontSize',8)
                axis equal 
            end
            caxis([0 maxDen3DCOM])
            colormap(jet)
            cB = colorbar;
            set(gca,'fontsize', 28)
            axis off

        end
        
    end

    figure()
    hold on
    scatter(g, datAll,'blue','filled')
    plot(3:6,mDen(3:6))
    xlim([2.5 6.5])


    [~,~,stats] = anovan(datAll,g);
    
        results = multcompare(stats);




%%
    
    
    
    
    




    % % % % % % 
    % % % % % % 
    % % % % % % 
    % % % % % % % Continue temporal analysis but over weeks
    % % % % % % uCurrOverall = unique(cell2mat(uCurrAll));
    % % % % % % uCurrOverall(uCurrOverall>6)=[];
    % % % % % % uChanOverall = unique(cell2mat(uChanAll));
    % % % % % % daysTrainedOverall = unique(ceil(unique(cell2mat(daysTrainedAll'))/dDiv));
    % % % % % % 
    % % % % % % weekDensityData = NaN(numel(uCurrOverall),numel(daysTrainedOverall),10);
    % % % % % % weekDensityDataSTD = NaN(numel(uCurrOverall),numel(daysTrainedOverall),10);
    % % % % % % 
    % % % % % % 
    % % % % % % 
    % % % % % % for cInd = 1:numel(uCurrOverall) % process each current
    % % % % % %     for di = 1:10 % process each distance from electrode
    % % % % % %         % Structures for averaging current specific response
    % % % % % %         mergingTraces = NaN(numel(daysTrainedOverall),numel(anID)*numel(uChanOverall));
    % % % % % %         mergingCells = cell(numel(daysTrainedOverall),numel(anID)*numel(uChanOverall));
    % % % % % %         traceNum=0;
    % % % % % % 
    % % % % % %         % Process each animal
    % % % % % %         for aInd = 1:numel(anID)
    % % % % % %             uCurr = uCurrAll{aInd};
    % % % % % %             uChan = uChanAll{aInd};
    % % % % % %             daysTrained = daysTrainedAll{aInd};
    % % % % % %             uChanCurrs = uChanCurrsRef{aInd};
    % % % % % %             popDates = popDatesAll{aInd};
    % % % % % %             populationDistances = populationDistancesAll{aInd};
    % % % % % % 
    % % % % % % 
    % % % % % % 
    % % % % % % 
    % % % % % %             curInd = find(uCurr==uCurrOverall(cInd));
    % % % % % %             if(curInd>0) % if current current is present for given animal, continue analysis
    % % % % % %                 % for chInd = 1:numel(uChan) % process each channel
    % % % % % %                 %     % if(aInd==5 && chInd==3)
    % % % % % %                 %     %     continue;
    % % % % % %                 %     %     % Skip ICMS 100 channel 6
    % % % % % %                 %     % end
    % % % % % %                 curCurr = uCurr(curInd);
    % % % % % %                 for p = 1:size(uChanCurrs,1)
    % % % % % %                     if(uChanCurrs(p,2)==curCurr) % if valid entry exists
    % % % % % %                         curDistances = populationDistances{p};
    % % % % % %                         nDates = popDates{p};
    % % % % % %                         curDays = daysTrained(nDates);  
    % % % % % %                         curDays = ceil(curDays/dDiv);
    % % % % % %                         for da = 1:numel(curDays)
    % % % % % % 
    % % % % % %                             % Add trace to current average
    % % % % % %                             traceNum = traceNum+1;
    % % % % % %                             alignedD = find(daysTrainedOverall==curDays(da)); % Align days
    % % % % % % 
    % % % % % %                             % Count number of neurons in target
    % % % % % %                             % range
    % % % % % %                             nInRangeA = (100*(di-1))<curDistances{da};
    % % % % % %                             nInRangeB = (100*di)>=curDistances{da};
    % % % % % %                             nInRange = nInRangeA & nInRangeB;
    % % % % % % 
    % % % % % %                             mergingTraces(alignedD,traceNum) = sum(nInRange) / (((4/3) * pi * (di*100).^3)-((4/3) * pi * ((di-1)*100).^3));
    % % % % % %                             mergingCells{alignedD,traceNum} = curDistances{da};
    % % % % % %                         end
    % % % % % %                     end
    % % % % % %                     % end
    % % % % % %                 end 
    % % % % % %             end
    % % % % % %         end
    % % % % % % 
    % % % % % %         % Average response for current and plot
    % % % % % %         curAve = mean(mergingTraces,2,'omitnan');
    % % % % % %         curSTD = std(mergingTraces,[],2,'omitnan');
    % % % % % %         curDays = daysTrainedOverall;
    % % % % % %         curDays = curDays(~isnan(curAve));
    % % % % % %         curAve = curAve(~isnan(curAve));
    % % % % % %         curSTD = curSTD(~isnan(curSTD));
    % % % % % % 
    % % % % % % 
    % % % % % %         weekInds = zeros(numel(curDays),1);
    % % % % % %         for k = 1:numel(curDays)
    % % % % % %             weekInds(k) = find(daysTrainedOverall==curDays(k));
    % % % % % %         end
    % % % % % % 
    % % % % % %         weekDensityData(cInd,weekInds,di) = curAve;
    % % % % % %         weekDensityDataSTD(cInd,weekInds,di) = curSTD;
    % % % % % %     end
    % % % % % % end
    % % % % % % 





    
% 
% 
%     %% Plot activation density with respect to weeks as traces rather than
%     % distance
%     currColorsT  = jet(numel(daysTrainedOverall));
%     for c = 1:numel(uCurrOverall)
%         figure
%         hold on
%         for t = 1:numel(daysTrainedOverall)
%             threshAves = squeeze(weekDensityData(c,t,:));
%             p1 = plot(1:10,threshAves,'-o');
%             tcolor = currColorsT(t,:);
%             p1.Color = tcolor(1:3);
%             p1.LineWidth=2;
%         end
%         xticks(1:10)
%         xticklabels({'0-100','100-200','200-300','300-400','400-500','500-600','600-700','700-800','800-900','900-1000'})
%         ylabel('Neural Activation Density (neurons \mum^3)')
%         xlabel('Distance from Electrode (\mum)')
%         xtickangle(45)
%         title(strcat('Activation Density at',{' '},num2str(uCurrOverall(c)),'\muA'))
% 
% 
%         lgd = legend('','Location', 'eastoutside');
%         col_names = num2str(daysTrainedOverall');
%         for j =1:length(col_names)
%             plot([NaN NaN], [NaN NaN], 'Color', currColorsT(j,:), 'DisplayName', col_names(j))
%         end
%         title(lgd,'Weeks of Training')
%     end
% 
% 
% 
% 
% 
%     %% Ring Plots for each current level over the course of weeks
% %     ringText = {'0-100','100-200','200-300','300-400','400-500','500-600','600-700','700-800','800-900','900-1000'};
%     cmap = jet(256);
%     for c = 1:numel(uCurrOverall)
%         maxDen3DCOM = max(weekDensityData(c,1:5,:),[],'all');
%         for t = 1:5%numel(daysTrainedOverall) % For each week of training
% 
%             threshAves = squeeze(weekDensityData(c,t,:));
%             threshAves(threshAves==0) = 1e-15;
%             if(~isnan(threshAves(1)))
%                 figName = strcat('Current: ',{' '},num2str(uCurrOverall(c)),' uA - Week:',{' '},num2str(daysTrainedOverall(t)));
%                 figure('Name',figName{1},'NumberTitle','off');
%                 hold on
%                 for k = 10:-1:1
%                     x1 = 0.5*abs(k-11);
%                     pos = [x1 x1 k k]; 
%                     if(threshAves(k)==0)
%                         cmapS = [1,1,1];
%                     else
%                         cmapS = cmap(ceil(256*threshAves(k)/maxDen3DCOM),:);
%                     end
%                     rectangle('Position',pos,'Curvature',[1 1],'FaceColor',cmapS)
% %                     text(5.5,x1+0.25,ringText{k},'FontSize',8)
%                     axis equal 
%                 end
%                 caxis([0 maxDen3DCOM])
%                 colormap(jet)
%                 cB = colorbar;
%                 set(gca,'fontsize', 28)
%                 axis off
% 
%             end
%         end
%     end
% 
% 
% 
    
% 
%     cmap = jet(256);
%     mDen = NaN(numel(uCurrOverall),1);
%     sDen = NaN(numel(uCurrOverall),1);
%     datAll = [];
%     g = [];
%     for c = 1:numel(uCurrOverall)
%         maxDen3DCOM = max(weekDensityData(:,2,:),[],'all');
% 
%         threshAves = squeeze(weekDensityData(c,2,:));
%         mDen(c) = mean(threshAves(1:5));
%         sDen(c) = std(threshAves(1:5));
%         datAll = [datAll; threshAves(1:5)];
%         g = [g; ones(5,1)*c];
% 
%         threshAves(threshAves==0) = 1e-15;
%         if(~isnan(threshAves(1)))
%             figName = strcat('Current: ',{' '},num2str(uCurrOverall(c)),' uA');
%             figure('Name',figName{1},'NumberTitle','off');
%             hold on
%             for k = 10:-1:1
%                 x1 = 0.5*abs(k-11);
%                 pos = [x1 x1 k k]; 
%                 if(threshAves(k)==0)
%                     cmapS = [1,1,1];
%                 else
%                     cmapS = cmap(ceil(256*threshAves(k)/maxDen3DCOM),:);
%                 end
%                 rectangle('Position',pos,'Curvature',[1 1],'FaceColor',cmapS)
% %                     text(5.5,x1+0.25,ringText{k},'FontSize',8)
%                 axis equal 
%             end
%             caxis([0 maxDen3DCOM])
%             colormap(jet)
%             cB = colorbar;
%             set(gca,'fontsize', 28)
%             axis off
% 
%         end
% 
%     end
% 
%     figure
%     errorbar(uCurrOverall(3:7),mDen(3:7),sDen(3:7))
%     xlim([1.5 6.5])
% 
    
%     %% Plot the density of activation for first 500 um and radius of
%     % activation versus time
%     for cInd = 1:numel(uCurrOverall) % for each current
%         if(uCurrOverall(cInd)<3)
%             continue
%         end
%         figName = strcat('Desnsity Vs Enclosing Radius - Current: ',{' '},num2str(uCurrOverall(cInd)));
%         figure('Name',figName{1},'NumberTitle','off');
% 
%         densityAves = mean(squeeze(weekDensityData(cInd,:,1:3)),2);
%         densityAveSTD = mean(squeeze(weekDensityDataSTD(cInd,:,1:3)),2);
%         plot(daysTrainedOverall(1:5),densityAves(1:5));
% 
% %         errorbar(daysTrainedOverall(1:5),densityAves(1:5),densityAveSTD(1:5));
%         hold on
%         yyaxis right
% 
% 
%         % Structures for averaging current specific response
%         mergingTraces = NaN(numel(daysTrainedOverall),numel(anID)*numel(uChanOverall));
%         traceNum=0;
% 
%         % Process each animal
%         for aInd = 1:numel(anID)
%             uChan = uChanAll{aInd};
%             uCurr = uCurrAll{aInd};
%             groupedRadii = enclosingRadiusAll{aInd};
%             daysTrained = ceil(daysTrainedAll{aInd}/dDiv);
%             uChanCurrs = uChanCurrsRef{aInd};
%             popDates = popDatesAll{aInd};
% 
%             curInd = find(uCurr==uCurrOverall(cInd));
%             if(curInd>0) % if current current is present for given animal, continue analysis
%                 for chInd = 1:numel(uChan) % process each channel
%                     if(aInd==5 && chInd==3)
%                         continue;
%                         % Skip ICMS 100 channel 6
%                     end
%                     curCh = uChan(chInd);
%                     curCurr = uCurr(curInd);
%                     for p = 1:size(uChanCurrs,1)
%                         if(uChanCurrs(p,2)==curCurr && uChanCurrs(p,1)==curCh) % if valid entry exists
%                             curRadii = groupedRadii{p};
%                             nDates = popDates{p};
%                             curDays = daysTrained(nDates);  
%                             curDays(curRadii==0)=[];
%                             curRadii(curRadii==0)=[];
% 
%                             scatter(curDays,curRadii,[],'filled');
% 
%                             % Add trace to current average
%                             traceNum = traceNum+1;
%                             alignedDays = zeros(numel(curDays),1);
%                             for d = 1:numel(curDays)
%                                 alignedDays(d) = find(daysTrainedOverall==curDays(d)); % Align days
%                             end
%                             mergingTraces(alignedDays,traceNum) = curRadii;
%                         end
%                     end
%                 end 
%             end
%         end
%         % Average response for current and plot
%         curAve = mean(mergingTraces,2,'omitnan');
%         curDays = daysTrainedOverall;
%         % curDays = curDays(~isnan(curAve));
%         % curAve = curAve(~isnan(curAve));
% 
%         curDays = curDays(1:5);
%         curAve = curAve(1:5);
% 
%         plot(curDays',curAve);
%         % p1.LineWidth=2;
%         xlim([-0.5, 4.5])
% 
% 
%         % ranksum(reshape(mergingTraces(1:2,:),1,[]),reshape(mergingTraces(3:5,:),1,[]))
% 
% 
%     end
% 
% 
    
    
  %% Radius of activation across currents
 uCurrSel = 3:6;
    g = [];
    allDistances = [];
    allDistRaw = [];
    allG = [];
    for cInd = 1:numel(uCurrSel) % process each current over weeks


        % Process each animal
        for aInd = 1:numel(anID)
            uChan = uChanAll{aInd};
            uCurr = uCurrAll{aInd};
            daysTrained = daysTrainedAll{aInd};
            uChanCurrs = uChanCurrsRef{aInd};
            popDates = popDatesAll{aInd};
            populationDistances = populationDistancesAll{aInd};

            curInd = find(uCurr==uCurrSel(cInd));
            if(curInd>0) % if current current is present for given animal, continue analysis
                for chInd = 1:numel(uChan) % process each channel
                     
                    curCurr = uCurr(curInd);
                    for p = 1:size(uChanCurrs,1)
                        if(uChanCurrs(p,2)==curCurr) % if valid entry exists
                            if(uChanCurrs(p,1)==uChan(chInd)) % if valid entry exists
                                curDistances = populationDistances{p};
                                nDates = popDates{p};
                                curDays = daysTrained(nDates);  
                                curDays = ceil(curDays/dDiv);
                                for di = 1:numel(curDays)
                                    if(curDays(di)==1)
                                        % curDave =  prctile(curDistances{di},80);
                                        curDave = mean(curDistances{di});
                                        g = [g; curCurr];
                                        allDistances = [allDistances; curDave];
                                        allDistRaw = [allDistRaw;curDistances{di}];
                                        allG = [allG ; ones(numel(curDistances{di}),1)*curCurr];
                                    end
                                end
                            end
                        end
                    end
                end 
            end
        end
    end
   
    figure
    scatter(g,allDistances)

    pAves = zeros(4,1);
    pSTDs = zeros(4,1);
    for h= 3:6
        pAves(h-2) = mean(allDistances(g==h),'omitnan');
        pSTDs(h-2) = std(allDistances(g==h),'omitnan');
    end
    figure()
    errorbar(3:6,pAves,pSTDs)

    % aDays = unique(g);
    title(strcat('Neural Activation Distance Vs Time : Curr',{' '},num2str(curCurr),' \muA'))
    % ylabel('Activation Distance (\mum)')
    % xlabel('Weeks of Training')
    % xticklabels(num2str(aDays))


    [~,~,stats] = anovan(allDistances,g);
    
        results = multcompare(stats);
        
    
      
    
    
    
    
    
    %% ENCLOSING RADIUS vs. Time, seperated ----------------------------------------------
    figure()
    hold on
    uCurrOverall = unique(cell2mat(uCurrAll));
    uCurrOverall(uCurrOverall>6)=[];
    uCurrOverall(uCurrOverall<2)=[];
    uChanOverall = unique(cell2mat(uChanAll));
    daysTrainedOverall = unique(ceil(cell2mat(daysTrainedAll')/dDiv));
    
    for cInd = 2:numel(uCurrOverall) % process each current

        % Structures for averaging current specific response
        mergingTraces = NaN(numel(daysTrainedOverall),numel(anID)*numel(uChanOverall));
        traceNum=0;

        % Process each animal
        for aInd = 1:numel(anID)
            uChan = uChanAll{aInd};
            uCurr = uCurrAll{aInd};
            groupedRadii = enclosingRadiusAll{aInd};
            daysTrained = ceil(daysTrainedAll{aInd}/dDiv);
            uChanCurrs = uChanCurrsRef{aInd};
            popDates = popDatesAll{aInd};

            curInd = find(uCurr==uCurrOverall(cInd));
            if(curInd>0) % if current current is present for given animal, continue analysis
                for chInd = 1:numel(uChan) % process each channel
                    if(aInd==5 && chInd==3)
                        continue;
                        % Skip ICMS 100 channel 6
                    end
                    curCh = uChan(chInd);
                    curCurr = uCurr(curInd);
                    for p = 1:size(uChanCurrs,1)
                        if(uChanCurrs(p,2)==curCurr && uChanCurrs(p,1)==curCh) % if valid entry exists
                            curRadii = groupedRadii{p};
                            nDates = popDates{p};
                            curDays = daysTrained(nDates);  
                            curDays(curRadii==0)=[];
                            curRadii(curRadii==0)=[];
                            
                            scatter(curDays,curRadii,[],currColors{cInd},'filled');

                            % Add trace to current average
                            traceNum = traceNum+1;
                            alignedDays = zeros(numel(curDays),1);
                            for d = 1:numel(curDays)
                                alignedDays(d) = find(daysTrainedOverall==curDays(d)); % Align days
                            end
                            mergingTraces(alignedDays,traceNum) = curRadii;
                        end
                    end
                end 
            end
        end
        % Average response for current and plot
        curAve = mean(mergingTraces,2,'omitnan');
        curDays = daysTrainedOverall;
        curDays = curDays(~isnan(curAve));
        curAve = curAve(~isnan(curAve));

        p1 = plot(curDays',curAve);
        tcolor = currColors{cInd};
        p1.Color = tcolor(1:3);
        p1.LineWidth=2;
    end
    title('Radius of Neural Activation over Time')
    ylabel('Enclosing Radius (\mum)')
    xlabel('Weeks of Training')
    
    lgd = legend('','Location', 'eastoutside');
    col_names = num2str(uCurrOverall);
    for j =1:length(col_names)
        plot([NaN NaN], [NaN NaN], 'Color', currColors{j}, 'DisplayName', strcat(col_names(j),' \muA'))
    end
    title(lgd,'Current')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %% ACTIVATED DENSITY vs Threshold  
    currColorsTurbo = jet(10);
    uCurrOverall = unique(cell2mat(uCurrAll));
    uCurrOverall(uCurrOverall>6)=[];
    uChanOverall = unique(cell2mat(uChanAll));
    daysThresholdOverall = unique(unique(round(cell2mat(daysThresholdsAll'))));
    daysThresholdOverall(daysThresholdOverall==0)=[];
    daysThresholdOverall(daysThresholdOverall>7)=[];
    
    
    
    threshDensityData = NaN(numel(uCurrOverall),numel(daysThresholdOverall),10);
    
    for cInd = 1:numel(uCurrOverall) % process each current
%         figure
%         hold on
        
        for di = 1:10 % process each distance from electrode
            % Structures for averaging current specific response
            mergingTraces = NaN(numel(daysThresholdOverall),numel(anID)*3*numel(uChanOverall));
            traceNum=0;

            % Process each animal
            for aInd = 1:numel(anID)
                uCurr = uCurrAll{aInd};
                uChan = uChanAll{aInd};
                daysThresholds = round(daysThresholdsAll{aInd});
                allDensityThresh = allDensityThreshAll{aInd};
                
                uThresholds = unique(unique(daysThresholds));
                uThresholds(uThresholds==0)=[];

                curInd = find(uCurr==uCurrOverall(cInd));
                if(curInd>0) % if current current is present for given animal, continue analysis
                    for chInd = 1:numel(uChan)% Iterate through each channel
                        if(aInd==5 && chInd==3)
                            continue;
                            % Skip ICMS 100 channel 6
                        end
                        curDat = squeeze(allDensityThresh(chInd,curInd,:,di));
                        curThresh = uThresholds;

                        % Only plot if there is data
                        if(~isempty(curDat))


                            
                            scatter(curThresh,curDat,[],currColorsTurbo(di,:),'filled');

                            % Add trace to current average
                            for d = 1:numel(curThresh)
                                % employ lazy abuse of memory to maintain data
                                % accuracy and avoid convoluted averaging
                                % algorithm.....  heh it works
                                alignedThresh = find(daysThresholdOverall==curThresh(d)); % Align thresholds
                                traceNum = traceNum+1;
                                mergingTraces(alignedThresh,traceNum) = curDat(d);
                            end
                        end 
                    end
                end
            end
            
            % Average response for current and plot
            curAve = mean(mergingTraces,2,'omitnan');
            curThreshs = daysThresholdOverall;
            curThreshs = curThreshs(~isnan(curAve));
            curAve = curAve(~isnan(curAve));
            
            % Add means to matrix for alternative plotting
            threshInds = zeros(numel(curThreshs),1);
            for k = 1:numel(curThreshs)
                threshInds(k) = find(daysThresholdOverall==curThreshs(k));
            end
            
            threshDensityData(cInd,threshInds,di) = curAve;

%             p1 = plot(curThreshs',curAve);
%             tcolor = currColorsTurbo(di,:);
%             p1.Color = tcolor(1:3);
%             p1.LineWidth=2;
        end

%         curCurr = uCurrOverall(cInd);
%         title(strcat('Activation Density Vs Threshold at',{' '},num2str(curCurr),'\muA'))
%         ylabel('Neural Activation Density (neurons \mum^3)')
%         xlabel('Detection Thresholds (\muA)')
%         
%         lgd = legend('','Location', 'eastoutside');
%         col_names = {'0-100','100-200','200-300','300-400','400-500','500-600','600-700','700-800','800-900','900-1000'};
%         for j =1:length(col_names)
%             plot([NaN NaN], [NaN NaN], 'Color', currColorsTurbo(j,:), 'DisplayName', col_names{j})
%         end
%         title(lgd,'Distance Bin (\mum)')
    end

    



%% Density at threshold - single plot
 figure()
    hold on
    g = [];
    threshDensities = [];
    threshDensitiesGrpA = [];
    threshDensitiesGrpB = [];
    weekThreshDense = cell(10,1);
    

    % Process each animal
    for aInd = 1:numel(anID)
        uChan = uChanAll{aInd};
        uCurr = uCurrAll{aInd};
        daysTrained = daysTrainedAll{aInd};
        weeksTrained = ceil(daysTrained./7);
        daysThresholds = daysThresholdsAll{aInd};
        curDensity = densityAll{aInd};

        for chInd = 1:numel(uChan) % process each channel
            longDens = NaN(numel(daysTrained),1);
            for cDay = 1:numel(daysTrained) % process each day
                curThresh = round(daysThresholds(chInd,cDay));
                % curInd = find(uCurr==curThresh);
                % if(~isempty(curInd))
                if(curThresh>0)
                    curDat = mean(curDensity(curThresh,chInd,cDay,1:5),'omitnan');
                    threshDensities = [threshDensities, curDat];
                    g = [g, weeksTrained(cDay)];
                    longDens(cDay) = curDat;

                    wTemp = weekThreshDense{weeksTrained(cDay)+1};
                    wTemp = [wTemp; curDensity(curThresh,chInd,cDay,:)];
                    weekThreshDense{weeksTrained(cDay)+1} = wTemp;
                end
            end

            threshDensitiesGrpA = [threshDensitiesGrpA, mean(longDens(weeksTrained<2),'omitnan')];
            longB = weeksTrained>1 & weeksTrained<5;
            threshDensitiesGrpB = [threshDensitiesGrpB, mean(longDens(3:5),'omitnan')];
        end
    end

    g(threshDensities==0)=[];
    threshDensities(threshDensities==0)=[];
    threshDensities(g>5)=[];
    g(g>5)=[];
    % g2 = double(g>1);
    % scatter(g,threshDensities,'filled')

    pAves = zeros(5,1);
    pSTDs = zeros(5,1);
    for h= 1:5
        datP = threshDensities(g==(h-1));
        datP(isnan(datP))=[];
        don = ones(numel(datP),1)*(h-1);
        scatter(don,datP,'blue','filled')
        pAves(h) = mean(datP,'omitnan');
        pSTDs(h) = std(datP,'omitnan');
    end
    plot(0:4,pAves)
    % figure
    % errorbar(0:4,pAves,pSTDs)

    % aDays = unique(g);
    title('Neural Activation density (initial 500um) At threshold')
    ylabel('Activation Density ')
    xlabel('Weeks of Training')
    xticks(0:4)
    xlim([-0.5 4.5])
    % xticklabels(num2str(aDays'))


    ranksum(threshDensities(g<2),threshDensities(g>1))
    % [~,~,stats] = anovan(threshDensities,g2);

%% Threshold  ring plots


    maxDen3DCOM = [];
    for t = 1:5 % For each week of training
        weekDensityData = squeeze(weekThreshDense{t});
        threshAves = mean(weekDensityData,1,'omitnan');
        maxDen3DCOM = max([maxDen3DCOM, threshAves]);
    end
    maxDen3DCOM = maxDen3DCOM*1.2;

    cmap = jet(256);
    % maxDen3DCOM = max(threshDensities,[],'all');
    for t = 1:5 % For each week of training
        weekDensityData = squeeze(weekThreshDense{t});
        threshAves = mean(weekDensityData,1,'omitnan');
        threshAves(threshAves==0) = 1e-15;
        if(~isnan(threshAves(1)))
            figName = strcat('Threshold Week:',{' '},num2str(t-1));
            figure('Name',figName{1},'NumberTitle','off');
            hold on
            for k = 10:-1:1
                x1 = 0.5*abs(k-11);
                pos = [x1 x1 k k]; 
                if(threshAves(k)==0)
                    cmapS = [1,1,1];
                else
                    cmapS = cmap(ceil(256*threshAves(k)/maxDen3DCOM),:);
                end
                rectangle('Position',pos,'Curvature',[1 1],'FaceColor',cmapS)
                axis equal 
            end
            caxis([0 maxDen3DCOM])
            colormap(jet)
            cB = colorbar;
            set(gca,'fontsize', 28)
            axis off
    
        end
    end
    






%%

    % % % % % 
    % % % % % 
    % % % % % 
    % % % % % % Plot activation density with respect to threshold (traces)
    % % % % % currColorsT  = jet(numel(daysThresholdOverall));
    % % % % % for c = 1:numel(uCurrOverall)
    % % % % %     figure
    % % % % %     hold on
    % % % % %     for t = 1:numel(daysThresholdOverall)
    % % % % %         threshAves = squeeze(threshDensityData(c,t,:));
    % % % % %         p1 = plot(1:10,threshAves,'-o');
    % % % % %         tcolor = currColorsT(t,:);
    % % % % %         p1.Color = tcolor;
    % % % % %         p1.LineWidth=2;
    % % % % %     end
    % % % % %     xticks(1:10)
    % % % % %     xticklabels({'0-100','100-200','200-300','300-400','400-500','500-600','600-700','700-800','800-900','900-1000'})
    % % % % %     ylabel('Neural Activation Density (neurons \mum^3)')
    % % % % %     xlabel('Distance from Electrode (\mum)')
    % % % % %     xtickangle(45)
    % % % % %     title(strcat('Activation Density Vs Threshold at',{' '},num2str(uCurrOverall(c)),'\muA'))
    % % % % % 
    % % % % % 
    % % % % %     lgd = legend('','Location', 'eastoutside');
    % % % % %     col_names = num2str(daysThresholdOverall);
    % % % % %     for j =1:length(col_names)
    % % % % %         plot([NaN NaN], [NaN NaN], 'Color', currColorsT(j,:), 'DisplayName', col_names(j))
    % % % % %     end
    % % % % %     title(lgd,'Threshold (\muA)')
    % % % % % end
    % % % % % 
    % % % % % 
    % % % % % 
    % % % % % % 
    % 
    % %% density at threshold  single
    %     % Structures for averaging current specific response
    %     % mergingTraces = NaN(numel(daysThresholdOverall),numel(anID)*3*numel(uChanOverall));
    %     % traceNum=0;
    % 
    %     grpA = [];
    %     grpB = [];
    %     % Process each animal
    %     for aInd = 1:numel(anID)
    %         uCurr = uCurrAll{aInd};
    %         uChan = uChanAll{aInd};
    %         daysThresholds = round(daysThresholdsAll{aInd});
    %         allDensityThresh = allDensityThreshAll{aInd};
    % 
    %         uThresholds = unique(unique(daysThresholds));
    %         uThresholds(uThresholds==0)=[];
    %         for cInd = 1:numel(uCurrOverall) % process each current
    %             curInd = find(uCurr==uCurrOverall(cInd));
    %             if(curInd>0) % if current current is present for given animal, continue analysis
    %                 for chInd = 1:numel(uChan)% Iterate through each channel
    %                     if(aInd==5 && chInd==3)
    %                         continue;
    %                         % Skip ICMS 100 channel 6
    %                     end
    %                     curDat = squeeze(allDensityThresh(chInd,curInd,:,1:5));
    %                     curThresh = uThresholds;
    % 
    %                     % Only plot if there is data
    %                     if(~isempty(curDat))
    % 
    % 
    %                         scatter(curThresh,curDat,[],currColorsTurbo(1:5,:),'filled');
    % 
    %                         % Add trace to current average
    %                         for d = 1:numel(curThresh)
    %                             % employ lazy abuse of memory to maintain data
    %                             % accuracy and avoid convoluted averaging
    %                             % algorithm.....  heh it works
    %                             alignedThresh = find(daysThresholdOverall==curThresh(d)); % Align thresholds
    %                             traceNum = traceNum+1;
    %                             mergingTraces(alignedThresh,traceNum) = curDat(d);
    %                         end
    %                     end 
    %                 end
    %             end
    %         end
    % 
    %         % Average response for current and plot
    %         curAve = mean(mergingTraces,2,'omitnan');
    %         curThreshs = daysThresholdOverall;
    %         curThreshs = curThreshs(~isnan(curAve));
    %         curAve = curAve(~isnan(curAve));
    % 
    %         % Add means to matrix for alternative plotting
    %         threshInds = zeros(numel(curThreshs),1);
    %         for k = 1:numel(curThreshs)
    %             threshInds(k) = find(daysThresholdOverall==curThreshs(k));
    %         end
    % 
    %         threshDensityData(cInd,threshInds,di) = curAve;
    %     end
    % end
    % 
    % 




    
    
    
    
    
    
    
% % % %     
% % % %     % Process behavioral thresholds for each channel over time and plot trend
% % % %     figure()
% % % %     hold on
% % % %     weeksTrainedOverall = unique(ceil(daysTrainedOverall/dDiv));
% % % % 
% % % %     mergingThreshold = NaN(numel(daysTrainedOverall),numel(anID)*3);
% % % %     mtC = 0;
% % % %     for aInd = 1:numel(anID)
% % % %         uChan = uChanAll{aInd};
% % % %         daysThresholds = daysThresholdsAll{aInd};
% % % %         weeksTrained = ceil(daysTrainedAll{aInd}/dDiv);
% % % %         alignedDays = zeros(numel(weeksTrained),1);
% % % %         for d = 1:numel(weeksTrained)
% % % %             alignedDays(d) = find(weeksTrainedOverall==weeksTrained(d)); % Align days
% % % %         end
% % % % 
% % % %         for chInd = 1:numel(uChan) % process each channel
% % % %             if(aInd==5 && chInd==3)
% % % %                 continue;
% % % %                 % Skip ICMS 100 channel 6
% % % %             end
% % % %             mtC = mtC+1;
% % % %             mergingThreshold(alignedDays,mtC) = daysThresholds(chInd,:);
% % % % 
% % % %             p2 = plot(weeksTrained,daysThresholds(chInd,:),'-',"Color", [0, 0, 0, 0.3]);
% % % %             lgN = lgN + 1;
% % % %             legText{lgN} = strcat('Animal',num2str(anID(aInd)),'Ch',num2str(uChan(chInd)),'Threshold');
% % % % 
% % % % 
% % % %         end
% % % %     end
% % % %     % Average response for threshold and plot
% % % %     aveThresh = mean(mergingThreshold,2,'omitnan');
% % % %     curDays = daysTrainedOverall;
% % % %     curDays = curDays(~isnan(aveThresh));
% % % %     aveThresh = aveThresh(~isnan(aveThresh));
% % % % 
% % % % 
% % % %     % Set up fittype and options.
% % % % %     ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
% % % % %     opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% % % % %     opts.Display = 'Off';
% % % % %     opts.MaxFunEvals = 20;
% % % % %     opts.MaxIter = 20;
% % % % %     opts.Robust = 'LAR';
% % % % %     opts.StartPoint = [100 0.001 2];
% % % % 
% % % %     % Fit model to data.
% % % % %     [fitresult, ~] = fit(curDays',aveThresh, ft, opts );
% % % % %     p1 = plot(fitresult);
% % % %     p1 = plot(weeksTrainedOverall,aveThresh);
% % % %     p1.Color = 'k';
% % % %     p1.LineWidth=2;
% % % % 
% % % %     ylabel('Behavioral Threshold (\muA)')
% % % %     xlabel('Weeks of Training')
% % % % 
% % % %     
    
    
    
    
    

  %% Process behavioral thresholds for each channel over weeks and plot trend
    figure()
    hold on
    weeksTrainedOverall = unique(ceil(daysTrainedOverall/dDiv));
    mergingThreshold = NaN(numel(weeksTrainedOverall),numel(anID)*3*numel(weeksTrainedOverall));
    mtC = 0;
    for aInd = 1:numel(anID)
        uChan = uChanAll{aInd};
        daysThresholds = daysThresholdsAll{aInd};
        if(aInd==1)
            daysThresholds = daysThresholds([1,3],:);
        end
        daysThresholds(daysThresholds==0)=NaN;
        daysTrained = daysTrainedAll2{aInd};
        weeksTrained = ceil(daysTrained/dDiv);
        
        for chInd = 1:numel(uChan) % process each channel
            if(aInd==5 && chInd==3)
                continue;
                % Skip ICMS 100 channel 6
            end
            for d = 1:numel(weeksTrained)% Align weeks and add to merging dataset
                mtC = mtC+1; 
                mergingThreshold(weeksTrainedOverall==weeksTrained(d),mtC) = daysThresholds(chInd,d);
            end
            
            scatter(weeksTrained,squeeze(daysThresholds(chInd,:)),[],[0.5 0.5 0.5],'filled');
        end
    end
    % Average response for threshold and plot
    mergingThreshold(mergingThreshold==0) = NaN; % remove 0's generated by excess cell addition
    aveThresh = mean(mergingThreshold,2,'omitnan');
    threshWeeks = weeksTrainedOverall;
    threshWeeks = threshWeeks(~isnan(aveThresh));
    threshWeeks(threshWeeks>4) = NaN;
    aveThresh = aveThresh(~isnan(aveThresh));

    p1 = plot(threshWeeks',aveThresh);
    p1.Color = 'k';
    p1.LineWidth=2;
    ylabel('Behavioral Threshold (\muA)')
    xlabel('Weeks of Training')



    
    
    
    
    
    
    % NEURAL ACTIVATION AT THREHOLD OVER TIME ---------------------------------------------- 
    figure()
    hold on
    
    % Structures for averaging current specific response
    daysTrainedOverall = unique(cell2mat(daysTrainedAll'));
    weeksTrainedOverall = unique(ceil(daysTrainedOverall/dDiv));
    mergingTraces = NaN(numel(weeksTrainedOverall),numel(anID)*3);
    mergingThresholds = NaN(numel(weeksTrainedOverall),numel(anID)*3);
    mt=0;

    % Process each animal
    for aInd = 1:numel(anID)
        uChan = uChanAll{aInd};
        uCurr = uCurrAll{aInd};
        groupedCurrChans = CurrChansAll2{aInd};
        weeksTrained = ceil(daysTrainedAll2{aInd}/dDiv);
        daysThresholds = daysThresholdsAll{aInd};
        if(aInd==1)
            daysThresholds = daysThresholds([1,3],:);
        end
        daysThresholds(daysThresholds==0)=NaN;
        for chInd = 1:numel(uChan) % Process each channel
            if(aInd==5 && chInd==3)
                continue;
                % Skip ICMS 100 channel 6
            end
            chanNeuCnt = squeeze(groupedCurrChans(:,chInd,:));
            chanThresh = daysThresholds(chInd,:);
            threshNeuCnts = NaN(size(chanThresh));

            % calculate response at each day's detection threshold  and
            % plot, nearest neighbor
            for d = 1:numel(weeksTrained)
                ccInd = find(uCurr == round(chanThresh(d)));
                if(~isempty(ccInd))
                    threshNeuCnts(d) = chanNeuCnt(d,ccInd);
                end
            end
             
            % Plot points for activation at threshold
            yyaxis left
%             plot(weeksTrained,threshNeuCnts)

            s = scatter(weeksTrained,threshNeuCnts,[],'b','filled');
            s.MarkerFaceAlpha = 0.5;
            
            
            % Plot points for threshold
            yyaxis right
            chanThresh(isnan(threshNeuCnts))=NaN;
            s = scatter(weeksTrained+0.1,chanThresh,[],[0.8500 0.3250 0.0980],'filled');
            s.MarkerFaceAlpha = 0.5;
            
            % Save data to merging structure to examine effect over all
            % weeks, animals, and channels
            for d = 1:numel(weeksTrained)
                mt = mt+1;
                mergingTraces(weeksTrainedOverall==weeksTrained(d),mt) = threshNeuCnts(d);
                mergingThresholds(weeksTrainedOverall==weeksTrained(d),mt) = chanThresh(d);
            end
        end
    end
    
    % Average response for each week and plot neural population
    yyaxis left
    mergingTraces(mergingTraces==0)=NaN;
    curAve = mean(mergingTraces,2,'omitnan');
    curWeeks = weeksTrainedOverall;
    curWeeks(curWeeks>4)=NaN;
%     curWeeks(curWeeks<1)=NaN;

    
    
    plot(curWeeks',curAve,'b','LineWidth',2);

    title('Neural Activation at Threshold - NN')
    ylabel('Number of  Neurons')
    xlabel('Weeks of Training')
    
    % Change axis and plot average threshold
    yyaxis right
    mergingThresholds(mergingThresholds==0)=NaN;
    curThreshAve = mean(mergingThresholds,2,'omitnan');
    plot(curWeeks',curThreshAve,'Color',[0.8500 0.3250 0.0980],'LineWidth',2);
    ylabel('Detection Threshold (\muA)')
    
    
    
    
    
    
    
    
    
    
    
    
    %% NEURAL ACTIVATION AT THREHOLD OVER TIME - With traces ---------------------------------------------- 
    figure()
    hold on
    
    % Structures for averaging current specific response
    daysTrainedOverall = unique(cell2mat(daysTrainedAll2'));
    weeksTrainedOverall = unique(ceil(daysTrainedOverall/dDiv));
    mergingTraces = NaN(numel(weeksTrainedOverall),numel(anID)*3);
    mergingThresholds = NaN(numel(weeksTrainedOverall),numel(anID)*3);
    mt=0;

    % Process each animal
    for aInd = 1:numel(anID)
        uChan = uChanAll{aInd};
        uCurr = uCurrAll{aInd};
        groupedCurrChans = CurrChansAll2{aInd};
        weeksTrained = ceil(daysTrainedAll2{aInd}/dDiv);
        daysTrained = daysTrainedAll2{aInd}/dDiv;
        daysThresholds = daysThresholdsAll{aInd};
        if(aInd==1)
            daysThresholds = daysThresholds([1,3],:);
        end
        daysThresholds(daysThresholds==0)=NaN;
        for chInd = 1:numel(uChan) % Process each channel
            if(aInd==5 && chInd==3)
                continue;
                % Skip ICMS 100 channel 6
            end
            chanNeuCnt = squeeze(groupedCurrChans(:,chInd,:));
            chanThresh = daysThresholds(chInd,:);
            threshNeuCnts = NaN(size(chanThresh));

            % calculate response at each day's detection threshold  and
            % plot, nearest neighbor
            for d = 1:numel(weeksTrained)
                ccInd = find(uCurr == round(chanThresh(d)));
                if(~isempty(ccInd))
                    threshNeuCnts(d) = chanNeuCnt(d,ccInd);
                end
            end
             
            % Plot points for activation at threshold
            yyaxis left
            plot(daysTrained,threshNeuCnts,'-','Color',[0 0 1 0.3]);
            
            
            % Plot points for threshold
            yyaxis right
            chanThresh(isnan(threshNeuCnts))=NaN;
            plot(daysTrained,chanThresh,'-','Color',[0.8500 0.3250 0.0980 0.3]);
            
            % Save data to merging structure to examine effect over all
            % weeks, animals, and channels
            for d = 1:numel(weeksTrained)
                mt = mt+1;
                mergingTraces(weeksTrainedOverall==weeksTrained(d),mt) = threshNeuCnts(d);
                mergingThresholds(weeksTrainedOverall==weeksTrained(d),mt) = chanThresh(d);
            end
        end
    end
    
    % Average response for each week and plot neural population
    yyaxis left
    mergingTraces(mergingTraces==0)=NaN;
    curAve = mean(mergingTraces,2,'omitnan');
    curWeeks = weeksTrainedOverall;
    curWeeks(curWeeks>4)=NaN;
%     curWeeks(curWeeks<1)=NaN;

    
    
    plot(curWeeks',curAve,'b','LineWidth',2);

    title('Neural Activation at Threshold - NN')
    ylabel('Number of  Neurons')
    xlabel('Weeks of Training')
    
    % Change axis and plot average threshold
    yyaxis right
    mergingThresholds(mergingThresholds==0)=NaN;
    curThreshAve = mean(mergingThresholds,2,'omitnan');
    plot(curWeeks',curThreshAve,'Color',[0.8500 0.3250 0.0980],'LineWidth',2);
    ylabel('Detection Threshold (\muA)')
    
    
    
    
    
    
    
    
    
    
    
    
%     
%     % NEURAL ACTIVATION AT THREHOLD OVER TIME ---------------------------------------------- 
%     figure()
%     hold on
%     
%     % Structures for averaging current specific response
%     daysTrainedOverall = unique(cell2mat(daysTrainedAll2'));
%     weeksTrainedOverall = unique(ceil(daysTrainedOverall/dDiv));
%     mergingTraces = NaN(numel(weeksTrainedOverall),numel(anID)*3);
%     mergingThresholds = NaN(numel(weeksTrainedOverall),numel(anID)*3);
%     mt=0;
% 
%     % Process each animal
%     for aInd = 1:numel(anID)
%         uChan = uChanAll{aInd};
%         uCurr = uCurrAll{aInd};
%         groupedCurrChans = CurrChansAll2{aInd};
%         weeksTrained = ceil(daysTrainedAll2{aInd}/dDiv);
%         daysThresholds = daysThresholdsAll{aInd};
%         if(aInd==1)
%             daysThresholds = daysThresholds([1,3],:);
%         end
%         daysThresholds(daysThresholds==0)=NaN;
%         for chInd = 1:numel(uChan) % Process each channel
%             if(aInd==5 && chInd==3)
%                 continue;
%                 % Skip ICMS 100 channel 6
%             end
%             chanNeuCnt = squeeze(groupedCurrChans(:,chInd,:));
%             chanThresh = daysThresholds(chInd,:);
%             threshNeuCnts = NaN(size(chanThresh));
% 
%             % calculate response at each day's detection threshold  and
%             % plot
%             for d = 1:numel(weeksTrained)
%                 x = uCurr;
%                 v = squeeze(chanNeuCnt(d,:));
%                 x(isnan(v))=[];
%                 v(isnan(v))=[];
%                 xq = chanThresh(d);
%                 if(numel(v)>1)
%                     threshNeuCnts(d) = interp1(x,v,xq,'linear');
%                 end
%             end
%              
%             % Plot points for activation at threshold
%             yyaxis left
%             
%             s = scatter(weeksTrained,threshNeuCnts,[],'b','filled');
%             s.MarkerFaceAlpha = 0.5;
%             
%             
%             % Plot points for threshold
%             yyaxis right
%             chanThresh(isnan(threshNeuCnts))=NaN;
%             s = scatter(weeksTrained+0.1,chanThresh,[],[0.8500 0.3250 0.0980],'filled');
%             s.MarkerFaceAlpha = 0.5;
%             
%             % Save data to merging structure to examine effect over all
%             % weeks, animals, and channels
%             for d = 1:numel(weeksTrained)
%                 mt = mt+1;
%                 mergingTraces(weeksTrainedOverall==weeksTrained(d),mt) = threshNeuCnts(d);
%                 mergingThresholds(weeksTrainedOverall==weeksTrained(d),mt) = chanThresh(d);
%             end
%         end
%     end
%     
%     % Average response for each week and plot neural population
%     yyaxis left
%     mergingTraces(mergingTraces==0)=NaN;
%     curAve = mean(mergingTraces,2,'omitnan');
%     curWeeks = weeksTrainedOverall;
%     curWeeks(curWeeks>4)=NaN;
% %     curWeeks(curWeeks<1)=NaN;
% 
%     
%     
%     plot(curWeeks',curAve,'b','LineWidth',2);
% 
%     title('Neural Activation at Threshold')
%     ylabel('Number of  Neurons')
%     xlabel('Weeks of Training')
%     
%     % Change axis and plot average threshold
%     yyaxis right
%     mergingThresholds(mergingThresholds==0)=NaN;
%     curThreshAve = mean(mergingThresholds,2,'omitnan');
%     plot(curWeeks',curThreshAve,'Color',[0.8500 0.3250 0.0980],'LineWidth',2);
%     ylabel('Detection Threshold (\muA)')
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     



   %% NEURAL ACTIVATION DENSITY AT THREHOLD OVER TIME ---------------------------------------------- 
    figure()
    hold on

    % Structures for averaging current specific response
    daysTrainedOverall = unique(cell2mat(daysTrainedAll2'));
    weeksTrainedOverall = unique(ceil(daysTrainedOverall/dDiv));
    mergingTraces = NaN(numel(weeksTrainedOverall),numel(anID)*3);
    mergingThresholds = NaN(numel(weeksTrainedOverall),numel(anID)*3);
    mt=0;

    % Process each animal
    for aInd = 1:numel(anID)
        uChan = uChanAll{aInd};
        uCurr = uCurrAll{aInd};
        groupedCurrChans = CurrChansAll2{aInd};
        weeksTrained = ceil(daysTrainedAll2{aInd}/dDiv);
        daysThresholds = daysThresholdsAll{aInd};
        if(aInd==1)
            daysThresholds = daysThresholds([1,3],:);
        end
        daysThresholds(daysThresholds==0)=NaN;
        for chInd = 1:numel(uChan) % Process each channel
            if(aInd==5 && chInd==3)
                continue;
                % Skip ICMS 100 channel 6
            end
            chanNeuCnt = squeeze(groupedCurrChans(:,chInd,:));
            chanThresh = daysThresholds(chInd,:);
            threshNeuCnts = NaN(size(chanThresh));

            % calculate response at each day's detection threshold  and
            % plot
            for d = 1:numel(weeksTrained)
                x = uCurr;
                v = squeeze(chanNeuCnt(d,:));
                x(isnan(v))=[];
                v(isnan(v))=[];
                xq = chanThresh(d);
                if(numel(v)>1)
                    threshNeuCnts(d) = interp1(x,v,xq,'linear');
                end
            end

            % Plot points for activation at threshold
            s = scatter(weeksTrained,threshNeuCnts,[],'b','filled');
            s.MarkerFaceAlpha = 0.5;


            % Save data to merging structure to examine effect over all
            % weeks, animals, and channels
            for d = 1:numel(weeksTrained)
                mt = mt+1;
                mergingTraces(weeksTrainedOverall==weeksTrained(d),mt) = threshNeuCnts(d);
                mergingThresholds(weeksTrainedOverall==weeksTrained(d),mt) = chanThresh(d);
            end
        end
    end

    % Average response for each week and plot neural population
    yyaxis left
    mergingTraces(mergingTraces==0)=NaN;
    curAve = mean(mergingTraces,2,'omitnan');
    curWeeks = weeksTrainedOverall;
    curWeeks(curWeeks>4)=NaN;


    plot(curWeeks',curAve,'b','LineWidth',2);

    title('Neural Activation at Threshold')
    ylabel('Number of  Neurons')
    xlabel('Weeks of Training')

    % Change axis and plot average threshold
    yyaxis right
    mergingThresholds(mergingThresholds==0)=NaN;
    curThreshAve = mean(mergingThresholds,2,'omitnan');
    plot(curWeeks',curThreshAve,'Color',[0.8500 0.3250 0.0980],'LineWidth',2);
    ylabel('Detection Threshold (\muA)')


    
    


    %% NEURAL ACTIVATION DENSITY AT THREHOLD OVER TIME ---------------------------------------------- 
    figure()
    hold on

    % Structures for averaging current specific response
    daysTrainedOverall = unique(cell2mat(daysTrainedAll2'));
    weeksTrainedOverall = unique(ceil(daysTrainedOverall/dDiv));
    mergingTraces = NaN(numel(weeksTrainedOverall),numel(anID)*3);
    mergingThresholds = NaN(numel(weeksTrainedOverall),numel(anID)*3);
    mt=0;

    % Process each animal
    for aInd = 1:numel(anID)
        uChan = uChanAll{aInd};
        uCurr = uCurrAll{aInd};
        groupedCurrChans = CurrChansAll2{aInd};
        weeksTrained = ceil(daysTrainedAll2{aInd}/dDiv);
        daysThresholds = daysThresholdsAll{aInd};
        if(aInd==1)
            daysThresholds = daysThresholds([1,3],:);
        end
        daysThresholds(daysThresholds==0)=NaN;
        for chInd = 1:numel(uChan) % Process each channel
            if(aInd==5 && chInd==3)
                continue;
                % Skip ICMS 100 channel 6
            end
            chanNeuCnt = squeeze(groupedCurrChans(:,chInd,:));
            chanThresh = daysThresholds(chInd,:);
            threshNeuCnts = NaN(size(chanThresh));

            % calculate response at each day's detection threshold  and
            % plot
            for d = 1:numel(weeksTrained)
                x = uCurr;
                v = squeeze(chanNeuCnt(d,:));
                x(isnan(v))=[];
                v(isnan(v))=[];
                xq = chanThresh(d);
                if(numel(v)>1)
                    threshNeuCnts(d) = interp1(x,v,xq,'linear');
                end
            end

            % Plot points for activation at threshold
            yyaxis left

            s = scatter(weeksTrained,threshNeuCnts,[],'b','filled');
            s.MarkerFaceAlpha = 0.5;


            % Plot points for threshold
            yyaxis right
            chanThresh(isnan(threshNeuCnts))=NaN;
            s = scatter(weeksTrained+0.1,chanThresh,[],[0.8500 0.3250 0.0980],'filled');
            s.MarkerFaceAlpha = 0.5;

            % Save data to merging structure to examine effect over all
            % weeks, animals, and channels
            for d = 1:numel(weeksTrained)
                mt = mt+1;
                mergingTraces(weeksTrainedOverall==weeksTrained(d),mt) = threshNeuCnts(d);
                mergingThresholds(weeksTrainedOverall==weeksTrained(d),mt) = chanThresh(d);
            end
        end
    end

    % Average response for each week and plot neural population
    yyaxis left
    mergingTraces(mergingTraces==0)=NaN;
    curAve = mean(mergingTraces,2,'omitnan');
    curWeeks = weeksTrainedOverall;
    curWeeks(curWeeks>4)=NaN;


    plot(curWeeks',curAve,'b','LineWidth',2);

    title('Neural Activation at Threshold')
    ylabel('Number of  Neurons')
    xlabel('Weeks of Training')

    % Change axis and plot average threshold
    yyaxis right
    mergingThresholds(mergingThresholds==0)=NaN;
    curThreshAve = mean(mergingThresholds,2,'omitnan');
    plot(curWeeks',curThreshAve,'Color',[0.8500 0.3250 0.0980],'LineWidth',2);
    ylabel('Detection Threshold (\muA)')


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %% ACTIVATED DENSITY vs Threshold  
    currColorsTurbo = jet(10);
    uCurrOverall = unique(cell2mat(uCurrAll));
    uCurrOverall(uCurrOverall>6)=[];
    uChanOverall = unique(cell2mat(uChanAll));
    daysThresholdOverall = unique(unique(round(cell2mat(daysThresholdsAll'))));
    daysThresholdOverall(daysThresholdOverall==0)=[];
    daysThresholdOverall(daysThresholdOverall>7)=[];
    
    
    
    threshDensityData = NaN(numel(uCurrOverall),numel(daysThresholdOverall),10);
    
    for cInd = 1:numel(uCurrOverall) % process each current
%         figure
%         hold on
        
        for di = 1:10 % process each distance from electrode
            % Structures for averaging current specific response
            mergingTraces = NaN(numel(daysThresholdOverall),numel(anID)*3*numel(uChanOverall));
            traceNum=0;

            % Process each animal
            for aInd = 1:numel(anID)
                uCurr = uCurrAll{aInd};
                uChan = uChanAll{aInd};
                daysThresholds = round(daysThresholdsAll{aInd});
                allDensityThresh = allDensityThreshAll{aInd};
                uThresholds = unique(unique(daysThresholds));
                uThresholds(uThresholds==0)=[];

                curInd = find(uCurr==uCurrOverall(cInd));
                if(curInd>0) % if current current is present for given animal, continue analysis
                    for chInd = 1:numel(uChan)% Iterate through each channel
                        if(aInd==5 && chInd==3)
                            continue;
                            % Skip ICMS 100 channel 6
                        end
                        curDat = squeeze(allDensityThresh(chInd,curInd,:,di));
                        curThresh = uThresholds;

                        % Only plot if there is data
                        if(~isempty(curDat))


                            
                            scatter(curThresh,curDat,[],currColorsTurbo(di,:),'filled');

                            % Add trace to current average
                            for d = 1:numel(curThresh)
                                % employ lazy abuse of memory to maintain data
                                % accuracy and avoid convoluted averaging
                                % algorithm.....  heh it works
                                alignedThresh = find(daysThresholdOverall==curThresh(d)); % Align thresholds
                                traceNum = traceNum+1;
                                mergingTraces(alignedThresh,traceNum) = curDat(d);
                            end
                        end 
                    end
                end
            end
            
            % Average response for current and plot
            curAve = mean(mergingTraces,2,'omitnan');
            curThreshs = daysThresholdOverall;
            curThreshs = curThreshs(~isnan(curAve));
            curAve = curAve(~isnan(curAve));
            
            % Add means to matrix for alternative plotting
            
            threshInds = zeros(numel(curThreshs),1);
            for k = 1:numel(curThreshs)
                threshInds(k) = find(daysThresholdOverall==curThreshs(k));
            end
            
            threshDensityData(cInd,threshInds,di) = curAve;
%             p1 = plot(curThreshs',curAve);
%             tcolor = currColorsTurbo(di,:);
%             p1.Color = tcolor(1:3);
%             p1.LineWidth=2;
        end
%         curCurr = uCurrOverall(cInd);
%         title(strcat('Activation Density Vs Threshold at',{' '},num2str(curCurr),'\muA'))
%         ylabel('Neural Activation Density (neurons \mum^3)')
%         xlabel('Detection Thresholds (\muA)')
%         
%         lgd = legend('','Location', 'eastoutside');
%         col_names = {'0-100','100-200','200-300','300-400','400-500','500-600','600-700','700-800','800-900','900-1000'};
%         for j =1:length(col_names)
%             plot([NaN NaN], [NaN NaN], 'Color', currColorsTurbo(j,:), 'DisplayName', col_names{j})
%         end
%         title(lgd,'Distance Bin (\mum)')
    end

    
    
    % Plot activation density with respect to threshold (traces)
    currColorsT  = jet(numel(daysThresholdOverall));
    for c = 1:numel(uCurrOverall)
        figure
        hold on
        for t = 1:numel(daysThresholdOverall)
            threshAves = squeeze(threshDensityData(c,t,:));
            p1 = plot(1:10,threshAves,'-o');
            tcolor = currColorsT(t,:);
            p1.Color = tcolor(1:3);
            p1.LineWidth=2;
        end
        xticks(1:10)
        xticklabels({'0-100','100-200','200-300','300-400','400-500','500-600','600-700','700-800','800-900','900-1000'})
        ylabel('Neural Activation Density (neurons \mum^3)')
        xlabel('Distance from Electrode (\mum)')
        xtickangle(45)
        title(strcat('Activation Density Vs Threshold at',{' '},num2str(uCurrOverall(c)),'\muA'))

        
        lgd = legend('','Location', 'eastoutside');
        col_names = num2str(daysThresholdOverall);
        for j =1:length(col_names)
            plot([NaN NaN], [NaN NaN], 'Color', currColorsT(j,:), 'DisplayName', col_names(j))
        end
        title(lgd,'Threshold (\muA)')
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
    
    
    %% Data Averaged over weeks
    % ACTIVATED POPULATION --------------------------------------
    figure()
    hold on
    legText = cell(1);
    lgN = 0;
    hLeg = legend('example');


    uCurrOverall = unique(cell2mat(uCurrAll));
    uCurrOverall(uCurrOverall>6)=[];
        uCurrOverall(uCurrOverall<2)=[];
    uChanOverall = unique(cell2mat(uChanAll));
    daysTrainedOverall = unique(cell2mat(daysTrainedAll'));
    daysTrainedOverall = unique(ceil(daysTrainedOverall/dDiv));
    for cInd = 1:numel(uCurrOverall) % process each current

        % Structures for averaging current specific response
        mergingTraces = NaN(numel(daysTrainedOverall),numel(anID)*numel(uChanOverall));
        traceNum=0;

        % Process each animal
        for aInd = 1:numel(anID)
            uChan = uChanAll{aInd};
            uCurr = uCurrAll{aInd};
            groupedCurrChans = CurrChansAll{aInd};
            daysTrained = daysTrainedAll{aInd};

            curInd = find(uCurr==uCurrOverall(cInd));
            if(curInd>0) % if current current is present for given animal, continue analysis
                for chInd = 1:numel(uChan) % process each channel
                    if(aInd==5 && chInd==3)
                        continue;
                        % Skip ICMS 100 channel 6
                    end
                    allChNeuCnt = groupedCurrChans(:,chInd,:);
                    mxChCnt = max(allChNeuCnt,[],'all');
                    neuCnt = groupedCurrChans(:,chInd,curInd)./mxChCnt;
                    curDays = daysTrained(~isnan(neuCnt));
                    neuCnt = neuCnt(~isnan(neuCnt));

                    % Only plot if there is data
                    if(~isempty(neuCnt))
                        %update neuron counts relative to maximum activation
                        %for each given channel
                        curWeeks = ceil(curDays/dDiv);

                        % Group values from same week together
                        uW = unique(curWeeks);
                        tNeuCnt = numel(uW);
                        for w = 1:numel(uW)
                           tNeuCnt(w) = mean(neuCnt(curWeeks==uW(w)));
                        end
                        neuCnt = tNeuCnt;
                        curDays = uW;
                                                
                        scatter(curDays,neuCnt,[],currColors{cInd},'filled');
%                         p1 = plot(curDays, neuCnt,'-'); % plot the number of neurons active for this channel and current
%                         p1.Color = currColors{cInd};
                        lgN = lgN + 1;
                        legText{lgN} = strcat('Animal',num2str(anID(aInd)),'Ch',num2str(uChan(chInd)),'--Current',num2str(uCurr(curInd)));

                        % Add trace to current average
                        traceNum = traceNum+1;
                        alignedDays = zeros(numel(curDays),1);
                        for d = 1:numel(curDays)
                            alignedDays(d) = find(daysTrainedOverall==curDays(d)); % Align days
                        end
                        mergingTraces(alignedDays,traceNum) = neuCnt;
                    end 
                end 
            end
        end
        % Average response for current and plot
        curAve = mean(mergingTraces,2,'omitnan');
        curDays = daysTrainedOverall;
        curDays = curDays(~isnan(curAve));
        curAve = curAve(~isnan(curAve));

    %     f=fit(curDays',curAve,'exp1');
    %     p1 = plot(f);
        p1 = plot(curDays',curAve);
        tcolor = currColors{cInd};
        p1.Color = tcolor(1:3);
        p1.LineWidth=2;
        lgN = lgN + 1;
        legText{lgN} = strcat('Fitted Trend for Current =\ ',num2str(uCurrOverall(cInd)),'\muA');

    end
    title('Neural Activation over Time')
    ylabel('Number of  Neurons')
    xlabel('Weeks of Training')
    legend(legText)
    set(hLeg,'visible','off')

   
    

end
