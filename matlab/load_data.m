% NOTE: This script is for regenerating figure_source_data.mat from the raw
% pipeline. It requires access to the lab network drive (S:\Roy) and raw
% .mat files from the imaging pipeline.
%
% For figure reproduction, this script is NOT needed. Instead, download
% figure_source_data.mat from Zenodo and place it in data/matlab/.
% The figure scripts (fig2/, fig3/, fig_s2/) load from that file directly.
%
% Load processed data (raw pipeline - requires network access)
cd('S:\Roy')
load('2P Image Segmentation_12-20-24\2P Image Segmentation\debug_3-16-24.mat');

sourceFolders = {sourceFolders83, sourceFolders92, sourceFolders93, sourceFolders98, sourceFolders100, sourceFolders101};
vertStep = 25; % electrode spacing information
contactPos = {contactPos83,contactPos92,contactPos93,contactPos98,contactPos100,contactPos101};

pathR{1} = 'robin\data\raster83.mat';
pathR{2} = 'robin\data\raster92.mat';
pathR{3} = 'robin\data\raster93.mat';
pathR{4} = 'robin\data\raster98.mat';
pathR{5} = 'robin\data\raster100.mat';
pathR{6} = 'robin\data\raster101.mat';

% pathF and pathC mat files contain BOTH hit and miss 
pathF{1} = 'AlternativeConcistencyDataMiss_ver6MISS83.mat';
pathF{2} = 'AlternativeConcistencyDataMiss_ver6MISS92.mat';
pathF{3} = 'AlternativeConcistencyDataMiss_ver6MISS93.mat';
pathF{4} = 'AlternativeConcistencyDataMiss_ver6MISS98.mat';
pathF{5} = 'AlternativeConcistencyDataMiss_ver6MISS100.mat';
pathF{6} = 'AlternativeConcistencyDataMiss_ver6MISS101.mat';

pathC{1} = 'AlternativeConcistencyDataMiss_ver6MISS83.mat';
pathC{2} = 'AlternativeConcistencyDataMiss_ver6MISS92.mat';
pathC{3} = 'AlternativeConcistencyDataMiss_ver6MISS93.mat';
pathC{4} = 'AlternativeConcistencyDataMiss_ver6MISS98.mat';
pathC{5} = 'AlternativeConcistencyDataMiss_ver6MISS100.mat';
pathC{6} = 'AlternativeConcistencyDataMiss_ver6MISS101.mat';

% load data 
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
    uChanCurrsAll{i} = tempData.uChanCurrs;
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

daysThresholdsAll{1} = daysThresholds83; % Manually loaded day thresholds
daysThresholdsAll{2} = daysThresholds92;
daysThresholdsAll{3} = daysThresholds93;
daysThresholdsAll{4} = daysThresholds98;
daysThresholdsAll{5} = daysThresholds100;
daysThresholdsAll{6} = daysThresholds101;

%%
uCurrOverallExaust = unique(cell2mat(uCurrAll));
densityAll = cell(numel(anID),1); 
densityAllAgr = cell(numel(anID),1);
densityAllAgr50 = cell(numel(anID),1);
r50Array = cell(numel(anID),1);
 
% Process each animal
for aInd = 1:numel(anID)
    % if aInd == 1; continue; end;
    uChan = uChanAll{aInd};
    uCurr = uCurrAll{aInd};
    daysTr = daysTrainedAll{aInd};
    uChanCurrs = uChanCurrsRef{aInd};
    popDates = popDatesAll{aInd};
    populationDistances = populationDistancesAll{aInd};

    weekGrid0 = unique(ceil(daysTr./7));  
    anVols = allVols{aInd};
    densityCur = densityAll{aInd};
    densityAgr = densityAllAgr{aInd};
    densityAgr50 = densityAllAgr50{aInd};
    r50s = r50Array{aInd};

    if(isempty(densityCur))
        % initialize variable if needed
        densityCur = NaN(numel(uCurrOverallExaust),3,numel(daysTr),10);
        densityAgr = NaN(numel(uCurrOverallExaust),3,numel(daysTr));
        densityAgr50 = NaN(numel(uCurrOverallExaust),3,numel(daysTr));
        r50s = NaN(numel(uCurrOverallExaust),3,numel(daysTr));
    end
    

    for cInd = 1:numel(uCurrOverallExaust) % process each current over weeks
        curInd = find(uCurr==uCurrOverallExaust(cInd), 1);
        if isempty(curInd), continue; end
        for chInd = 1:numel(uChan) % process each channel
            % Skip ICMS 100 channel 6 (error was found in this dataset)
            if(aInd==5 && chInd==3); continue; end
                
            curCurr = uCurr(curInd);
            for p = 1:size(uChanCurrs,1)
                if uChanCurrs(p,2)==curCurr && uChanCurrs(p,1)==uChan(chInd) % if valid entry exists
                    curDistances = populationDistances{p};
                    nDates = popDates{p};
                    curWeeks0 = ceil(daysTr(nDates)./7);
                    [tf, wkRows] = ismember(curWeeks0, weekGrid0);

                    curDays = daysTr(nDates);  
                    curDays = ceil(curDays/7);
                    for wkRow = 1:numel(curDays)
                        d = curDays(wkRow);
                        activeDists = curDistances{wkRow};
                        activeDists(activeDists==0)=[];                  
                        mergedDensity = zeros(10,1);
                        if(sum(activeDists<1000)>1)
                            for m = 1:10
                                curcnt = activeDists>=(m-1)*100 & activeDists<m*100;
                                mergedDensity(m) = sum(curcnt)/(((4/3)*pi*(0.1*m)^3)-((4/3)*pi*(0.1*(m-1))^3));  % currently calculated in mm rather than um
                            end
                            densityCur(cInd,chInd,wkRow,:) = mergedDensity;
                            
                            % density within 500 um
                            in500um = (activeDists >= 0 & activeDists < 500);
                            n500 = sum(in500um);
                            vol500  = (4/3)*pi*(0.5)^3; 
                            densityAgr(cInd, chInd, wkRow) = n500 / vol500;

                            % density within r50, radius that encloses 50% of activated neurons
                            r50_um = median(activeDists);
                            r50_mm = r50_um / 1000; % convert µm -> mm
                            inR50 = (activeDists >= 0 & activeDists <= r50_um);
                            nR50 = sum(inR50);    
                            volR50  = (4/3)*pi*(r50_mm)^3; 
                            densityAgr50(cInd, chInd, wkRow) = nR50 / volR50;

                            % store r50_um 
                            r50s(cInd, chInd, wkRow) = r50_um;
                        end
                    end
                end
            end 
        end
    end  
    densityAllAgr{aInd} = densityAgr;
    densityAllAgr50{aInd} = densityAgr50;
    densityAll{aInd} = densityCur; 
    r50Array{aInd} = r50s;
end

%% Visualize histogram of r50s
figure;
for aInd = 1:6
    subplot(3, 2, aInd);
    histogram(r50Array{aInd}, 10);
    title('ICMS' + string(anID(aInd)));
    xlim([0, 1000])
    xlabel('Distance (um)')
    ylabel('Count')
    ylim([0, 30])
end
sgtitle('r50 histograms')

%%
%% Density at threshold - single plot
figure()
hold on
g = [];
threshDensities = [];
threshDensitiesGrpA = [];
threshDensitiesGrpB = [];
weekThreshDense = cell(10,1);

for aInd = 1:numel(anID)
    if aInd == 1; continue; end
    uChan = uChanAll{aInd};
    uCurr = uCurrAll{aInd};
    daysTrained = daysTrainedAll{aInd};
    weeksTrained = ceil(daysTrained./7);
    daysThresholds = daysThresholdsAll{aInd};
    curDensity = densityAllAgr50{aInd};

    for chInd = 1:numel(uChan) % process each channel
        longDens = NaN(numel(daysTrained),1);
        for cDay = 1:numel(daysTrained) % process each day
            curThresh = round(daysThresholds(chInd,cDay));
            % curInd = find(uCurr==curThresh);
            % if(~isempty(curInd))
            if(curThresh>0)
                % curDat = mean(curDensity(curThresh,chInd,cDay,1:5),'omitnan');
                curDat = curDensity(curThresh, chInd, cDay);
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
% title('Neural Activation density (initial 500um) At threshold')

title('Neural activation density (with r50) at threshold')

ylabel('Activation Density ')
xlabel('Weeks of Training')
xticks(0:4)
xlim([-0.5 4.5])
% xticklabels(num2str(aDays'))


ranksum(threshDensities(g<2),threshDensities(g>1))
% [~,~,stats] = anovan(threshDensities,g2);

%% ========== SAVE INTERMEDIARY DATA FOR REPRODUCIBILITY ==========
% Saves all workspace variables needed by figure scripts (Fig2C, Fig2D,
% Fig2E, Fig2F, Fig2G, SFig2D, SFig2E, SFig2F, etc.) into a single .mat.
% Upload this file for reviewers / collaborators to reproduce figures
% without needing access to the raw data pipeline.

savePath = fullfile(fileparts(mfilename('fullpath')), 'data');
if ~exist(savePath, 'dir'), mkdir(savePath); end

saveFile = fullfile(savePath, sprintf('figure_source_data_%s.mat', datestr(now, 'yyyy_mm_dd')));

save(saveFile, ...
    ... % --- Animal identifiers ---
    'anID', ...
    ... % --- Channel / current metadata ---
    'uChanAll', 'uCurrAll', 'uChanCurrsRef', 'uChanCurrsAll', ...
    'uCurrOverallExaust', ...
    ... % --- Neuron counts (unconstrained, from original pipeline) ---
    'CurrChansAll', 'CurrChansAll2', ...
    'groupedCurrChansAll', 'groupedCurrChansThreshAll', ...
    'groupedCurrChansConsistencyAll', ...
    'groupedCurrChansBaseMeanAll', 'groupedCurrChansBaseSTDAll', ...
    'groupedCurrChansBaseMeanUpper10All', ...
    'groupedCurrChansThreshBaseMeanAll', 'groupedCurrChansThreshBaseSTDAll', ...
    'groupedCurrChansThreshBaseMeanUpper10All', ...
    ... % --- Temporal data ---
    'daysTrainedAll', 'daysTrainedAll2', 'daysThresholdsAll', ...
    ... % --- Population distances (needed for 800 µm constraint) ---
    'populationDistancesAll', 'populationDistanceAll', ...
    'populationHistogramAll', 'popDatesAll', ...
    ... % --- Density data ---
    'densityAll', 'densityAllAgr', 'densityAllAgr50', 'r50Array', ...
    'allDensityAll', 'allDensityThreshAll', ...
    'allPopulationAll', 'enclosingRadiusAll', 'allRadiiThreshAll', ...
    ... % --- Correlation / structural data ---
    'correlationDates', 'correlationMatricies', 'correlationMatricies3D', ...
    'correlationMatriciesBase', 'correlationMatricies3DBase', ...
    'correlationMatriciesROI', 'correlationMatriciesROI3D', ...
    'correlationMatriciesNeuronCount', 'correlationMatriciesNeuronCountNorm', ...
    'diffMatricies', 'diffMatriciesSTD', 'diffMatriciesVals', 'diffMatricies3D', ...
    ... % --- Region / raster / volume data ---
    'allRegions', 'allVols', 'allRasters', ...
    'sharedROIMappingAll', 'neighborComparisions', ...
    ... % --- Electrode geometry ---
    'contactPos', 'vertStep', 'sourceFolders', ...
    '-v7.3');  % Use v7.3 for large files

fprintf('\n=== REPRODUCIBILITY DATA SAVED ===\n');
fprintf('  File: %s\n', saveFile);
d = dir(saveFile);
fprintf('  Size: %.1f MB\n', d.bytes / 1e6);
fprintf('  Contains all variables needed by figure scripts.\n');
fprintf('  To reload: load(''%s'')\n', saveFile);