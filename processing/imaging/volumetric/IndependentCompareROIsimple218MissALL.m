                        

% Process all datasets listed in results files and quantify neural response
function IndependentCompareROIsimple218MissALL(sourceFolders,daysTrained,daysThresholds,elecPts,icms)
%     javaaddpath 'C:\Program Files\MATLAB\R2019a\java\mij.jar'
%     javaaddpath 'C:\Program Files\MATLAB\R2019a\java\ij.jar'
%     addpath(genpath('./mij'));
	roiMode=1;
    verbose=0;
    
    switch icms
        case 83
            refInd = 8;
            ad = 3;
%             vOFF = [1 1 1 1 1 1 0 1 1];
            vOFF = [0 0 0 0 0 0 0 0 0];
        case 92
            refInd = 6;
            ad = 3;
            vOFF = [0 0 0 0 0 0 0 0 0];
        case 93
            refInd = 6;
            ad = 3;
            vOFF = [0 0 0 0 0 0 0 0 0 0];
        case 98
            refInd = 5; %need to change
            ad = 5;
            vOFF = [0 0 0 0 0 0 0];
        case 100
            refInd = 1;
            ad = 3;
            vOFF = [0 0 0 0 0];
        case 101 % issue with days thresholds - ln 1647
            refInd = 4;
            ad = 10;
            vOFF = [0 0 0 0 0 0];
    end
    
    
    % check that the analysis will be performed over multiple trials.
    % Remove options that do not have sufficent number of entries per
    % parameter set
    paramMin = 4;
    dfThresh = 5;

    weeksTrained = ceil(daysTrained/7);
    weekLim = max(weeksTrained);
    daysThresholdsOrig = daysThresholds;
    daysThresholds = round(daysThresholds);
    
    
    % electrode spacing information
    vertStep = 25;
    contactsOrdered = [9,13,11,29,13,27,15,25,7,17,5,19,3,21,1,23,10,32,12,30,14,28,16,26,8,18,6,20,4,22,2,24];
    if(icms>90)
        contactsOrdered = 1:32;
    end


    
    % Setup Output Folder
    selpathOutDist = '\\10.129.151.108\nei1\xieluanlabs\xl_stimulation\ICMS00\';
    selpathOut = strcat(selpathOutDist,'ROIAnalysis_',datestr(datetime('now'),'mm-dd-yyyy_HH-MM'));
    mkdir(selpathOut); %Create output folder
   

    % find datafolders for each dataset
    trueSourceFolders = cell(numel(sourceFolders) ,1);
    cN = 0;
    for curFile = 1:numel(sourceFolders) 
        % find analysis folders
        selpath = sourceFolders{curFile};
        idcs = strfind(selpath,'\');
        pullPath = selpath(1:idcs(end)-1);
        segListing = dir(pullPath);
        for s = 1:numel(segListing)
            segName = segListing(s).name;
%             if(contains(segName, 'A1'))
            if(contains(segName, 'D4'))
                cN = cN+1;
                trueSourceFolders{curFile} = strcat(pullPath,'\',segName);
                cN = cN+1;
                trueSourceFolders{curFile} = strcat(pullPath,'\',segName);
                break
            end
        end
    end




    %% Iterater through all datasets and align each to the first dataset
    %  While aligning data also check check what channels and currents were
    %  examined
    transforms = cell(numel(trueSourceFolders),1);
    chanCurrPerms = zeros(numel(trueSourceFolders),5);
    ccpC = 1;
    imgDepths = zeros(numel(trueSourceFolders),1);
    
    % Manually Select which of the datasets will serve as the reference
    % image for alignment (does not have to be the first)
    prestimMeans = zeros(numel(trueSourceFolders),1);
    prestimMeanVols = cell(numel(trueSourceFolders),1);
    for curFile = 1:numel(trueSourceFolders) 
        % Load baseline mean volume for registration 
        selpath = sourceFolders{curFile}; 
        idcs = strfind(selpath,'\');
        pullPath = selpath(1:idcs(end)-1);
        trialPath = strcat(pullPath,'\TRIALS2'); %construct folder for trial specific data
        if(~exist(trialPath, 'dir'))
            trialPath = strcat(pullPath,'\TRIALS5'); %construct folder for trial specific data, older format
        end
        destFile = strcat(trialPath,'\baselineSTD_Overall.mat');
        temp = load(destFile,'overallBaselineMean');
        curBaselineMean = temp.overallBaselineMean;
        imgDepths(curFile) = size(curBaselineMean,3);
        prestimMeans(curFile) = mean(curBaselineMean,'all');
        prestimMeanVols{curFile} = curBaselineMean;
    end
    sharedDepth = min(imgDepths);
    
    
    
    
    

    % Iterate through each valid output folder
    sourcePerms = cell(numel(trueSourceFolders),1);
    
    for curFile = 1:numel(trueSourceFolders) 
        curpath = trueSourceFolders{curFile}; 
        ccpCI = ccpC;
        
        channels = load(strcat(curpath,'\channels'));
        if(isfield(channels,'uChans'))
            channels = channels.uChans;
        else
            channels = channels.channels;
        end
        channels(channels==0)=[]; % remove 0 channels
    
        % select overall analysis output folder
        caListing = dir(curpath);
        caListing(1:2)=[];

        dirFlags = [caListing.isdir];
        if(mod(curFile,2)==1)
            HSS = 'HIT';
        else
            HSS = 'MISS';
        end
        caListing = caListing(dirFlags);

        
        % find currents analyzed on each channel
        for chI = 1:numel(channels)
            curCh = channels(chI);
            currs = zeros(numel(caListing),1);
            freqs = zeros(numel(caListing),1);
            pulses = zeros(numel(caListing),1);
            n=0;
            for j = 1:numel(caListing)
                temp = caListing(j).name;
                if(contains(temp,strcat('Ch',num2str(curCh))))
                    if(contains(temp,HSS))
                        % extract current value
                        parseName = strsplit(temp,'CUR');

                        if(contains(parseName{2},'FREQ'))
                            currName = strsplit(parseName{2},'FREQ');
                            freqName = strsplit(currName{2},'PULSE');
                            pulseName = strsplit(freqName{2},HSS);

                            n = n+1;
                            currs(n) = str2double(currName{1}); % current
                            freqs(n) = str2double(freqName{1}); 
                            pulses(n) = str2double(pulseName{1});
                        elseif(contains(parseName{2},'Freq'))
                            currName = strsplit(parseName{2},'Freq');
                            freqName = strsplit(currName{2},'Dur');
                            pulseName = strsplit(freqName{2},HSS);

                            n = n+1;
                            currs(n) = str2double(currName{1}); % current
                            freqs(n) = str2double(freqName{1}); 
                            pulses(n) = str2double(pulseName{1});
                        else
                            currName = strsplit(parseName{2},'.mat');
                            n = n+1;
                            currs(n) = str2double(currName{1}); % current
                            freqs(n) = 100; 
                            pulses(n) = 167;
                        end
                    end
                end
            end
            currs = currs(1:n);% trim array
            freqs = freqs(1:n);% trim array
            pulses = pulses(1:n);% trim array

            chanCurrPerms(ccpC:ccpC+(n-1),1) = curCh;
            chanCurrPerms(ccpC:ccpC+(n-1),2) = currs;
            chanCurrPerms(ccpC:ccpC+(n-1),3) = curFile;
            chanCurrPerms(ccpC:ccpC+(n-1),4) = freqs;
            chanCurrPerms(ccpC:ccpC+(n-1),5) = pulses;

            ccpC = ccpC + n;
        end    
        chanCurrPermsSubset = chanCurrPerms(ccpCI:ccpC-1,:);
        sourcePerms{curFile} = chanCurrPermsSubset;
    end

    
%% Perform secondary alignment - patch warp
    % Load baseline mean volume for registration 
    selpath = sourceFolders{refInd}; 
    idcs = strfind(selpath,'\');
    pullPath = selpath(1:idcs(end)-1);
    trialPath = strcat(pullPath,'\TRIALS2'); %construct folder for trial specific data
    if(~exist(trialPath, 'dir'))
        trialPath = strcat(pullPath,'\TRIALS5'); %construct folder for trial specific data, older format
    end
    destFile = strcat(trialPath,'\baselineSTD_Overall.mat');
    temp = load(destFile,'overallBaselineMean');
    
    
    
    curBaselineMean = temp.overallBaselineMean;
    image1_max = squeeze(max(curBaselineMean(:,:,1:ad),[],3));
    image1_mean = squeeze(mean(curBaselineMean(:,:,1:ad),3));

    

    ptForms = cell(numel(trueSourceFolders),1);
    for cd = 1:numel(trueSourceFolders) 

        curSource = sourceFolders{cd};
        fprintf('Aligning Dataset: ');
        pathDispName = erase(curSource,'\\10.129.151.108\nei1\xieluanlabs\xl_stimulation\');
        pathDispName = erase(pathDispName,'\RAW');
        disp(pathDispName)

        % Load baseline mean volume for registration 
        selpath = sourceFolders{cd}; 
        idcs = strfind(selpath,'\');
        pullPath = selpath(1:idcs(end)-1);
        trialPath = strcat(pullPath,'\TRIALS2'); %construct folder for trial specific data
        if(~exist(trialPath, 'dir'))
            trialPath = strcat(pullPath,'\TRIALS5'); %construct folder for trial specific data, older format
        end
        
        destFile = strcat(trialPath,'\baselineSTD_Overall.mat');
        temp = load(destFile,'overallBaselineMean');
        curBaselineMean = temp.overallBaselineMean;
        image2_max = squeeze(max(curBaselineMean(:,:,1:ad),[],3));
        image2_mean = squeeze(mean(curBaselineMean(:,:,1:ad),3));

        %% Normalize image intensity
        image1_mean = 1000*reshape(normalize(image1_mean(:), 'range'), size(image1_mean));
        image2_mean = 1000*reshape(normalize(image2_mean(:), 'range'), size(image2_mean));
        image1_max = 1000*reshape(normalize(image1_max(:), 'range'), size(image1_max));
        image2_max = 1000*reshape(normalize(image2_max(:), 'range'), size(image2_max));
        image1_all = cat(3, image1_mean, image1_max);
        image2_all = cat(3, image2_mean, image2_max);
        
        
        
        transform1 = 'euclidean';
        transform2 = 'affine';
        warp_blocksize = 6;
        warp_overlap_pix_frac = 0.25;
        norm_radius = 0;     % Set to 0 when the signals are sparse or the result does not look good.
        alignVal = inf;
        % run multiple iteration of the fitting and accept the one with the
        % smallest mis-alignment score
        for t = 1:10
        
            patchwarp_results = patchwarp_across_sessions(image1_all, image2_all,transform1, transform2, warp_blocksize, warp_overlap_pix_frac, norm_radius);
            shiftedplane = spatial_interp_patchwarp(image2_mean, patchwarp_results.warp1_cell{1}, 'euclidean', 1:512, 1:512);
            
            FOV = ones(512,512);
            shiftedFOV = spatial_interp_patchwarp(FOV, patchwarp_results.warp1_cell{1}, 'euclidean', 1:512, 1:512);

            refComp = image1_mean./max(image1_mean,[],'all');
            curShift = shiftedplane./max(shiftedplane,[],'all');
            compVol = refComp - curShift;
            compVol(shiftedFOV==0)=NaN;
            curAV = abs(mean(compVol,'all','omitnan'));
            
            if(t==1)
                ptForms{cd} = patchwarp_results.warp1_cell{1};
                alignVal = curAV;
            else
                % perform optimization check
                if(alignVal>curAV)
                    alignVal = curAV;
                    ptForms{cd} = patchwarp_results.warp1_cell{1};
                end
            end
        
        end        

        % figure('NumberTitle', 'off', 'Name', strcat('Sess',num2str(refInd),'--vs--Sess',num2str(cd)))
        % shiftedplane = spatial_interp_patchwarp(image2_mean, ptForms{cd}, 'euclidean', 1:512, 1:512);
        % imshow(imfuse(shiftedplane./max(shiftedplane,[],'all'),image1_mean./max(image1_mean,[],'all'),'falsecolor','Scaling','joint','ColorChannels','red-cyan'));
  
        shiftedVol = zeros(size(curBaselineMean));
        for d = 1:size(curBaselineMean,3)
            shiftedVol(:,:,d) = spatial_interp_patchwarp(curBaselineMean(:,:,d), ptForms{cd}, 'euclidean', 1:512, 1:512);
        end
        
%         figure('NumberTitle', 'off', 'Name', strcat('Sess',num2str(cd)));
%         imshow3D(shiftedVol)
    end

    
    
   
   
    
    
    

    %% Remove any sessions that were collected with varied pulse duration or
    %frequency
    chanCurrPerms(chanCurrPerms(:,4)~=100,:)=[];
    chanCurrPerms(chanCurrPerms(:,5)~=167,:)=[];
    
    % Remove channels that are only examined for a single session
    chDays = chanCurrPerms(:,[1,3]);
    chDays = unique(chDays,'rows');
    chCounts = chDays(:,1);
    possChans = unique(chCounts);
    for c = 1:numel(possChans)
        if(sum(possChans(c)==chCounts)<2)
            % Remove Invalid Channel
            chanCurrPerms(chanCurrPerms(:,1)==possChans(c),:)=[];
        end
    end
    chanCurrPermsOrig = chanCurrPerms; % use this in case I need to change naming for alignmnet but don't want to redo all the datastructures in the server
    
    
    uChanCurrs = unique(chanCurrPerms(:,1:2),'rows');
    invalidUPerm = false(size(chanCurrPerms,1),1);
    % Iterate through each unique combination counting repeats
    for curPerm = 1:size(uChanCurrs)
        [~, inds]=ismember(chanCurrPerms(:,1:2),uChanCurrs(curPerm,:),'rows');
        if(sum(inds)<paramMin)
            invalidUPerm(curPerm)=1;
        end
    end
    uChanCurrs(uChanCurrs(:,1)>100,:)=[]; % Remove repeat data



    % genrate overlay mask that all datasets are co-registered datasets
    imgSz = [512,512];
    overlay = true(imgSz);
	for t = 1:numel(transforms)
        moving = ones(imgSz);
        shiftedOverlay = spatial_interp_patchwarp(moving, ptForms{t}, 'euclidean', 1:512, 1:512);
        overlay(shiftedOverlay==0)=0;
	end


    overlay3D = false(512,512,sharedDepth);
    for d = 1:sharedDepth
       overlay3D(:,:,d)=overlay;
    end
    overlay3D = overlay3D(:,:,1:sharedDepth);
    fOverlay2D = reshape(overlay,1,[]);
    fOverlay3D = reshape(overlay3D,1,[]);
    vInds2D = find(fOverlay2D); % Start analysis off with averaged 2D
    vInds3D = find(fOverlay3D); % Extend to 3D as additional option

   
 %% Check alignment between datasets
    
    % Collect metrics to compare from across all datasets
    allMeans = cell(numel(trueSourceFolders),1);
    allROI = cell(numel(trueSourceFolders),1);
    for cd = 1:numel(trueSourceFolders)
        if(mod(cd,2)==1)
            HSS = 'HIT';
        else
            HSS = 'MISS';
        end


        % Load datasets from each aligned dataset to see alignment
         activeRows = chanCurrPermsOrig;
         vInd = find(activeRows(:,3)==cd);
         vInd = vInd(1);

        % Open ROI file for each datset
        curpath = trueSourceFolders{cd};
        fileNameROIH = strcat(curpath,'\Ch',num2str(activeRows(vInd,1)),'-CUR',num2str(activeRows(vInd,2)),HSS,'\ROI.mat'); 
        fileNameP = strcat(curpath,'\Ch',num2str(activeRows(vInd,1)),'-CUR',num2str(activeRows(vInd,2)),HSS,'\PreStim.mat'); 

        if(~isfile(fileNameROIH)) % check if older filename style is used or new format
            fileNameROIH = strcat(curpath,'\Ch',num2str(activeRows(vInd,1)),'-CUR',num2str(activeRows(vInd,2)),'FREQ',num2str(activeRows(vInd,4)),'PULSE',num2str(activeRows(vInd,5)),HSS,'\ROI.mat');
            fileNameP = strcat(curpath,'\Ch',num2str(activeRows(vInd,1)),'-CUR',num2str(activeRows(vInd,2)),'FREQ',num2str(activeRows(vInd,4)),'PULSE',num2str(activeRows(vInd,5)),HSS,'\PreStim.mat');
        end
        roiData = load(fileNameROIH);
        preData = load(fileNameP);

        curBaseMean = preData.meanPreStim; % Midline Flourescence Volume


        % All possible ROIs
        midActive = roiData.midActive;
        allROIs = false(size(midActive));
        regions = roiData.regions;
        for r = 1:numel(regions)
            curPixels = regions(r).PixelList;
            for m = 1:size(curPixels,1)
                allROIs(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
            end   
        end


        
        % Perform volumetric transformation
        for d = 1:sharedDepth
            allROIs(:,:,d) = spatial_interp_patchwarp(allROIs(:,:,d), ptForms{cd}, 'euclidean', 1:512, 1:512);
            curBaseMean(:,:,d) = spatial_interp_patchwarp(curBaseMean(:,:,d), ptForms{cd}, 'euclidean', 1:512, 1:512);
        end
        curBaseMean = curBaseMean(:,:,1:sharedDepth);        
        curBaseMean(overlay3D==0)=0;
        allROIs(overlay3D==0)=0;
        allMeans{cd} = curBaseMean;
        allROI{cd} = allROIs;
    end
    
    
    % % Perform comparisons of X, Y, and Z alignment of datasets 
    % currColors = {[1 0 0],[1 0 1],[0 0 1],[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0 0 0],[1 1 0],[0.4940 0.1840 0.5560],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]}; % Hopefully these sets of colors and symbols will be enough
    % figure('NumberTitle', 'off', 'Name', 'Baseline Mean Flourescence Alignment')
    % hold on
    % for cd = 1:numel(trueSourceFolders)
    % 
    %     subplot(2,2,1) % x alignment
    %     hold on
    %     p1 = plot(mean(mean(allMeans{cd},3),2));
    %     tcolor = currColors{cd};
    %     p1.Color = tcolor(1:3);
    %     p1.LineWidth=1;
    %     title('X Alignment')
    %     xlim([0 512])
    % 
    %     subplot(2,2,2) % y alignment
    %     hold on
    %     p1 = plot(mean(mean(allMeans{cd},3),1));
    %     tcolor = currColors{cd};
    %     p1.Color = tcolor(1:3);
    %     p1.LineWidth=1;
    %     title('Y Alignment')
    %     xlim([0 512])
    % 
    %     subplot(2,2,3) % z alignment
    %     hold on
    %     p1 = plot(squeeze(mean(mean(allMeans{cd},2),1)));
    %     tcolor = currColors{cd};
    %     p1.Color = tcolor(1:3);
    %     p1.LineWidth=1;
    %     title('Z Alignment')
    % end
    % sessN = 1:numel(trueSourceFolders);
    % legend('','Location', 'eastoutside');
    % col_names = num2str(sessN');
    % for j =1:length(col_names)
    %     plot([NaN NaN], [NaN NaN], 'Color', currColors{j}, 'DisplayName', strcat('Session', col_names(j)))
    % end
    
    
    if(verbose==1)
        figure('NumberTitle', 'off', 'Name', 'Baseline Normalized Flourescence Alignment')
        hold on
        for cd = 1:numel(trueSourceFolders)

            subplot(2,2,1) % x alignment
            hold on
            p1 = plot(mean(mean(allMeans{cd},3),2)./max(allMeans{cd},[],'all'));
            tcolor = currColors{cd};
            p1.Color = tcolor(1:3);
            p1.LineWidth=1;
            title('X Alignment')
            xlim([0 512])

            subplot(2,2,2) % y alignment
            hold on
            p1 = plot(mean(mean(allMeans{cd},3),1)./max(allMeans{cd},[],'all'));
            tcolor = currColors{cd};
            p1.Color = tcolor(1:3);
            p1.LineWidth=1;
            title('Y Alignment')
            xlim([0 512])

            subplot(2,2,3) % z alignment
            hold on
            p1 = plot(squeeze(mean(mean(allMeans{cd},2),1))./max(allMeans{cd},[],'all'));
            tcolor = currColors{cd};
            p1.Color = tcolor(1:3);
            p1.LineWidth=1;
            title('Z Alignment')
        end
        sessN = 1:numel(trueSourceFolders);
        legend('','Location', 'eastoutside');
        col_names = num2str(sessN');
        for j =1:length(col_names)
            plot([NaN NaN], [NaN NaN], 'Color', currColors{j}, 'DisplayName', strcat('Session', col_names(j)))
        end


        figure('NumberTitle', 'off', 'Name', 'ROI Distribution Alignment')
        hold on
        for cd = 1:numel(trueSourceFolders)

            subplot(2,2,1) % x alignment
            hold on
            p1 = plot(mean(mean(allROI{cd},3),2));
            tcolor = currColors{cd};
            p1.Color = tcolor(1:3);
            p1.LineWidth=1;
            title('X Alignment')
            xlim([0 512])

            subplot(2,2,2) % y alignment
            hold on
            p1 = plot(mean(mean(allROI{cd},3),1));
            tcolor = currColors{cd};
            p1.Color = tcolor(1:3);
            p1.LineWidth=1;
            title('Y Alignment')
            xlim([0 512])

            subplot(2,2,3) % z alignment
            hold on
            p1 = plot(squeeze(mean(mean(allROI{cd},2),1)));
            tcolor = currColors{cd};
            p1.Color = tcolor(1:3);
            p1.LineWidth=1;
            title('Z Alignment')
        end
        sessN = 1:numel(trueSourceFolders);
        legend('','Location', 'eastoutside');
        col_names = num2str(sessN');
        for j =1:length(col_names)
            plot([NaN NaN], [NaN NaN], 'Color', currColors{j}, 'DisplayName', strcat('Session', col_names(j)))
        end


        figure('NumberTitle', 'off', 'Name', 'ROI Distribution Alignment Smoothed')
        hold on
        for cd = 1:numel(trueSourceFolders)
            subplot(2,2,1) % x alignment
            hold on
            p1 = plot(movmean(mean(mean(allROI{cd},3),2),20));
            tcolor = currColors{cd};
            p1.Color = tcolor(1:3);
            p1.LineWidth=1;
            title('X Alignment')
            xlim([0 512])

            subplot(2,2,2) % y alignment
            hold on
            p1 = plot(movmean(mean(mean(allROI{cd},3),1),20));
            tcolor = currColors{cd};
            p1.Color = tcolor(1:3);
            p1.LineWidth=1;
            title('Y Alignment')
            xlim([0 512])

            subplot(2,2,3) % z alignment
            hold on
            p1 = plot(movmean(squeeze(mean(mean(allROI{cd},2),1)),2));
            tcolor = currColors{cd};
            p1.Color = tcolor(1:3);
            p1.LineWidth=1;
            title('Z Alignment')
        end
        sessN = 1:numel(trueSourceFolders);
        legend('','Location', 'eastoutside');
        col_names = num2str(sessN');
        for j =1:length(col_names)
            plot([NaN NaN], [NaN NaN], 'Color', currColors{j}, 'DisplayName', strcat('Session', col_names(j)))
        end
    end
    
    

    
    
    
    
    
    
    
    
    
    %% Process all datasets for current animal, both for HIT (1) and MISS (2) scenarios
    %   Processing is idential for both Hit and Miss scenarios, simply
    %   pulling from different datasets
%     for HITvMISS = 1:2
        HITvMISS = 0;
         HSS = 'HIT';
         % HSS2 = 'HIT';
        % if(HITvMISS==1) % Assign Hit or Miss trial analysis string (Hit String Select)
        %     HSS = 'HIT';
        % else
        %     HSS = 'MISS';
        % end

        % Iterate through each unique combination
        deltaOverall = 0;
        deltaInhibit = 0;
        deltaG = 0;
        dC = 0;
        correlationMatriciesNeuronCount = cell(size(uChanCurrs,1),1);
        allRasters = cell(size(uChanCurrs,1),1);
        correlationMatriciesNeuronCountNorm = cell(size(uChanCurrs,1),1);

        correlationMatricies = cell(size(uChanCurrs,1),1);
        correlationMatricies3D = cell(size(uChanCurrs,1),1);
        correlationMatriciesBase = cell(size(uChanCurrs,1),1);
        correlationMatricies3DBase = cell(size(uChanCurrs,1),1);
        diffMatricies = cell(size(uChanCurrs,1),1);
        diffMatriciesSTD = cell(size(uChanCurrs,1),1);
        diffMatriciesVals = cell(size(uChanCurrs,1),1);
        diffMatricies3D = cell(size(uChanCurrs,1),1);
        correlationMatriciesROI = cell(size(uChanCurrs,1),1);
        correlationMatriciesROI3D = cell(size(uChanCurrs,1),1);
        correlationDates = cell(size(uChanCurrs,1),1);
        overallBaseMean = cell(size(uChanCurrs,1),1);
        overallBaseMin = cell(size(uChanCurrs,1),1);
        overallBaseMax = cell(size(uChanCurrs,1),1);
        overallBaseSTD = cell(size(uChanCurrs,1),1);
        overallBaseMeanUpper10= cell(size(uChanCurrs,1),1);
        populationDistance = cell(size(uChanCurrs,1),1);
        enclosingRadius = cell(size(uChanCurrs,1),1);
        populationDistances = cell(size(uChanCurrs,1),1);
        densityDistance = cell(size(uChanCurrs,1),1);
        popDates =  cell(size(uChanCurrs,1),1);
        populationHistogram =  cell(size(uChanCurrs,1),1);
        
%         masterROISum = zeros(512,512,sharedDepth);
        
        
        
        for curPerm = 1:size(uChanCurrs,1)


            % Calculate electrode position
            curElecPos = [0,0,0];
            for k = 1:numel(contactsOrdered)
                if(contactsOrdered(k)==uChanCurrs(curPerm,1))
                    % Current channel matches stimulating channel,
                    % pull position and use for density analysis
                    curElecPos = elecPts(k,:);
                    break
                end
                if(uChanCurrs(curPerm,1)==29)
                    curElecPos = elecPts(3,:);
                    break
                end
                if(uChanCurrs(curPerm,1)==10)
                    curElecPos = elecPts(4,:);
                    break
                end
                if(uChanCurrs(curPerm,1)==27)
                    curElecPos = elecPts(5,:);
                    break
                end
            end



             % find population of sessions which had trials at the specified
             % parameter set
             [~, inds]=ismember(chanCurrPerms(:,1:2),uChanCurrs(curPerm,:),'rows');
             inds2 = find(inds);
             activeRows = chanCurrPermsOrig(inds2,:);




             %check if each value exists for hits (they might not)
             validFile = true(size(activeRows,1),1);
             for r = 1:size(activeRows,1) % For each active row in current set, pull values
                % Load baseline mean volume for registration

                % Open ROI file for each channel/current
                curpath = trueSourceFolders{activeRows(r,3)};
                if(mod(activeRows(r,3),2)==1)
                    HSS = 'HIT';
                else
                    HSS = 'MISS';
                end
                fileNameH = strcat(curpath,'\Ch',num2str(activeRows(r,1)),'-CUR',num2str(activeRows(r,2)),HSS,'\DeltaVol.mat'); 
                if(~isfile(fileNameH)) % check if older filename style is used or new format
                    fileNameH = strcat(curpath,'\Ch',num2str(activeRows(r,1)),'-CUR',num2str(activeRows(r,2)),'FREQ',num2str(activeRows(r,4)),'PULSE',num2str(activeRows(r,5)),HSS,'\DeltaVol.mat');
                end
                if(~isfile(fileNameH)) % check file existence
                    validFile(r)=0;
                end
             end
             activeRows = activeRows(validFile,:); % only process where there is data....  this is going to be hell for the miss and comparitive analysis



             % construct matrix to hold values from each entry for comparision
             sessData = zeros(size(activeRows,1),numel(vInds2D));% 2D comparision dataset
             sessData3D = zeros(size(activeRows,1),numel(vInds3D));% 3D comparision dataset
             sessDataBase = zeros(size(activeRows,1),numel(vInds2D));% 2D comparision dataset
             sessData3DBase = zeros(size(activeRows,1),numel(vInds3D));% 3D comparision dataset
             sessROI = zeros(size(activeRows,1),numel(vInds2D));% 2D ROI comparision dataset
             sessROI3D = zeros(size(activeRows,1),numel(vInds3D));% 3D ROI comparision dataset
             sessActive = zeros(size(activeRows,1),numel(vInds2D));% 2D ROI comparision dataset
             sessActive3D = zeros(size(activeRows,1),numel(vInds3D));% 3D ROI comparision dataset
             sessNeuronCount = zeros(size(activeRows,1),1); % Number of active neurons
			 sessNeuronCountNormalized = zeros(numel(activeRows,1),1); % Number of active neurons normalized
             HitDensity = zeros(size(activeRows,1),10);
             HitDistance = cell(size(activeRows,1),1);
             HitRadius = zeros(size(activeRows,1),1);
             poplationDistHIT= zeros(size(activeRows,1),10);
             poplationDistHITFine= zeros(size(activeRows,1),20);
             activeRowDates = zeros(size(activeRows,1),1);
             sessBaseMean = zeros(size(activeRows,1),1); 
             sessBaseSTD = zeros(size(activeRows,1),1); 
             sessBaseMin = zeros(size(activeRows,1),1);
             sessBaseMax = zeros(size(activeRows,1),1);
             sessBaseMeanUpper10 = zeros(size(activeRows,1),1); 
             sessRaster = cell(size(activeRows,1),1);


             %Add new analysis values here to collect across iterations of
             %combination... 

             for r = 1:size(activeRows,1) % For each active row in current set, pull values
                % Load baseline mean volume for registration
                curRowFile = activeRows(r,3);




                % Open ROI file for each channel/current

                if(mod(activeRows(r,3),2)==1)
                    HSS = 'HIT';
                else
                    HSS = 'MISS';
                end


                curpath = trueSourceFolders{activeRows(r,3)};
                fileNameP = strcat(curpath,'\Ch',num2str(activeRows(r,1)),'-CUR',num2str(activeRows(r,2)),HSS,'\PreStim.mat'); 
                fileNameH = strcat(curpath,'\Ch',num2str(activeRows(r,1)),'-CUR',num2str(activeRows(r,2)),HSS,'\DeltaVol.mat'); 
                fileNameROIH = strcat(curpath,'\Ch',num2str(activeRows(r,1)),'-CUR',num2str(activeRows(r,2)),HSS,'\ROI.mat'); 

                if(~isfile(fileNameH)) % check if older filename style is used or new format
                    fileNameP = strcat(curpath,'\Ch',num2str(activeRows(r,1)),'-CUR',num2str(activeRows(r,2)),'FREQ',num2str(activeRows(r,4)),'PULSE',num2str(activeRows(r,5)),HSS,'\PreStim.mat');
                    fileNameH = strcat(curpath,'\Ch',num2str(activeRows(r,1)),'-CUR',num2str(activeRows(r,2)),'FREQ',num2str(activeRows(r,4)),'PULSE',num2str(activeRows(r,5)),HSS,'\DeltaVol.mat');
                    fileNameROIH = strcat(curpath,'\Ch',num2str(activeRows(r,1)),'-CUR',num2str(activeRows(r,2)),'FREQ',num2str(activeRows(r,4)),'PULSE',num2str(activeRows(r,5)),HSS,'\ROI.mat');
%                     trialPath = strcat(pullPath,'\TRIALS2'); %construct folder for trial specific data, older format
                end
                hitData = load(fileNameH);
                roiData = load(fileNameROIH);
                preData = load(fileNameP);

                
                
                




                curDeltaMean = hitData.midDelta; % Midline Flourescence Volume
                curDeltaMean = curDeltaMean(:,:,1:sharedDepth);
                curBaseMean = preData.meanPreStim; % Midline Flourescence Volume
                curBaseMean = curBaseMean(:,:,1:sharedDepth);
                activeVolume = roiData.midActive;
                activeVolume = activeVolume(:,:,1:sharedDepth);

                % optionally pull whole possible ROIs
                midActive = roiData.midActive;
                possibleROI = false(size(midActive));
                regions = roiData.regions;
                
                if(roiMode==1)
                    sessRaster{r} = roiData.midRaster;
                    roiVolume = roiData.validROI; % Detected ROI volume
                    roiVolume = roiVolume(:,:,1:sharedDepth);
                else
                    curRaster = false(size(regions,1),1);
                    for i = 1:numel(regions)
                        % Pull sum of activations data from overall dataset
                        curPixels = regions(i).PixelList;
                        curMidVals = zeros(size(curPixels,1),1);
                        curVol = false(size(midActive));
                        for m = 1:size(curPixels,1)
                            curMidVals(m) = midActive(curPixels(m,2),curPixels(m,1),curPixels(m,3));
                            curVol(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
                        end   

                        if(sum(curMidVals)>=1)
                            possibleROI(curVol)=1;
                            curRaster(i)=1;
                        end
                    end
                    roiVolume = possibleROI(:,:,1:sharedDepth);
                    sessRaster{r} = curRaster;
                end
                
                

                % Normalize baseline data so data can be compared between days
    %             curBaseMean = curBaseMean/prestimMeans(activeRows(r,3)); % Normalize session baseline volumes to overall baseline mean 
                refBase = prestimMeanVols{activeRows(r,3)};
                curBaseMean = curBaseMean./refBase(:,:,1:sharedDepth);



                rasterMidActive = roiData.midRaster;


                if(isfield(roiData,'rasterPreStim'))
                    rasterPreActive = roiData.rasterPreStim;
                else
                    rasterPreActive = rasterMidActive;
                end
                
                
                % generate prestimulation active volume
                roiVolumePre = false(size(roiVolume));
                for i = 1:numel(regions)
                    if(rasterPreActive(i)==1)
                        curPixels = regions(i).PixelList;
                        for m = 1:size(curPixels,1)
                            roiVolumePre(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
                        end 
                    end
                end


                preMidOverlap = rasterPreActive;
                preMidOverlap(rasterMidActive==0)=0;


                if(isfield(roiData,'validROIInhibit'))
                    inhibitROI = roiData.validROIInhibit;
                else
                    inhibitROI = roiData.validROI;
                end

                % Collected delta of flourescence values for the inhbited
                % regions (Chong's analysis)
                inhibit1D = reshape(inhibitROI,1,[]);
                delta1D = reshape(hitData.midDelta,1,[]);
                inhibitDelta = delta1D(inhibit1D==1); 
                % Collect baseline and inhibited neural flourescence
                dC = dC+1;
                deltaInhibit(dC) = mean(inhibitDelta);
                deltaOverall(dC) = mean(delta1D);
                deltaG(dC) = curPerm;

                midActiveCount(dC) = sum(rasterMidActive);
                preActiveCount(dC) = sum(rasterPreActive);
                overlapCount(dC) = sum(preMidOverlap);


                % Flatten image
                curFlat = squeeze(mean(curDeltaMean,3)); % average for 2D analysis
                curFlatB = squeeze(mean(curBaseMean,3)); 
                roiFlat = squeeze(mean(roiVolume,3));
                activeFlat = squeeze(mean(activeVolume,3)); 

                %transform image - 2D
                shiftedImg = spatial_interp_patchwarp(curFlat, ptForms{curRowFile}, 'euclidean', 1:512, 1:512);
                shiftedImgB = spatial_interp_patchwarp(curFlatB, ptForms{curRowFile}, 'euclidean', 1:512, 1:512);
                shiftedROI = spatial_interp_patchwarp(roiFlat, ptForms{curRowFile}, 'euclidean', 1:512, 1:512);
                shiftedActive = spatial_interp_patchwarp(activeFlat, ptForms{curRowFile}, 'euclidean', 1:512, 1:512);



                shiftedVol = zeros(size(curDeltaMean));
                shiftedVolB = zeros(size(curBaseMean));
                shiftedRoiVol = zeros(size(roiVolume));
                shiftedActiveVol = zeros(size(activeVolume)); 
                shiftedRoiVolPre = zeros(size(roiVolumePre));


                for d = 1:sharedDepth
                    shiftedVol(:,:,d) = spatial_interp_patchwarp(curDeltaMean(:,:,d), ptForms{curRowFile}, 'euclidean', 1:512, 1:512);
                    shiftedVolB(:,:,d) = spatial_interp_patchwarp(curBaseMean(:,:,d), ptForms{curRowFile}, 'euclidean', 1:512, 1:512);
                    shiftedRoiVol(:,:,d) = spatial_interp_patchwarp(roiVolume(:,:,d), ptForms{curRowFile}, 'euclidean', 1:512, 1:512);
                    shiftedActiveVol(:,:,d) = spatial_interp_patchwarp(activeVolume(:,:,d), ptForms{curRowFile}, 'euclidean', 1:512, 1:512);
                    shiftedRoiVolPre(:,:,d) = spatial_interp_patchwarp(roiVolumePre(:,:,d), ptForms{curRowFile}, 'euclidean', 1:512, 1:512);
                end
                
                shiftedVol = shiftedVol(:,:,1:sharedDepth);
                shiftedVolB = shiftedVolB(:,:,1:sharedDepth);
                shiftedRoiVol = shiftedRoiVol(:,:,1:sharedDepth);
                shiftedActiveVol = shiftedActiveVol(:,:,1:sharedDepth);
                shiftedRoiVolPre = shiftedRoiVolPre(:,:,1:sharedDepth);
                

                % reduce dimensionality and extract valid indicies (aligned) for comparision
                % Raw Flouescence Image
                shiftedVector = reshape(shiftedImg,1,[]);
                validPixels = shiftedVector(vInds2D);
                sessData(r,:) = validPixels;

                shiftedVector3D = reshape(shiftedVol,1,[]);
                validVoxels = shiftedVector3D(vInds3D);
                sessData3D(r,:) = validVoxels;


                shiftedVector = reshape(shiftedImgB,1,[]);
                validPixels = shiftedVector(vInds2D);
                sessDataBase(r,:) = validPixels;

                shiftedVector3D = reshape(shiftedVolB,1,[]);
                validVoxels = shiftedVector3D(vInds3D);
                sessData3DBase(r,:) = validVoxels;

                % Calculate baseline metrics
                sessBaseMean(r) = mean(validVoxels);
                sessBaseSTD(r) = std(validVoxels);
                sessBaseMin(r) = min(validVoxels);
                sessBaseMax(r) = max(validVoxels);
                upper10 = sort(validVoxels);
                upper10 = upper10(round(0.9*numel(upper10)):end);
                sessBaseMeanUpper10(r) = mean(upper10);

                % ROI volume
                shiftedVector = reshape(shiftedROI,1,[]);
                validPixels = shiftedVector(vInds2D);
                sessROI(r,:) = validPixels;

                shiftedVector3D = reshape(shiftedRoiVol,1,[]);
                validVoxels = shiftedVector3D(vInds3D);
                sessROI3D(r,:) = validVoxels;
%                 masterROISum = masterROISum + shiftedRoiVol;
                
				% calculate the number of active regions based on the shared volume between sessions. Alternative calculates based on overall 
				shiftedRoiVolMasked = shiftedRoiVol;
				shiftedRoiVolMasked(overlay3D==0) = 0; %remove regions that are not shared between sessions
				curMaskedRegions = bwconncomp(shiftedRoiVolMasked);
				roiCount = curMaskedRegions.NumObjects; 
				
				
				shiftedBaselineVolMasked = shiftedRoiVolPre;
				shiftedBaselineVolMasked(overlay3D==0) = 0; %remove regions that are not shared between sessions
				curMaskedRegions = bwconncomp(shiftedBaselineVolMasked);
				baselineRoiCount = curMaskedRegions.NumObjects; 
				baselineRoiCount = max([baselineRoiCount, 1]);
				
				
				sessNeuronCount(r) = roiCount;
				sessNeuronCountNormalized(r) = roiCount/baselineRoiCount;

				
				
				
				
				


                % Active volume
                shiftedVector = reshape(shiftedActive,1,[]);
                validPixels = shiftedVector(vInds2D);
                sessActive(r,:) = validPixels;

                shiftedVector3D = reshape(shiftedActiveVol,1,[]);
                validVoxels = shiftedVector3D(vInds3D);
                sessActive3D(r,:) = validVoxels;



                % Perform density calculation for neural populations near 
                % For each detected neuron, calculate the distance from
                % the electrode that is was, in 3D space
                curSegData = bwareaopen(shiftedRoiVol,5); % remove small noise
                neuronRegs = regionprops(curSegData,'centroid');
                centroids = cat(1,neuronRegs.Centroid);
                if(numel(centroids)==0)
                    % If no centroids exist just skip this dataset
                    allDist3DHIT = 0;
                else
                    allDist3DHIT = zeros(size(centroids,1),1);
                    for k = 1:size(centroids,1) % Perform 2D & 3D distance calculations
                        allDist3DHIT(k) = sqrt(((centroids(k,1)-curElecPos(1))*2.187)^2 + ((centroids(k,2)-curElecPos(2))*2.187)^2  + ((centroids(k,3)-(curElecPos(3)/25))*vertStep)^2 );
                    end
                end



                % Calculate density for channel/current
                pop3DHIT = histcounts(allDist3DHIT,0:100:1000);
                pop3DHITfine = histcounts(allDist3DHIT,0:50:1000);
                temp = zeros(10,1);
                for k = 1:10
                    % calculate 3D density for area surveyed
                    temp(k) = pop3DHIT(k)/(((4/3)*pi()*(100*k)^2)-((4/3)*pi()*(100*(k-1))^2));
                end

                if(~isempty(centroids))
                    centroids(:,3) = centroids(:,3)*vertStep;
                    [minRad,~,~]=ExactMinBoundSphere3D(centroids);
                else
                    minRad=0;
                end
                HitRadius(r) = minRad;
                HitDistance{r} = allDist3DHIT;
                HitDensity(r,:) = temp;
                activeRowDates(r) = activeRows(r,3);
                poplationDistHIT(r,:) = pop3DHIT;
                poplationDistHITFine(r,:) =  pop3DHITfine;
             end


             
             
             
             
             populationHistogram{curPerm} = poplationDistHITFine;
             popDates{curPerm} = activeRowDates;
             populationDistance{curPerm} = poplationDistHIT;
             densityDistance{curPerm} =    HitDensity;
             enclosingRadius{curPerm} = HitRadius;
             populationDistances{curPerm} = HitDistance;
             % Analyze all session values and determine regions with low
             % intensity and little variability 



             % For current set of session values perform correlation analysis
             % pixel by pixel (or voxel by voxel)
             curCorrMatrix = zeros(size(activeRows,1),size(activeRows,1));
             curCorrMatrix3D = zeros(size(activeRows,1),size(activeRows,1));
             curCorrMatrixBase = zeros(size(activeRows,1),size(activeRows,1));
             curCorrMatrix3DBase = zeros(size(activeRows,1),size(activeRows,1));
             curDiffMatrix = zeros(size(activeRows,1),size(activeRows,1));
             curDiffMatrixSTD = zeros(size(activeRows,1),size(activeRows,1));
             curDiffMatrixVals = cell(size(activeRows,1),size(activeRows,1));

             curDiffMatrix3D = zeros(size(activeRows,1),size(activeRows,1));
             curROICorrMatrix = zeros(size(activeRows,1),size(activeRows,1));
             curROICorrMatrix3D = zeros(size(activeRows,1),size(activeRows,1));
             for i = 1:size(activeRows,1)
                 for j = 1:size(activeRows,1)
                     % Analyze correlation of flourescence


                     % Evaluate correlation, using only pairs where both are nonzero
                     chan1 = sessData(i,:)';
                     chan2 = sessData(j,:)';
                     base1 = sessDataBase(i,:)';
                     base2 = sessDataBase(j,:)';
                     ch1Valid = (chan1 > dfThresh) | (chan1 < -dfThresh);
                     ch2Valid = (chan2 > dfThresh) | (chan2 < -dfThresh);
                     idxBothNonzero = ch1Valid;
                     idxBothNonzero(ch2Valid==0)=0;


                     if(sum(idxBothNonzero)>0)
                         curCorrMatrix(i,j) = corr(chan1(idxBothNonzero),chan2(idxBothNonzero));
                         curDiffMatrix(i,j) = mean(chan1(idxBothNonzero)-chan2(idxBothNonzero)); 
                         curDiffMatrixSTD(i,j) = std(chan1(idxBothNonzero)-chan2(idxBothNonzero)); 
                         curDiffMatrixVals{i,j} = chan1(idxBothNonzero)-chan2(idxBothNonzero);
                     end



                     curCorrMatrixBase(i,j) = corr(base1,base2);



                     chan1 = sessData3D(i,:)';
                     chan2 = sessData3D(j,:)';
                     base1 = sessData3DBase(i,:)';
                     base2 = sessData3DBase(j,:)';
                     ch1Valid = (chan1 > dfThresh) | (chan1 < -dfThresh);
                     ch2Valid = (chan2 > dfThresh) | (chan2 < -dfThresh);
                     idxBothNonzero = ch1Valid & ch2Valid;
                     if(sum(idxBothNonzero)>0)
                         curCorrMatrix3D(i,j) = corr(chan1(idxBothNonzero),chan2(idxBothNonzero));
                         curDiffMatrix3D(i,j) = mean(chan1(idxBothNonzero)-chan2(idxBothNonzero));
                     end
                     curCorrMatrix3DBase(i,j) = corr(base1,base2);



                     % Analyze alignment of ROIs detected - since masks are
                     % binary we will use overlap % to determine fit
                     maskAll = sessROI(i,:);
                     maskAll(sessROI(j,:)==1)=1;
                     maskBoth = sessROI(i,:);
                     maskBoth(sessROI(j,:)==0)=0;
                     curROICorrMatrix(i,j) = maskBoth/maskAll;

                     maskAll = sessROI3D(i,:);
                     maskAll(sessROI3D(j,:)==1)=1;
                     maskBoth = sessROI3D(i,:);
                     maskBoth(sessROI3D(j,:)==0)=0;
                     curROICorrMatrix3D(i,j) = maskBoth/maskAll;
                 end
             end

             % Save resulting correlation matricies
             correlationDates{curPerm} = activeRows(:,3);
             correlationMatricies{curPerm} = curCorrMatrix;
             correlationMatricies3D{curPerm} = curCorrMatrix3D;
             correlationMatriciesBase{curPerm} = curCorrMatrixBase;
             correlationMatricies3DBase{curPerm} = curCorrMatrix3DBase;
             diffMatricies{curPerm} = curDiffMatrix;
             diffMatriciesSTD{curPerm} = curDiffMatrixSTD;
             diffMatriciesVals{curPerm} = curDiffMatrixVals;
             diffMatricies3D{curPerm} = curDiffMatrix3D;
             correlationMatriciesROI{curPerm} = curROICorrMatrix;
             correlationMatriciesROI3D{curPerm} = curROICorrMatrix3D;
             correlationMatriciesNeuronCount{curPerm} = sessNeuronCount; % track the number of active neruons
             correlationMatriciesNeuronCountNorm{curPerm} = sessNeuronCountNormalized; % track the number of active neruons
             allRasters{curPerm} = sessRaster;

             overallBaseMean{curPerm} = sessBaseMean;
             overallBaseMin{curPerm} = sessBaseMin;
             overallBaseMax{curPerm} = sessBaseMax;
             overallBaseSTD{curPerm} = sessBaseSTD;
             overallBaseMeanUpper10{curPerm} = sessBaseMeanUpper10;
        end

  
        
    
    
    
  




        %% Analyze the number of active neurons over time
        figure()
        hold on
        legText = cell(size(uChanCurrs,1),1);
        lgN = 0;
        uCurr = unique(uChanCurrs(:,2));
        uChan = unique(uChanCurrs(:,1));
        uThresh = unique(reshape(daysThresholds,1,[]));
        uThresh(uThresh==0)=[]; % remove thresholds at 0

        groupedCurrChans = NaN(numel(daysTrained),numel(uChan),numel(uCurr));
        groupedRasters = cell(numel(daysTrained),numel(uChan),numel(uCurr));
        currColors = {'r','m','b','y','g','k','#D95319','#7E2F8E','#77AC30','#4DBEEE','#EDB120','#7E2F8E','#FF00FF','#0072BD'}; % Hopefully these sets of colors and symbols will be enough
        chanSyms = {'o','+','square','x','diamond','*','^','<','>','|'};
        for cInd = 1:numel(uCurr)
            curCurr = uCurr(cInd); % process each current
            groupedCurrs = NaN(numel(daysTrained),size(uChanCurrs,1));

            for p = 1:size(uChanCurrs,1)
                if(uChanCurrs(p,2)==curCurr) % if valid current include in average

                    chanInd = find(uChan==uChanCurrs(p,1));
                    dateInds = correlationDates{p};
                    if(~isempty(dateInds))
                        curDays = daysTrained(dateInds);
                        neuCnt = correlationMatriciesNeuronCount{p};
                        rasters = allRasters{p};

                        p1 = plot(curDays,neuCnt,'--'); % plot the number of neurons active for this channel and current
                        p1.Marker = chanSyms{chanInd};
                        p1.Color = currColors{cInd};
                        groupedCurrs(dateInds,p)=neuCnt;
                        lgN = lgN + 1;
                        legText{lgN} = strcat('CH',num2str(uChanCurrs(p,1)),'--CURR',num2str(uChanCurrs(p,2)));
                        groupedCurrChans(dateInds,chanInd,cInd) = neuCnt;
                        groupedRasters(dateInds,chanInd,cInd) = rasters;
                    end
                end
            end    
            % plot grouped neural activation for this current
            grpTrace = mean(groupedCurrs,2,'omitnan');
            curTrace2 = grpTrace;
            curTrace2(isnan(curTrace2))=0;
            ctX = find(curTrace2);
            ctY = curTrace2(ctX);
            if(~isempty(ctY))
                p2 = plot(daysTrained(ctX),ctY,'linewidth',2);
                p2.Color = currColors{cInd};
                lgN = lgN + 1;
                legText{lgN} = strcat('Grouped CURR',num2str(curCurr));
            end
        end
        toKeep = false(numel(legText),1);
        for k = 1:numel(legText)
            if(~isempty(legText{k}))
                toKeep(k)=1;
            end
        end
        legText = legText(toKeep);
        legend(legText)
        title(strcat('Neural Activation Vs Time'))
        ylabel('Number of  Neurons')
        xlabel('Days of Training')

        
        
        
        
        
         %% Analyze the number of active neurons over Weeks
        figure()
        hold on
        legText = cell(size(uChanCurrs,1),1);
        lgN = 0;
        uCurr = unique(uChanCurrs(:,2));
        uChan = unique(uChanCurrs(:,1));
        uThresh = unique(reshape(daysThresholds,1,[]));
        uThresh(uThresh==0)=[]; % remove thresholds at 0

        weeksTrained = unique(ceil(daysTrained/7));
        chanSyms = {'o','+','square','x','diamond','*','^','<','>','|'};
        for cInd = 1:numel(uCurr)
            curCurr = uCurr(cInd); % process each current
            groupedCurrs = NaN(numel(weeksTrained),size(uChanCurrs,1));

            for p = 1:size(uChanCurrs,1)
                if(uChanCurrs(p,2)==curCurr) % if valid current include in average

                    chanInd = find(uChan==uChanCurrs(p,1));
                    dateInds = correlationDates{p};
                    curDays = ceil(daysTrained(dateInds)/7);
                    neuCnt = correlationMatriciesNeuronCount{p};

                    % check if multiple sessions occur within the same
                    % week, if so, merge and take average
                    uCurDays = unique(curDays);
                    uNeuCnt = zeros(size(uCurDays));
                    dateInds2 = zeros(size(uCurDays));
                    for uc = 1:numel(uCurDays)
                        ucInds = find(curDays==uCurDays(uc));
                        uNeuCnt(uc) = mean(neuCnt(ucInds));
                        dateInds2(uc) = find(weeksTrained == uCurDays(uc));
                    end
                    if(~isempty(uNeuCnt))      
                        p1 = plot(uCurDays,uNeuCnt,'--'); % plot the number of neurons active for this channel and current
                        p1.Marker = chanSyms{chanInd};
                        p1.Color = currColors{cInd};
                        groupedCurrs(dateInds2,p)=uNeuCnt;
                        lgN = lgN + 1;
                        legText{lgN} = strcat('CH',num2str(uChanCurrs(p,1)),'--CURR',num2str(uChanCurrs(p,2)));
                    end
                end
            end    
            % plot grouped neural activation for this current
            grpTrace = mean(groupedCurrs,2,'omitnan');
            curTrace2 = grpTrace;
            curTrace2(isnan(curTrace2))=0;
            ctX = find(curTrace2);
            ctY = curTrace2(ctX);
            if(~isempty(ctX))
                p2 = plot(weeksTrained(ctX),ctY,'linewidth',2);
                p2.Color = currColors{cInd};
                lgN = lgN + 1;
                legText{lgN} = strcat('Grouped CURR',num2str(curCurr));
            end
        end
        legend(legText)
        title('Neural Activation Vs Weeks')
        ylabel('Number of  Neurons')
        xlabel('Weeks of Training')
        
        
        
        
        % % plot the threshold of current animal over time
        % figure()
        % hold on
        % for ch = 1:numel(uChan)
        %     thresh = daysThresholdsOrig(ch,:);
        %     curDays = daysTrained/7;
        %     curDays(thresh==0)=[];
        %     thresh(thresh==0)=[];
        %     p1 = plot(curDays,thresh,'linewidth',2);
        %     p1.Marker = chanSyms{ch};
        %     p1.Color = [0 0 0];
        % end
        % ylim([0 10])
        % [hleg,~] = legend(num2str(uChan));
        % title(hleg,'Channel')
        % xlabel('Weeks of Training')
        % ylabel('Threshold (\muA)')
        % title('Detection Threshold over Time')
        
        
        
         %% Analyze the number of active neurons over time
        figure()
        hold on
        legText = cell(size(uChanCurrs,1),1);
        lgN = 0;
        uCurr = unique(uChanCurrs(:,2));
        uChan = unique(uChanCurrs(:,1));
        uThresh = unique(reshape(daysThresholds,1,[]));
        uThresh(uThresh==0)=[]; % remove thresholds at 0

        groupedCurrChansNorm = NaN(numel(daysTrained),numel(uChan),numel(uCurr));
        for cInd = 1:numel(uCurr)
            curCurr = uCurr(cInd); % process each current
            groupedCurrs = NaN(numel(daysTrained),size(uChanCurrs,1));

            for p = 1:size(uChanCurrs,1)
                if(uChanCurrs(p,2)==curCurr) % if valid current include in average

                    chanInd = find(uChan==uChanCurrs(p,1));
                    dateInds = correlationDates{p};
                    curDays = daysTrained(dateInds);
                    neuCnt = correlationMatriciesNeuronCountNorm{p};

                    p1 = plot(curDays,neuCnt,'--'); % plot the number of neurons active for this channel and current
                    p1.Marker = chanSyms{chanInd};
                    p1.Color = currColors{cInd};
                    groupedCurrs(dateInds,p)=neuCnt;
                    lgN = lgN + 1;
                    legText{lgN} = strcat('CH',num2str(uChanCurrs(p,1)),'--CURR',num2str(uChanCurrs(p,2)));
                    groupedCurrChansNorm(dateInds,chanInd,cInd) = neuCnt;
                end
            end    
            % plot grouped neural activation for this current
            grpTrace = mean(groupedCurrs,2,'omitnan');
            curTrace2 = grpTrace;
            curTrace2(isnan(curTrace2))=0;
            ctX = find(curTrace2);
            ctY = curTrace2(ctX);
            if(~isempty(ctX))
                p2 = plot(daysTrained(ctX),ctY,'linewidth',2);
                p2.Color = currColors{cInd};
                lgN = lgN + 1;
                legText{lgN} = strcat('Grouped CURR',num2str(curCurr));
            end
        end
        legend(legText)
        title('Normalized Neural Activation Vs Time')
        ylabel('Number of  Neurons')
        xlabel('Days of Training')
        
        
        
        
        
        
        %% Analyze the normalized number of active neurons over Weeks
        figure()
        hold on
        legText = cell(size(uChanCurrs,1),1);
        lgN = 0;
        uCurr = unique(uChanCurrs(:,2));
        uChan = unique(uChanCurrs(:,1));
        uThresh = unique(reshape(daysThresholds,1,[]));
        uThresh(uThresh==0)=[]; % remove thresholds at 0

        weeksTrained = unique(ceil(daysTrained/7));
        for cInd = 1:numel(uCurr)
            curCurr = uCurr(cInd); % process each current
            groupedCurrs = NaN(numel(weeksTrained),size(uChanCurrs,1));

            for p = 1:size(uChanCurrs,1)
                if(uChanCurrs(p,2)==curCurr) % if valid current include in average

                    chanInd = find(uChan==uChanCurrs(p,1));
                    dateInds = correlationDates{p};
                    curDays = ceil(daysTrained(dateInds)/7);
                    neuCnt = correlationMatriciesNeuronCountNorm{p};

                    % check if multiple sessions occur within the same
                    % week, if so, merge and take average
                    uCurDays = unique(curDays);
                    uNeuCnt = zeros(size(uCurDays));
                    dateInds2 = zeros(size(uCurDays));
                    for uc = 1:numel(uCurDays)
                        ucInds = find(curDays==uCurDays(uc));
                        uNeuCnt(uc) = mean(neuCnt(ucInds));
                        dateInds2(uc) = find(weeksTrained == uCurDays(uc));
                    end
                                        
                    p1 = plot(uCurDays,uNeuCnt,'--'); % plot the number of neurons active for this channel and current
                    p1.Marker = chanSyms{chanInd};
                    p1.Color = currColors{cInd};
                    groupedCurrs(dateInds2,p)=uNeuCnt;
                    lgN = lgN + 1;
                    legText{lgN} = strcat('CH',num2str(uChanCurrs(p,1)),'--CURR',num2str(uChanCurrs(p,2)));
                end
            end    
            % plot grouped neural activation for this current
            grpTrace = mean(groupedCurrs,2,'omitnan');
            curTrace2 = grpTrace;
            curTrace2(isnan(curTrace2))=0;
                ctX = find(curTrace2);
                ctY = curTrace2(ctX);
            if(~isempty(ctX))
                p2 = plot(weeksTrained(ctX),ctY,'linewidth',2);
                p2.Color = currColors{cInd};
                lgN = lgN + 1;
                legText{lgN} = strcat('Grouped CURR',num2str(curCurr));
            end
        end
        legend(legText)
        title('Normalized Neural Activation Vs Weeks')
        ylabel('Number of  Neurons')
        xlabel('Weeks of Training')
        
        

        
        
        %% Analyze the number of active neurons over time at the threshold of detection
        figure()
        hold on
        legText = cell(numel(uChan),1);
        lgN = 0;
        uCurr = unique(uChanCurrs(:,2));
        uChan = unique(uChanCurrs(:,1));

        for chInd = 1:numel(uChan)
            curThreshes = daysThresholdsOrig(chInd,:); % Thresholds for current channel
            curNeuonCnts = squeeze(groupedCurrChans(:,chInd,:));

            interpNeuCnts = NaN(numel(curThreshes),1);
            for d = 1:numel(curThreshes)
                % find interpolated neural activation (if possible)
                interpNeuCnts(d) = interp1(uCurr,squeeze(curNeuonCnts(d,:)),curThreshes(d));
            end

            dateInds = correlationDates{p};
            curDays = daysTrained;

            p1 = plot(curDays,interpNeuCnts,'-'); % plot the number of neurons active at threshold for this channel 
            p1.Color = currColors{chInd};
            lgN = lgN + 1;
            legText{lgN} = strcat('CH',num2str(uChan(chInd)));
        end
        
        legend(legText)
        title('Neural Activation at Detection Threshold Vs Time')
        ylabel('Number of  Neurons')
        xlabel('Days of Training')
        
        
        
        
        figure()
        hold on
        legText = cell(numel(uChan),1);
        lgN = 0;
        uCurr = unique(uChanCurrs(:,2));
        uChan = unique(uChanCurrs(:,1));

        for chInd = 1:numel(uChan)
            curThreshes = daysThresholdsOrig(chInd,:); % Thresholds for current channel
            curNeuonCnts = squeeze(groupedCurrChans(:,chInd,:));

            interpNeuCnts = NaN(numel(curThreshes),1);
            for d = 1:numel(curThreshes)
                % find interpolated neural activation (if possible)
                interpNeuCnts(d) = interp1(uCurr,curNeuonCnts(d,:),curThreshes(d));
            end

            % prep data for week structure rather than by day
            dateInds = correlationDates{p};
            curDays = ceil(daysTrained/7);
            
            % check if multiple sessions occur within the same
            % week, if so, merge and take average
            uCurDays = unique(curDays);
            interpNeuCnts2 = zeros(size(uCurDays));
            dateInds2 = zeros(size(uCurDays));
            for uc = 1:numel(uCurDays)
                ucInds = find(curDays==uCurDays(uc));
                interpNeuCnts2(uc) = mean(interpNeuCnts(ucInds));
                dateInds2(uc) = find(weeksTrained == uCurDays(uc));
            end
            

            p1 = plot(uCurDays,interpNeuCnts2,'-'); % plot the number of neurons active at threshold for this channel 
            p1.Color = currColors{chInd};
            lgN = lgN + 1;
            legText{lgN} = strcat('CH',num2str(uChan(chInd)));
        end
        
        legend(legText)
        title('Neural Activation at Detection Threshold Vs Weeks')
        ylabel('Number of  Neurons')
        xlabel('Weeks of Training')
        
        
        


        % Overall neural activation with respect to threshold
        figure()
        hold on
        groupedCurrChansThresh = NaN(numel(uThresh),numel(uChan),numel(uCurr));
        for cInd = 1:numel(uCurr)
            curCurr = uCurr(cInd); % process each current
            groupedCurrsThresh = NaN(numel(uThresh),size(uChanCurrs,1));
            for p = 1:size(uChanCurrs,1)
                if(uChanCurrs(p,2)==curCurr) % if valid current include in average

                    chanInd = find(uChan==uChanCurrs(p,1));
                    dateInds = correlationDates{p};

                    curThresholds = daysThresholds(chanInd,dateInds);
                    threshInds = zeros(numel(curThresholds),1);
                    for c = 1:numel(curThresholds)
                        threshInds(c) = find(uThresh == curThresholds(c));
                    end
                    utInds = unique(threshInds);
                    neuCnt = correlationMatriciesNeuronCount{p};
                    neuCnts2 = zeros(numel(utInds),1);
                    neuCnts2Num = zeros(numel(utInds),1);
                    indsInds = zeros(numel(threshInds),1); % I've gone too deep with indexing, send help
                    for c = 1:numel(curThresholds)
                        indsInds(c) = find(utInds == threshInds(c));
                    end

                    for c = 1:numel(curThresholds)
                        if(neuCnts2(indsInds(c))==0)
                            neuCnts2(indsInds(c)) = neuCnt(c);
                            neuCnts2Num(indsInds(c)) = 1;
                        else
                            neuCnts2(indsInds(c)) = neuCnts2(indsInds(c)) + neuCnt(c);
                            neuCnts2Num(indsInds(c)) = neuCnts2Num(indsInds(c))+1;
                        end
                    end
                    neuCnts2 = neuCnts2./neuCnts2Num; % average for multiple sessions at same threshold

                    p1 = plot(uThresh(utInds),neuCnts2,'--'); % plot the number of neurons active for this channel and current
                    p1.Marker = chanSyms{chanInd};
                    p1.Color = currColors{cInd};
                    groupedCurrsThresh(utInds,p)=neuCnts2;
                    lgN = lgN + 1;
                    legText{lgN} = strcat('CH',num2str(uChanCurrs(p,1)),'--CURR',num2str(uChanCurrs(p,2)));
                    groupedCurrChansThresh(utInds,chanInd,cInd) = neuCnts2;
                end
            end    
            % plot grouped neural activation for this current
            grpTrace = mean(groupedCurrsThresh,2,'omitnan');
            curTrace2 = grpTrace;
            curTrace2(isnan(curTrace2))=0;
            ctX = find(curTrace2);
            ctY = curTrace2(ctX);
            if(~isempty(ctX))
                p2 = plot(uThresh(ctX),ctY,'linewidth',2);
                p2.Color = currColors{cInd};
                lgN = lgN + 1;
                legText{lgN} = strcat('Grouped CURR',num2str(curCurr));
            end
        end
        legend(legText)
        title('Neural Activation Vs Threshold')
        ylabel('Number of  Neurons')
        xlabel('Threshold (\muA)')




        %% Plot baseline activity across time and threshold
        uCurr = unique(uChanCurrs(:,2));
        uChan = unique(uChanCurrs(:,1));
        uThresh = unique(reshape(daysThresholds,1,[]));
        uThresh(uThresh==0)=[]; % remove thresholds at 0

        groupedCurrChansBaseMean = NaN(numel(daysTrained),numel(uChan),numel(uCurr));
        groupedCurrChansBaseMin = NaN(numel(daysTrained),numel(uChan),numel(uCurr));
        groupedCurrChansBaseMax = NaN(numel(daysTrained),numel(uChan),numel(uCurr));
        groupedCurrChansBaseSTD = NaN(numel(daysTrained),numel(uChan),numel(uCurr));
        groupedCurrChansBaseMeanUpper10 = NaN(numel(daysTrained),numel(uChan),numel(uCurr));
        groupedCurrChansThreshBaseMean = NaN(numel(uThresh),numel(uChan),numel(uCurr));
        groupedCurrChansThreshBaseMin = NaN(numel(uThresh),numel(uChan),numel(uCurr));
        groupedCurrChansThreshBaseMax = NaN(numel(uThresh),numel(uChan),numel(uCurr));
        groupedCurrChansThreshBaseSTD = NaN(numel(uThresh),numel(uChan),numel(uCurr));
        groupedCurrChansThreshBaseMeanUpper10 = NaN(numel(uThresh),numel(uChan),numel(uCurr));




        % Baseline mean intensity - with respect to time
        curFig = figure();
        hold on
        legText = cell(size(uChanCurrs,1),1);
        lgN = 0;    

        for cInd = 1:numel(uCurr)
            curCurr = uCurr(cInd); % process each current
            groupedCurrsResponse = NaN(numel(daysTrained),size(uChanCurrs,1));

            for p = 1:size(uChanCurrs,1)
                if(uChanCurrs(p,2)==curCurr) % if valid current include in average

                    chanInd = find(uChan==uChanCurrs(p,1));
                    dateInds = correlationDates{p};
                    curDays = daysTrained(dateInds);
                    bMean = overallBaseMean{p};

                    p1 = plot(curDays,bMean,'--'); % plot the number of neurons active for this channel and current
                    p1.Marker = chanSyms{chanInd};
                    p1.Color = currColors{cInd};
                    groupedCurrsResponse(dateInds,p)=bMean;
                    lgN = lgN + 1;
                    legText{lgN} = strcat('CH',num2str(uChanCurrs(p,1)),'--CURR',num2str(uChanCurrs(p,2)));
                    groupedCurrChansBaseMean(dateInds,chanInd,cInd) = bMean; % This gets sent along for cross animal analysis
                end
            end    
            % plot grouped neural activation for this current
            grpTrace = mean(groupedCurrsResponse,2,'omitnan');
            curTrace2 = grpTrace;
            curTrace2(isnan(curTrace2))=0;
            ctX = find(curTrace2);
            ctY = curTrace2(ctX);
            if(~isempty(ctX))
                p2 = plot(daysTrained(ctX),ctY,'linewidth',2);
                p2.Color = currColors{cInd};
                lgN = lgN + 1;
                legText{lgN} = strcat('Grouped CURR',num2str(curCurr));
            end
        end
        legend(legText)
        title('Baseline Mean Flourescence Vs Time')
        ylabel('Baseline Mean Flourescence intensity')
        xlabel('Days of Training')
        close(curFig)

        % Baseline mean intensity - with respect to threshold
        curFig = figure();
        hold on
        for cInd = 1:numel(uCurr)
            curCurr = uCurr(cInd); % process each current
            groupedCurrsThresh = NaN(numel(uThresh),size(uChanCurrs,1));
            for p = 1:size(uChanCurrs,1)
                if(uChanCurrs(p,2)==curCurr) % if valid current include in average

                    chanInd = find(uChan==uChanCurrs(p,1));
                    dateInds = correlationDates{p};

                    curThresholds = daysThresholds(chanInd,dateInds);
                    threshInds = zeros(numel(curThresholds),1);
                    for c = 1:numel(curThresholds)
                        threshInds(c) = find(uThresh == curThresholds(c));
                    end
                    utInds = unique(threshInds);
                    baseMetric = overallBaseMean{p}; %%%%%%%%%%%%%%%
                    baseMetric2 = zeros(numel(utInds),1);
                    neuCnts2Num = zeros(numel(utInds),1);
                    indsInds = zeros(numel(threshInds),1); % I've gone too deep with indexing, send help
                    for c = 1:numel(curThresholds)
                        indsInds(c) = find(utInds == threshInds(c));
                    end

                    for c = 1:numel(curThresholds)
                        if(baseMetric2(indsInds(c))==0)
                            baseMetric2(indsInds(c)) = baseMetric(c);
                            neuCnts2Num(indsInds(c)) = 1;
                        else
                            baseMetric2(indsInds(c)) = baseMetric2(indsInds(c)) + baseMetric(c);
                            neuCnts2Num(indsInds(c)) = neuCnts2Num(indsInds(c))+1;
                        end
                    end
                    baseMetric2 = baseMetric2./neuCnts2Num; % average for multiple sessions at same threshold

                    p1 = plot(uThresh(utInds),baseMetric2,'--'); % plot the baseline metricfor this channel and current
                    p1.Marker = chanSyms{chanInd};
                    p1.Color = currColors{cInd};
                    groupedCurrsThresh(utInds,p)=baseMetric2;
                    lgN = lgN + 1;
                    legText{lgN} = strcat('CH',num2str(uChanCurrs(p,1)),'--CURR',num2str(uChanCurrs(p,2)));
                    groupedCurrChansThreshBaseMean(utInds,chanInd,cInd) = baseMetric2;  % This gets sent along for cross animal analysis
                end
            end    
            % plot grouped neural activation for this current
            grpTrace = mean(groupedCurrsThresh,2,'omitnan');
            curTrace2 = grpTrace;
            curTrace2(isnan(curTrace2))=0;
            ctX = find(curTrace2);
            ctY = curTrace2(ctX);
            if(~isempty(ctX))
                p2 = plot(uThresh(ctX),ctY,'linewidth',2);
                p2.Color = currColors{cInd};
                lgN = lgN + 1;
                legText{lgN} = strcat('Grouped CURR',num2str(curCurr));
            end
        end
        legend(legText)
        title('Baseline Mean Flourescence Vs Threshold')
        ylabel('Baseline Mean Flourescence intensity')
        xlabel('Threshold (\muA)')
        close(curFig)



         % Baseline Upper 10% mean intensity - with respect to time
        curFig = figure();
        hold on
        legText = cell(size(uChanCurrs,1),1);
        lgN = 0;    

        for cInd = 1:numel(uCurr)
            curCurr = uCurr(cInd); % process each current
            groupedCurrsResponse = NaN(numel(daysTrained),size(uChanCurrs,1));

            for p = 1:size(uChanCurrs,1)
                if(uChanCurrs(p,2)==curCurr) % if valid current include in average

                    chanInd = find(uChan==uChanCurrs(p,1));
                    dateInds = correlationDates{p};
                    curDays = daysTrained(dateInds);
                    bMean10 = overallBaseMeanUpper10{p};

                    p1 = plot(curDays,bMean10,'--'); % plot the number of neurons active for this channel and current
                    p1.Marker = chanSyms{chanInd};
                    p1.Color = currColors{cInd};
                    groupedCurrsResponse(dateInds,p)=bMean10;
                    lgN = lgN + 1;
                    legText{lgN} = strcat('CH',num2str(uChanCurrs(p,1)),'--CURR',num2str(uChanCurrs(p,2)));
                    groupedCurrChansBaseMeanUpper10(dateInds,chanInd,cInd) = bMean10; % This gets sent along for cross animal analysis
                end
            end    
            % plot grouped neural activation for this current
            grpTrace = mean(groupedCurrsResponse,2,'omitnan');
            curTrace2 = grpTrace;
            curTrace2(isnan(curTrace2))=0;
            ctX = find(curTrace2);
            ctY = curTrace2(ctX);
            if(~isempty(ctX))
                p2 = plot(daysTrained(ctX),ctY,'linewidth',2);
                p2.Color = currColors{cInd};
                lgN = lgN + 1;
                legText{lgN} = strcat('Grouped CURR',num2str(curCurr));
            end
        end
        legend(legText)
        title('Upper 10% Baseline Mean Flourescence Vs Time')
        ylabel('Upper 10% Baseline Mean Flourescence intensity')
        xlabel('Days of Training')
        close(curFig)

        % Baseline mean intensity - with respect to threshold
        curFig = figure();
        hold on
        for cInd = 1:numel(uCurr)
            curCurr = uCurr(cInd); % process each current
            groupedCurrsThresh = NaN(numel(uThresh),size(uChanCurrs,1));
            for p = 1:size(uChanCurrs,1)
                if(uChanCurrs(p,2)==curCurr) % if valid current include in average

                    chanInd = find(uChan==uChanCurrs(p,1));
                    dateInds = correlationDates{p};

                    curThresholds = daysThresholds(chanInd,dateInds);
                    threshInds = zeros(numel(curThresholds),1);
                    for c = 1:numel(curThresholds)
                        threshInds(c) = find(uThresh == curThresholds(c));
                    end
                    utInds = unique(threshInds);
                    baseMetric = overallBaseMeanUpper10{p}; %%%%%%%%%%%%%%%
                    baseMetric2 = zeros(numel(utInds),1);
                    neuCnts2Num = zeros(numel(utInds),1);
                    indsInds = zeros(numel(threshInds),1); % I've gone too deep with indexing, send help
                    for c = 1:numel(curThresholds)
                        indsInds(c) = find(utInds == threshInds(c));
                    end

                    for c = 1:numel(curThresholds)
                        if(baseMetric2(indsInds(c))==0)
                            baseMetric2(indsInds(c)) = baseMetric(c);
                            neuCnts2Num(indsInds(c)) = 1;
                        else
                            baseMetric2(indsInds(c)) = baseMetric2(indsInds(c)) + baseMetric(c);
                            neuCnts2Num(indsInds(c)) = neuCnts2Num(indsInds(c))+1;
                        end
                    end
                    baseMetric2 = baseMetric2./neuCnts2Num; % average for multiple sessions at same threshold

                    p1 = plot(uThresh(utInds),baseMetric2,'--'); % plot the baseline metricfor this channel and current
                    p1.Marker = chanSyms{chanInd};
                    p1.Color = currColors{cInd};
                    groupedCurrsThresh(utInds,p)=baseMetric2;
                    lgN = lgN + 1;
                    legText{lgN} = strcat('CH',num2str(uChanCurrs(p,1)),'--CURR',num2str(uChanCurrs(p,2)));
                    groupedCurrChansThreshBaseMeanUpper10(utInds,chanInd,cInd) = baseMetric2;  % This gets sent along for cross animal analysis
                end
            end    
            % plot grouped neural activation for this current
            grpTrace = mean(groupedCurrsThresh,2,'omitnan');
            curTrace2 = grpTrace;
            curTrace2(isnan(curTrace2))=0;
            ctX = find(curTrace2);
            ctY = curTrace2(ctX);
            if(~isempty(ctX))
                p2 = plot(uThresh(ctX),ctY,'linewidth',2);
                p2.Color = currColors{cInd};
                lgN = lgN + 1;
                legText{lgN} = strcat('Grouped CURR',num2str(curCurr));
            end
        end
        legend(legText)
        title('Upper 10% Baseline Mean Flourescence Vs Threshold')
        ylabel('Upper 10% Baseline Mean Flourescence intensity')
        xlabel('Threshold (\muA)')
        close(curFig)





        % Baseline max intensity - with respect to time
    %     figure()
    %     hold on
        legText = cell(size(uChanCurrs,1),1);
        lgN = 0;    

        for cInd = 1:numel(uCurr)
            curCurr = uCurr(cInd); % process each current
            groupedCurrsResponse = NaN(numel(daysTrained),size(uChanCurrs,1));

            for p = 1:size(uChanCurrs,1)
                if(uChanCurrs(p,2)==curCurr) % if valid current include in average

                    chanInd = find(uChan==uChanCurrs(p,1));
                    dateInds = correlationDates{p};
                    curDays = daysTrained(dateInds);
                    bMax = overallBaseMax{p};

    %                 p1 = plot(curDays,bMax,'--'); % plot the number of neurons active for this channel and current
    %                 p1.Marker = chanSyms{chanInd};
    %                 p1.Color = currColors{cInd};
                    groupedCurrsResponse(dateInds,p)=bMax;
                    lgN = lgN + 1;
                    legText{lgN} = strcat('CH',num2str(uChanCurrs(p,1)),'--CURR',num2str(uChanCurrs(p,2)));
                    groupedCurrChansBaseMax(dateInds,chanInd,cInd) = bMax; % This gets sent along for cross animal analysis
                end
            end    
    %         % plot grouped neural activation for this current
    %         grpTrace = mean(groupedCurrsResponse,2,'omitnan');
    %         curTrace2 = grpTrace;
    %         curTrace2(isnan(curTrace2))=0;
    %         ctX = find(curTrace2);
    %         ctY = curTrace2(ctX);
    %         p2 = plot(daysTrained(ctX),ctY,'linewidth',2);
    %         p2.Color = currColors{cInd};
    %         lgN = lgN + 1;
    %         legText{lgN} = strcat('Grouped CURR',num2str(curCurr));
        end
    %     legend(legText)
    %     title('Baseline Max Flourescence Vs Time')
    %     ylabel('Baseline Max Flourescence intensity')
    %     xlabel('Days of Training')


        % Baseline max intensity - with respect to threshold
    %     figure()
    %     hold on
        for cInd = 1:numel(uCurr)
            curCurr = uCurr(cInd); % process each current
            groupedCurrsThresh = NaN(numel(uThresh),size(uChanCurrs,1));
            for p = 1:size(uChanCurrs,1)
                if(uChanCurrs(p,2)==curCurr) % if valid current include in average

                    chanInd = find(uChan==uChanCurrs(p,1));
                    dateInds = correlationDates{p};

                    curThresholds = daysThresholds(chanInd,dateInds);
                    threshInds = zeros(numel(curThresholds),1);
                    for c = 1:numel(curThresholds)
                        threshInds(c) = find(uThresh == curThresholds(c));
                    end
                    utInds = unique(threshInds);
                    baseMetric = overallBaseMax{p}; %%%%%%%%%%%%%%%
                    baseMetric2 = zeros(numel(utInds),1);
                    neuCnts2Num = zeros(numel(utInds),1);
                    indsInds = zeros(numel(threshInds),1); % I've gone too deep with indexing, send help
                    for c = 1:numel(curThresholds)
                        indsInds(c) = find(utInds == threshInds(c));
                    end

                    for c = 1:numel(curThresholds)
                        if(baseMetric2(indsInds(c))==0)
                            baseMetric2(indsInds(c)) = baseMetric(c);
                            neuCnts2Num(indsInds(c)) = 1;
                        else
                            baseMetric2(indsInds(c)) = baseMetric2(indsInds(c)) + baseMetric(c);
                            neuCnts2Num(indsInds(c)) = neuCnts2Num(indsInds(c))+1;
                        end
                    end
                    baseMetric2 = baseMetric2./neuCnts2Num; % average for multiple sessions at same threshold

    %                 p1 = plot(uThresh(utInds),baseMetric2,'--'); % plot the baseline metricfor this channel and current
    %                 p1.Marker = chanSyms{chanInd};
    %                 p1.Color = currColors{cInd};
                    groupedCurrsThresh(utInds,p)=baseMetric2;
                    lgN = lgN + 1;
                    legText{lgN} = strcat('CH',num2str(uChanCurrs(p,1)),'--CURR',num2str(uChanCurrs(p,2)));
                    groupedCurrChansThreshBaseMax(utInds,chanInd,cInd) = baseMetric2;  % This gets sent along for cross animal analysis
                end
            end    
    %         % plot grouped neural activation for this current
    %         grpTrace = mean(groupedCurrsThresh,2,'omitnan');
    %         curTrace2 = grpTrace;
    %         curTrace2(isnan(curTrace2))=0;
    %         ctX = find(curTrace2);
    %         ctY = curTrace2(ctX);
    %         p2 = plot(uThresh(ctX),ctY,'linewidth',2);
    %         p2.Color = currColors{cInd};
    %         lgN = lgN + 1;
    %         legText{lgN} = strcat('Grouped CURR',num2str(curCurr));
        end
    %     legend(legText)
    %     title('Baseline Max Flourescence Vs Time')
    %     ylabel('Baseline Max Flourescence intensity')
    %     xlabel('Threshold (\muA)')



        % Baseline min intensity - with respect to time
    %     figure()
    %     hold on
        legText = cell(size(uChanCurrs,1),1);
        lgN = 0;    

        for cInd = 1:numel(uCurr)
            curCurr = uCurr(cInd); % process each current
            groupedCurrsResponse = NaN(numel(daysTrained),size(uChanCurrs,1));

            for p = 1:size(uChanCurrs,1)
                if(uChanCurrs(p,2)==curCurr) % if valid current include in average

                    chanInd = find(uChan==uChanCurrs(p,1));
                    dateInds = correlationDates{p};
                    curDays = daysTrained(dateInds);
                    bMin = overallBaseMin{p};

    %                 p1 = plot(curDays,bMin,'--'); % plot the number of neurons active for this channel and current
    %                 p1.Marker = chanSyms{chanInd};
    %                 p1.Color = currColors{cInd};
                    groupedCurrsResponse(dateInds,p)=bMin;
                    lgN = lgN + 1;
                    legText{lgN} = strcat('CH',num2str(uChanCurrs(p,1)),'--CURR',num2str(uChanCurrs(p,2)));
                    groupedCurrChansBaseMin(dateInds,chanInd,cInd) = bMin; % This gets sent along for cross animal analysis
                end
            end    
    %         % plot grouped neural activation for this current
    %         grpTrace = mean(groupedCurrsResponse,2,'omitnan');
    %         curTrace2 = grpTrace;
    %         curTrace2(isnan(curTrace2))=0;
    %         ctX = find(curTrace2);
    %         ctY = curTrace2(ctX);
    %         p2 = plot(daysTrained(ctX),ctY,'linewidth',2);
    %         p2.Color = currColors{cInd};
    %         lgN = lgN + 1;
    %         legText{lgN} = strcat('Grouped CURR',num2str(curCurr));
        end
    %     legend(legText)
    %     title('Baseline Min Flourescence Vs Time')
    %     ylabel('Baseline Min Flourescence intensity')
    %     xlabel('Days of Training')


        % Baseline min intensity - with respect to threshold
    %     figure()
    %     hold on
        for cInd = 1:numel(uCurr)
            curCurr = uCurr(cInd); % process each current
            groupedCurrsThresh = NaN(numel(uThresh),size(uChanCurrs,1));
            for p = 1:size(uChanCurrs,1)
                if(uChanCurrs(p,2)==curCurr) % if valid current include in average

                    chanInd = find(uChan==uChanCurrs(p,1));
                    dateInds = correlationDates{p};

                    curThresholds = daysThresholds(chanInd,dateInds);
                    threshInds = zeros(numel(curThresholds),1);
                    for c = 1:numel(curThresholds)
                        threshInds(c) = find(uThresh == curThresholds(c));
                    end
                    utInds = unique(threshInds);
                    baseMetric = overallBaseMin{p}; %%%%%%%%%%%%%%%
                    baseMetric2 = zeros(numel(utInds),1);
                    neuCnts2Num = zeros(numel(utInds),1);
                    indsInds = zeros(numel(threshInds),1); % I've gone too deep with indexing, send help
                    for c = 1:numel(curThresholds)
                        indsInds(c) = find(utInds == threshInds(c));
                    end

                    for c = 1:numel(curThresholds)
                        if(baseMetric2(indsInds(c))==0)
                            baseMetric2(indsInds(c)) = baseMetric(c);
                            neuCnts2Num(indsInds(c)) = 1;
                        else
                            baseMetric2(indsInds(c)) = baseMetric2(indsInds(c)) + baseMetric(c);
                            neuCnts2Num(indsInds(c)) = neuCnts2Num(indsInds(c))+1;
                        end
                    end
                    baseMetric2 = baseMetric2./neuCnts2Num; % average for multiple sessions at same threshold
    % 
    %                 p1 = plot(uThresh(utInds),baseMetric2,'--'); % plot the baseline metricfor this channel and current
    %                 p1.Marker = chanSyms{chanInd};
    %                 p1.Color = currColors{cInd};
                    groupedCurrsThresh(utInds,p)=baseMetric2;
                    lgN = lgN + 1;
                    legText{lgN} = strcat('CH',num2str(uChanCurrs(p,1)),'--CURR',num2str(uChanCurrs(p,2)));
                    groupedCurrChansThreshBaseMin(utInds,chanInd,cInd) = baseMetric2;  % This gets sent along for cross animal analysis
                end
            end    
    %         % plot grouped neural activation for this current
    %         grpTrace = mean(groupedCurrsThresh,2,'omitnan');
    %         curTrace2 = grpTrace;
    %         curTrace2(isnan(curTrace2))=0;
    %         ctX = find(curTrace2);
    %         ctY = curTrace2(ctX);
    %         p2 = plot(uThresh(ctX),ctY,'linewidth',2);
    %         p2.Color = currColors{cInd};
    %         lgN = lgN + 1;
    %         legText{lgN} = strcat('Grouped CURR',num2str(curCurr));
        end
    %     legend(legText)
    %     title('Min Baseline Flourescence Vs Threshold')
    %     ylabel('Baseline Min Flourescence intensity')
    %     xlabel('Threshold (\muA)')




         % Baseline STD intensity - with respect to time
        curFig = figure();
        hold on
        legText = cell(size(uChanCurrs,1),1);
        lgN = 0;    
        for cInd = 1:numel(uCurr)
            curCurr = uCurr(cInd); % process each current
            groupedCurrsResponse = NaN(numel(daysTrained),size(uChanCurrs,1));
            for p = 1:size(uChanCurrs,1)
                if(uChanCurrs(p,2)==curCurr) % if valid current include in average

                    chanInd = find(uChan==uChanCurrs(p,1));
                    dateInds = correlationDates{p};
                    curDays = daysTrained(dateInds);
                    bSTD = overallBaseSTD{p};

                    p1 = plot(curDays,bSTD,'--'); % plot the number of neurons active for this channel and current
                    p1.Marker = chanSyms{chanInd};
                    p1.Color = currColors{cInd};
                    groupedCurrsResponse(dateInds,p)=bSTD;
                    lgN = lgN + 1;
                    legText{lgN} = strcat('CH',num2str(uChanCurrs(p,1)),'--CURR',num2str(uChanCurrs(p,2)));
                    groupedCurrChansBaseSTD(dateInds,chanInd,cInd) = bSTD; % This gets sent along for cross animal analysis
                end
            end    
            % plot grouped neural activation for this current
            grpTrace = mean(groupedCurrsResponse,2,'omitnan');
            curTrace2 = grpTrace;
            curTrace2(isnan(curTrace2))=0;
            ctX = find(curTrace2);
            ctY = curTrace2(ctX);
            if(~isempty(ctX))
                p2 = plot(daysTrained(ctX),ctY,'linewidth',2);
                p2.Color = currColors{cInd};
                lgN = lgN + 1;
                legText{lgN} = strcat('Grouped CURR',num2str(curCurr));
            end
        end
        legend(legText)
        title('Baseline STD Flourescence Vs Time')
        ylabel('Baseline STD Flourescence intensity')
        xlabel('Days of Training')
        close(curFig)

        % Baseline min intensity - with respect to threshold
        curFig = figure();
        hold on
        for cInd = 1:numel(uCurr)
            curCurr = uCurr(cInd); % process each current
            groupedCurrsThresh = NaN(numel(uThresh),size(uChanCurrs,1));
            for p = 1:size(uChanCurrs,1)
                if(uChanCurrs(p,2)==curCurr) % if valid current include in average
                    chanInd = find(uChan==uChanCurrs(p,1));
                    dateInds = correlationDates{p};
                    curThresholds = daysThresholds(chanInd,dateInds);
                    threshInds = zeros(numel(curThresholds),1);
                    for c = 1:numel(curThresholds)
                        threshInds(c) = find(uThresh == curThresholds(c));
                    end
                    utInds = unique(threshInds);
                    baseMetric = overallBaseSTD{p}; %%%%%%%%%%%%%%%
                    baseMetric2 = zeros(numel(utInds),1);
                    neuCnts2Num = zeros(numel(utInds),1);
                    indsInds = zeros(numel(threshInds),1); % I've gone too deep with indexing, send help
                    for c = 1:numel(curThresholds)
                        indsInds(c) = find(utInds == threshInds(c));
                    end

                    for c = 1:numel(curThresholds)
                        if(baseMetric2(indsInds(c))==0)
                            baseMetric2(indsInds(c)) = baseMetric(c);
                            neuCnts2Num(indsInds(c)) = 1;
                        else
                            baseMetric2(indsInds(c)) = baseMetric2(indsInds(c)) + baseMetric(c);
                            neuCnts2Num(indsInds(c)) = neuCnts2Num(indsInds(c))+1;
                        end
                    end
                    baseMetric2 = baseMetric2./neuCnts2Num; % average for multiple sessions at same threshold
                    p1 = plot(uThresh(utInds),baseMetric2,'--'); % plot the baseline metricfor this channel and current
                    p1.Marker = chanSyms{chanInd};
                    p1.Color = currColors{cInd};
                    groupedCurrsThresh(utInds,p)=baseMetric2;
                    lgN = lgN + 1;
                    legText{lgN} = strcat('CH',num2str(uChanCurrs(p,1)),'--CURR',num2str(uChanCurrs(p,2)));
                    groupedCurrChansThreshBaseSTD(utInds,chanInd,cInd) = baseMetric2;  % This gets sent along for cross animal analysis
                end
            end    
            % plot grouped neural activation for this current
            grpTrace = mean(groupedCurrsThresh,2,'omitnan');
            curTrace2 = grpTrace;
            curTrace2(isnan(curTrace2))=0;
            ctX = find(curTrace2);
            ctY = curTrace2(ctX);
            if(~isempty(ctX))
                p2 = plot(uThresh(ctX),ctY,'linewidth',2);
                p2.Color = currColors{cInd};
                lgN = lgN + 1;
                legText{lgN} = strcat('Grouped CURR',num2str(curCurr));
            end
        end
        legend(legText)
        title('Baseline STD Flourescence Vs Threshold')
        ylabel('Baseline STD Flourescence intensity')
        xlabel('Threshold (\muA)')
        close(curFig)



        %% Activation trends with regards to days of training
        uCurr = unique(uChanCurrs(:,2));
        uChan = unique(uChanCurrs(:,1));

        allDensity = cell(numel(uCurr),numel(daysTrained),10);
        allPopulation = cell(numel(uCurr),numel(daysTrained),10);
        for chInd = 1:numel(uChan)
            curCh = uChan(chInd);

            figure()
            hold on
            vCurr = false(numel(uCurr),1); 
            for cInd = 1:numel(uCurr)
                curCurr = uCurr(cInd); % process each current
                for p = 1:size(uChanCurrs,1)
                    if(uChanCurrs(p,2)==curCurr && uChanCurrs(p,1)==curCh ) % if valid current include in average
                        vCurr(cInd) = 1;
                        curRadii = enclosingRadius{p};
                        nDates = popDates{p};
                        curDays = daysTrained(nDates);  
                        p2 = plot(curDays,curRadii, '-o','linewidth',2);
                        p2.Color = currColors{cInd};
                    end
                end
            end
            legend(num2str(uCurr(vCurr)))
            title(strcat('Radius of Activation Vs Time  :  Ch',num2str(curCh), {'  '},HSS))
            ylabel('Activation Sphere Radius (\mum)')
            xlabel('Days of Training')




            for cInd = 1:numel(uCurr)
                curCurr = uCurr(cInd); % process each current
                for p = 1:size(uChanCurrs,1)
                    if(uChanCurrs(p,2)==curCurr && uChanCurrs(p,1)==curCh ) % if valid current include in average


                        curFig = figure();
                        hold on
                        curHistograms = populationDistance{p};
                        nDates = popDates{p};
                        curDays = daysTrained(nDates);
                        for d = 1:numel(nDates)
                            h = histogram('BinEdges',0:100:1000,'BinCounts',curHistograms(d,:)+(0.01*d),'DisplayStyle','stairs');
                            set(h,'LineWidth',3)
                        end
                        leg = legend(num2str(curDays'));
                        title(leg,'Days of Training')
                        title(strcat('Neural Activation Histogram Vs Time  :  Ch',num2str(curCh),'Curr',num2str(curCurr)))
                        ylabel('Neural Activation Counts')
                        xlabel('Days of Training')
                        close(curFig)

                        curFig = figure();
                        hold on
                        curDistances = populationDistances{p};
                        nDates = popDates{p};
                        curDays = daysTrained(nDates);  
                        g = [];
                        allDistances = [];
                        for d = 1:numel(nDates)
                            g = [g; d*ones(numel(curDistances{d}),1)];
                            allDistances = [allDistances; curDistances{d}];
                        end
                        boxplot(allDistances,g)
%                         title(leg,'Days of Training')
                        title(strcat('Neural Activation Distance Vs Time  :  Ch',num2str(curCh),'Curr',num2str(curCurr)))
                        ylabel('Activation Distance (\mum)')
                        xticklabels(num2str(curDays'))
                        xlabel('Days of Training')
                        close(curFig)

                        curFig = figure();
                        hold on
                        nDense = densityDistance{p}; 
                        nDates = popDates{p};
                        curDays = daysTrained(nDates);
                        for di = 1:10
                            plot(curDays,nDense(:,di),'-o','LineWidth',2)

                            for d = 1:numel(nDates)
                                cd = find(daysTrained == daysTrained(nDates(d)));
                                if(isempty(allDensity{cInd,cd,di}))
                                    allDensity{cInd,cd,di} = nDense(d,di);
                                else
                                    allDensity{cInd,cd,di} = [allDensity{cInd,cd,di}; nDense(d,di)];
                                end
                            end
                        end
                        legend({'0-100','100-200','200-300','300-400','400-500','500-600','600-700','700-800','800-900','900-1000'})
                        title(strcat('Neural Activation Density Vs Time  :  Ch',num2str(curCh),'Curr',num2str(curCurr)))
                        ylabel('Neural Activation Density (neurons/\mum^3)')
                        xlabel('Days of Training')
                        close(curFig)
                    end
                end 
            end
        end


        for cInd = 1:numel(uCurr)
            % Overall current response is not looking to clear, lots of popping
            % between dates ( force to weeks and rerun)
            curFig = figure();
            hold on
            for di = 1:10
                curDat = allDensity(cInd,:,di);
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

                plot(curDays,curMeans,'-o','LineWidth',2)
            end
            curCurr = uCurr(cInd);
            legend({'0-100','100-200','200-300','300-400','400-500','500-600','600-700','700-800','800-900','900-1000'})
            title(strcat('Overall Neural Activation Density Vs Time  :  Curr',num2str(curCurr)))
            ylabel('Neural Activation Density (neurons/\mum^3)')
            xlabel('Days of Training') 
            close(curFig)
        end



            % Overall Activation Distances per current
            for cInd = 1:numel(uCurr)
                curCurr = uCurr(cInd); % process each current
                curFig = figure();
                hold on
                g = [];
                allDistances = [];
                for p = 1:size(uChanCurrs,1)
                    if(uChanCurrs(p,2)==curCurr) % if valid current include in average

                        curDistances = populationDistances{p};
                        nDates = popDates{p};

                        for d = 1:numel(nDates)
                            g = [g; d*ones(numel(curDistances{d}),1)];
                            allDistances = [allDistances; curDistances{d}];
                        end
                    end
                end

                boxplot(allDistances,g)
                title(strcat('Neural Activation Distance Vs Time  : Curr',num2str(curCurr)))
                ylabel('Activation Distance (\mum)')

                curDays = daysTrained(unique(g));  
                xticklabels(num2str(curDays'))
                xlabel('Days of Training')
                close(curFig)
            end


            % Overall radius of activation vs current
            figure()
            hold on
            for cInd = 1:numel(uCurr)
                curCurr = uCurr(cInd); % process each current

                possibleRadii = NaN(numel(daysTrained),size(uChanCurrs,1));
                for p = 1:size(uChanCurrs,1)
                    if(uChanCurrs(p,2)==curCurr) % if valid current include in average
                        curRadii = enclosingRadius{p};
                        nDates = popDates{p};                    
                        possibleRadii(nDates,p) = curRadii;
                    end
                end
                aveRadii = mean(possibleRadii,2,'omitnan');
                p2 = plot(daysTrained,aveRadii, '-o','linewidth',2);
                p2.Color = currColors{cInd};
            end
            legend(num2str(uCurr))
            title('Overall Radius of Activation Vs Time')
            ylabel('Activation Sphere Radius (\mum)')
            xlabel('Days of Training')



        %% Activation trends with regards to days of detection threshold
        uCurr = unique(uChanCurrs(:,2));
        uChan = unique(uChanCurrs(:,1));


        uThresholds = unique(unique(daysThresholds));
        uThresholds(uThresholds==0)=[];
        allDensityThresh = NaN(numel(uChan),numel(uCurr),numel(uThresholds),10);
        allDensityThreshN = zeros(numel(uChan),numel(uCurr),numel(uThresholds),10);
        allRadiiThresh = NaN(numel(uChan),numel(uCurr),numel(uThresholds));


        for chInd = 1:numel(uChan)
            curCh = uChan(chInd);

            curFig = figure();
            hold on
            vCurr = false(numel(uCurr),1); 
            for cInd = 1:numel(uCurr)
                curCurr = uCurr(cInd); % process each current
                for p = 1:size(uChanCurrs,1)
                    if(uChanCurrs(p,2)==curCurr && uChanCurrs(p,1)==curCh ) % if valid current include in average
                        vCurr(cInd) = 1;
                        curRadii = enclosingRadius{p};
                        nDates = popDates{p};


                        % Covert data from date ordered to threshold grouped
                        curThresholds = daysThresholds(chInd,nDates);
                        threshInds = zeros(numel(curThresholds),1);
                        for c = 1:numel(curThresholds)
                            threshInds(c) = find(uThresh == curThresholds(c));
                        end
                        utInds = unique(threshInds);
                        var2 = zeros(numel(utInds),1);
                        var2Num = zeros(numel(utInds),1);
                        indsInds = zeros(numel(threshInds),1);
                        for c = 1:numel(curThresholds)
                            indsInds(c) = find(utInds == threshInds(c));
                        end
                        for c = 1:numel(curThresholds)
                            if(var2(indsInds(c))==0)
                                var2(indsInds(c)) = curRadii(c);
                                var2Num(indsInds(c)) = 1;
                            else
                                var2(indsInds(c)) = var2(indsInds(c)) + curRadii(c);
                                var2Num(indsInds(c)) = var2Num(indsInds(c))+1;
                            end
                        end
                        curRadii = var2./var2Num; % average for multiple sessions at same threshold
                        allRadiiThresh(chInd,cInd,utInds) = curRadii;


                        curThresh = uThresh(utInds);  
                        p2 = plot(curThresh,curRadii, '-o','linewidth',2);
                        p2.Color = currColors{cInd};
                    end
                end
            end
            legend(num2str(uCurr(vCurr)))
            title(strcat('Radius of Activation Vs Threshold  :  Ch',num2str(curCh)))
            ylabel('Activation Sphere Radius (\mum)')
            xlabel('Detection Threshold (\muA')
            close(curFig)



            for cInd = 1:numel(uCurr)
                curCurr = uCurr(cInd); % process each current
                for p = 1:size(uChanCurrs,1)
                    if(uChanCurrs(p,2)==curCurr && uChanCurrs(p,1)==curCh ) % if valid current include in average

                        nDates = popDates{p};
                        curHistograms = populationDistance{p};
                        curDistances = populationDistances{p};
                        nDense = densityDistance{p}; 
                        % Covert data from date ordered to threshold grouped
                        curThresholds = daysThresholds(chInd,nDates);
                        threshInds = zeros(numel(curThresholds),1);
                        for c = 1:numel(curThresholds)
                            threshInds(c) = find(uThresh == curThresholds(c));
                        end
                        utInds = unique(threshInds);
                        var2 = zeros(numel(utInds),10);
                        var4 = zeros(numel(utInds),10);

                        var2Num = zeros(numel(utInds),1);
                        indsInds = zeros(numel(threshInds),1);
                        for c = 1:numel(curThresholds)
                            indsInds(c) = find(utInds == threshInds(c));
                        end
                        for c = 1:numel(curThresholds)
                            if(var2(indsInds(c))==0)
                                var2(indsInds(c),:) = curHistograms(c,:);
                                var4(indsInds(c),:) = nDense(c,:);
                                var2Num(indsInds(c)) = 1;

                            else
                                var2(indsInds(c),:) = var2(indsInds(c),:) + curHistograms(c,:);
                                var4(indsInds(c),:) = var4(indsInds(c),:) + nDense(c,:);
                                var2Num(indsInds(c)) = var2Num(indsInds(c))+1;
                            end
                        end
                        curHistograms = var2./var2Num; % average for multiple sessions at same threshold
                        nDense = var4./var2Num; 
                        curThresh = uThresh(utInds);


                        curFig = figure();
                        hold on
                        for d = 1:numel(utInds)
                            h = histogram('BinEdges',0:100:1000,'BinCounts',curHistograms(d,:)+(0.01*d),'DisplayStyle','stairs');
                            set(h,'LineWidth',3)
                        end
                        leg = legend(num2str(curThresh'));
                        title(leg,'Detection Threshold (\muA)')
                        title(strcat('Neural Activation Histogram Vs Threshold  :  Ch',num2str(curCh),'Curr',num2str(curCurr)))
                        ylabel('Neural Activation Counts')
                        xlabel('Distance from Contact Site (\mum)')
                        close(curFig)


                        curFig = figure();
                        hold on
                        for di = 1:10
                            plot(curThresh,nDense(:,di),'-o','LineWidth',2)

                            for d = 1:numel(curThresh)

                                cd = find(uThresholds == curThresh(d));
                                if(isnan(allDensityThresh(chInd,cInd,cd,di)))
                                    allDensityThresh(chInd,cInd,cd,di) = nDense(d,di);
                                    allDensityThreshN(chInd,cInd,cd,di) = 1;
                                else
                                    allDensityThresh(chInd,cInd,cd,di) = allDensityThresh(chInd,cInd,cd,di)+ nDense(d,di);
                                    allDensityThreshN(chInd,cInd,cd,di) = allDensityThreshN(chInd,cInd,cd,di) + 1;
                                end
                            end
                        end
                        legend({'0-100','100-200','200-300','300-400','400-500','500-600','600-700','700-800','800-900','900-1000'})
                        title(strcat('Neural Activation Density Vs Threshold  :  Ch',num2str(curCh),'Curr',num2str(curCurr)))
                        ylabel('Neural Activation Density (neurons/\mum^3)')
                        xlabel('Distance from Contact Site (\mum)')
                        close(curFig)



                        curFig = figure();
                        hold on
                        curThresholds = daysThresholds(chInd,nDates);
                        uT = unique(curThresholds);
                        g = [];
                        allDistances = [];
                        for d = 1:numel(curThresholds)
                            g = [g; curThresholds(d)*ones(numel(curDistances{d}),1)];
                            allDistances = [allDistances; curDistances{d}];
                        end
                        boxplot(allDistances,g)
                        title(strcat('Neural Activation Distance Vs Threshold  :  Ch',num2str(curCh),'Curr',num2str(curCurr)))
                        ylabel('Activation Distance (\mum)')
                        xticklabels(num2str(uT'))
                        xlabel('Detection Thresholds (\muA)')
                        close(curFig)
                    end
                end 
            end
        end

        allDensityThresh = allDensityThresh./allDensityThreshN;



        for cInd = 1:numel(uCurr)
            curFig = figure();
            hold on
            for di = 1:10
                curDat = allDensityThresh(:,cInd,:,di);
                curMeans = squeeze(mean(curDat,1,'omitnan')); % average over channels
                plot(uThresholds,curMeans,'-o','LineWidth',2)
            end
            curCurr = uCurr(cInd);
            legend({'0-100','100-200','200-300','300-400','400-500','500-600','600-700','700-800','800-900','900-1000'})
            title(strcat('Overall Neural Activation Density Vs Threshold  :  Curr',num2str(curCurr)))
            ylabel('Neural Activation Density (neurons/\mum^3)')
            xlabel('Detection Thresholds (\muA)')
            close(curFig)
        end


        % Overall Activation Distances per current
        for cInd = 1:numel(uCurr)
            curCurr = uCurr(cInd); % process each current
            curFig = figure();
            hold on
            g = [];
            allDistances = [];
            for p = 1:size(uChanCurrs,1)
                if(uChanCurrs(p,2)==curCurr) % if valid current include in average
                    curDistances = populationDistances{p};
                    nDates = popDates{p};
                    curThresholds = daysThresholds(chInd,nDates);
                    
                    for d = 1:numel(curThresholds)
                        g = [g; curThresholds(d)*ones(numel(curDistances{d}),1)];
                        allDistances = [allDistances; curDistances{d}];
                    end
                end
            end
            
            boxplot(allDistances,g)
            title(strcat('Neural Activation Distance Vs Threshold  : Curr',num2str(curCurr)))
            ylabel('Activation Distance (\mum)')
            curThreshs = (unique(g));  
            xticklabels(num2str(curThreshs))
            xlabel('Detection Thresholds (\muA)')
            close(curFig)
        end

        
        % Overall radius of activation vs current
        curFig = figure();
        hold on
        for cInd = 1:numel(uCurr)
            curRadii = allRadiiThresh(:,cInd,:);
            curMean = squeeze(mean(curRadii,1,'omitnan'));
   
            p2 = plot(uThresholds,curMean, '-o','linewidth',2);
            p2.Color = currColors{cInd};
        end
        legend(num2str(uCurr))
        title(strcat('Overall Radius of Activation Vs Threshold'))
        ylabel('Activation Sphere Radius (\mum)')
        xlabel('Detection Thresholds (\muA)')
        close(curFig)
        
        
        
        
        
        
        
       %% Process neighbor contact site metrics
      
        
        % Update channel names/depths if not in the range of 1-9
        for k = 1:numel(uChan)
            if(uChan(k)==29)
                uChan(k) = 3;
            end
            if(uChan(k)==10)
                uChan(k)=4;
            end
            if(uChan(k)==27)
                uChan(k)=5;
            end
        end
        
        
        
        % Examine neural activation at threshold for neighboring contact
        % sites on given days
        for s = 1:numel(daysTrained) % iterate through each session
            curDat = NaN(numel(uChan),1);
            xlab = cell(numel(uChan),1);
            for c1 = 1:numel(uChan) % iterate through each channel
                curThresh1 = daysThresholds(c1,s); % Find day's threshold
                
                % Find matching current level & subsequent neural response near threshold 
                curInd1 = find(uCurr==curThresh1);
                if(~isempty(curInd1))
                    c1NeuCnt = groupedCurrChans(s,c1,curInd1);
                    curDat(c1) = c1NeuCnt;
                end
                xlab(c1) = strcat(num2str(uChan(c1)),{' '},'/',{' '},num2str(curThresh1));
            end
            
            figure()
            hold on
            bar(1:3,curDat)
            xticks(1:3)
            xticklabels(xlab')
            xlabel('Channel / Threshold Current (\muA)')
            ylabel('Activated Neurons')
            title(strcat('Channel Activation Comparison At Threshold - IMCS',{' '},num2str(icms),{' '},'Day',{' '},num2str(daysTrained(s))))
        end        
        
        
        
        % Examine regions of co-activation at threshold for neighboring contact
        % sites on given days
        for s = 1:numel(daysTrained) % iterate through each session
            curDat = NaN(3,1);
            xlab = cell(3,1);
            for c1 = 1:numel(uChan) % iterate through each channel
                curThresh1 = daysThresholds(c1,s); % Find day's threshold
                
                % Find matching current level & subsequent neural response near threshold 
                curInd1 = find(uCurr==curThresh1);
                if(~isempty(curInd1))
                    c1Raster = groupedRasters{s,c1,curInd1};
                    
                    for c2 = c1+1:numel(uChan) % iterate through each channel
                        curThresh2 = daysThresholds(c2,s); % Find day's threshold

                        % Find matching current level & subsequent neural response near threshold 
                        curInd2 = find(uCurr==curThresh2);
                        if(~isempty(curInd2))
                            c2Raster = groupedRasters{s,c2,curInd2};
                            rasterOverlap = c1Raster;
                            rasterOverlap(c2Raster==0)=0;
                            n = c1 + c2 - 2; % Silly but should work for small case of 3 channels max
                            curDat(n) = sum(rasterOverlap);
                            xlab(n) = strcat(num2str(uChan(c1)),{' '},'vs.',{' '},num2str(uChan(c2)));
                        end
                    end
                end
            end
            
            figure()
            hold on
            bar(1:3,curDat)
            xticks(1:3)
            xticklabels(xlab')
            xlabel('Compared Channels')
            ylabel('Co-Actived Regions')
            title(strcat('Channel Co-Activation Comparison At Threshold - IMCS',{' '},num2str(icms),{' '},'Day',{' '},num2str(daysTrained(s))))
        end  
        
        
        
        
        % Calculate the number of neurons relative between neighboring
        % contact sites - Comparisons are only made on the same day between
        % matching currents
        neighborComps = NaN(numel(daysTrained),3,8,numel(uCurr));
        for s = 1:numel(daysTrained) % iterate through each session
            for cInd = 1:numel(uCurr)
                for c1 = 1:numel(uChan) % 
                    % Load channel specific data
                    c1NeuCnt = groupedCurrChans(s,c1,cInd);
                    c1Thresh = daysThresholds(c1,s);
                    c1Raster = groupedRasters{s,c1,cInd};
                    if(isnan(c1NeuCnt)) % no data, go to next iteration
                        continue
                    end

                    % compare to other channels
                    for c2 = c1+1:numel(uChan)
                        c2NeuCnt = groupedCurrChans(s,c2,cInd);
                        c2Thresh = daysThresholds(c2,s);
                        c2Raster = groupedRasters{s,c2,cInd};
                        if(isnan(c2NeuCnt)) % no data, go to next iteration
                            continue
                        end
                        
                        % perform comparison
                        curDist = abs(uChan(c1) - uChan(c2))*60;
                        rasterOverlap = c1Raster;
                        rasterOverlap(c2Raster==0)=0;
                        absDiff = abs(c1NeuCnt - c2NeuCnt);
                        totalActive = (c1NeuCnt+c2NeuCnt);
                        relDiff = abs((c1NeuCnt - c2NeuCnt)/(c1NeuCnt+c2NeuCnt));
                        threshDiff = abs(c1Thresh - c2Thresh);
                        
                        n = c1 + c2 - 2; % Silly but should work for small case of 3 channels max
                        neighborComps(s,n,1,cInd) = curDist; % distance
                        neighborComps(s,n,2,cInd) = threshDiff; % difference in thresholds
                        neighborComps(s,n,3,cInd) = sum(rasterOverlap)/numel(c1Raster); % percentage of shared regions
                        neighborComps(s,n,4,cInd) = relDiff; % relative difference
                        neighborComps(s,n,5,cInd) = c1Thresh;
                        neighborComps(s,n,6,cInd) = c2Thresh;
                        neighborComps(s,n,7,cInd) = c1NeuCnt;
                        neighborComps(s,n,8,cInd) = c2NeuCnt;
                        neighborComps(s,n,9,cInd) = absDiff;
                        neighborComps(s,n,10,cInd) = totalActive;
                    end
                end
            end
        end
            
        
        % Use box and wisker plots for distances as they'll all be 60 or
        % 180
        
        % Overall - grouped by Threshold difference & Contact site Distance
        g1 = reshape(neighborComps(:,:,1,:),1,[]);
        g1Valid = ~isnan(g1);
        g1Valid(g1==0)=0;
        g2 = reshape(neighborComps(:,:,2,:),1,[]);
        g2Valid = ~isnan(g2);
        g = [g1',g2'];
        gValid = g1Valid & g2Valid;
        g = g(gValid,:);
        
        sharedRegs = 100*reshape(neighborComps(:,:,3,:),1,[]);
        sharedRegs = sharedRegs(gValid);
        relDiff = 100*reshape(neighborComps(:,:,4,:),1,[]);
        relDiff = relDiff(gValid);
        
        absDiff = reshape(neighborComps(:,:,9,:),1,[]);
        absDiff = absDiff(gValid);
        totAct = reshape(neighborComps(:,:,10,:),1,[]);
        totAct = totAct(gValid);
        
        % 
        % figure()% Shared Region percentage
        % hold on
        % boxplot(sharedRegs,g)
        % xlabel('Intercontact Distance (\mum) / Detection Threshold Difference')
        % ylabel('Percentage of Shared Regions')
        % title(strcat('ICMS',{' '},num2str(icms),{' '},'Neighboring Contact Analysis - Co-Activated Regions'))
        % 
        % figure()% Relative Difference
        % hold on
        % boxplot(relDiff,g)
        % xlabel('Intercontact Distance (\mum) / Detection Threshold Difference')
        % ylabel('Relative Difference in Region Activation')
        % title(strcat('ICMS',{' '},num2str(icms),{' '},'Neighboring Contact Analysis - Relative Difference'))
        % 
        % 
        % figure()% threshold variation ( violin plot)
        % hold on
        % g2v = g2(g2Valid);
        % g1v = g1(g1Valid);
        % violinplot(g2v, g1v);
        % xlabel('Intercontact Distance (\mum)')
        % ylabel('Detection Threshold Difference')
        % title(strcat('ICMS',{' '},num2str(icms),{' '},'Neighboring Contact Analysis - Relative Thresholds'))
        % 
        % 
        % figure()% absolute Difference - really needs stimulation current to be considered
        % hold on
        % boxplot(absDiff,g)
        % xlabel('Intercontact Distance (\mum) / Detection Threshold Difference')
        % ylabel('Absolute Difference in Region Activation')
        % title(strcat('ICMS',{' '},num2str(icms),{' '},'Neighboring Contact Analysis - Absolute Difference'))
        % 
        % figure()% absolute Difference - really needs stimulation current to be considered
        % hold on
        % boxplot(totAct,g)
        % xlabel('Intercontact Distance (\mum) / Detection Threshold Difference')
        % ylabel('Summed Total Region Activation')
        % title(strcat('ICMS',{' '},num2str(icms),{' '},'Neighboring Contact Analysis - Total Activation'))

%         % process for each contact site pair
%         for n= 1:3
%             switch n 
%                 case 1 
%                     ch1 = uChan(1);
%                     ch2 = uChan(2);
%                 case 2
%                     ch1 = uChan(1);
%                     ch2 = uChan(3);
%                 otherwise
%                     ch1 = uChan(2);
%                     ch2 = uChan(3);
%             end
%             chDist = 60*abs(ch1-ch2);
%             
%             
%             g = reshape(neighborComps(:,n,2),1,[]);
%             gValid = ~isnan(g);
%             g = g(gValid);
% 
%             sharedRegs = 100*reshape(neighborComps(:,n,3),1,[]);
%             sharedRegs = sharedRegs(gValid);
%             relDiff = 100*reshape(neighborComps(:,n,4),1,[]);
%             relDiff = relDiff(gValid);
%             absDiff = reshape(neighborComps(:,n,9,:),1,[]);
%             absDiff = absDiff(gValid);
%             totAct = reshape(neighborComps(:,n,10,:),1,[]);
%             totAct = totAct(gValid);
%         
%             figure()% Shared Region percentage
%             hold on
%             boxplot(sharedRegs,g)
%             xlabel('Intercontact Distance (\mum) / Detection Threshold Difference')
%             ylabel('Percentage of Shared Regions')
%             title(strcat('ICMS',{' '},num2str(icms),{' '},'Ch',{' '},num2str(ch1),{' '},'vs. Ch',{' '},num2str(ch2),{' '},'- Distance',{' '}, num2str(chDist) ,'\mum - Shared Regions'))
% 
%             figure()% Relative Difference
%             hold on
%             boxplot(relDiff,g)
%             xlabel('Intercontact Distance (\mum) / Detection Threshold Difference')
%             ylabel('Relative Difference in Regions')
%             title(strcat('ICMS',{' '},num2str(icms),{' '},'Ch',{' '},num2str(ch1),{' '},'vs. Ch',{' '},num2str(ch2),{' '},'- Distance',{' '}, num2str(chDist) ,'\mum - Relative Difference'))
%         
%         
%             figure()% absolute Difference - really needs stimulation current to be considered
%             hold on
%             boxplot(absDiff,g)
%             xlabel('Intercontact Distance (\mum) / Detection Threshold Difference')
%             ylabel('Absolute Difference in Region Activation')
%             title(strcat('ICMS',{' '},num2str(icms),{' '},'Ch',{' '},num2str(ch1),{' '},'vs. Ch',{' '},num2str(ch2),{' '},'- Distance',{' '}, num2str(chDist) ,'\mum - Absolute Difference'))
% 
%             figure()% absolute Difference - really needs stimulation current to be considered
%             hold on
%             boxplot(totAct,g)
%             xlabel('Intercontact Distance (\mum) / Detection Threshold Difference')
%             ylabel('Summed Total Region Activation')
%             title(strcat('ICMS',{' '},num2str(icms),{' '},'Ch',{' '},num2str(ch1),{' '},'vs. Ch',{' '},num2str(ch2),{' '},'- Distance',{' '}, num2str(chDist) ,'\mum - Total Activation'))
%         end
            














    
    %% Align ROIs between sessions
    
% % %     % calcuate roi centroid deviation, then attempt realignment of ROIs
% % %     HSS = 'HIT';
% % %     sharedROIs = cell(numel(trueSourceFolders),1);
% % %     for cd1 = 1:numel(trueSourceFolders)
% % %         
% % %         ind = find(chanCurrPermsOrig(:,3)==cd1,1);
% % %         activeRow = chanCurrPermsOrig(ind,:);
% % %         curpath = trueSourceFolders{cd1};
% % %         fileNameROIH = strcat(curpath,'\Ch',num2str(activeRow(1)),'-CUR',num2str(activeRow(2)),HSS,'\ROI.mat'); 
% % %         if(~isfile(fileNameROIH)) % check if older filename style is used or new format
% % %             fileNameROIH = strcat(curpath,'\Ch',num2str(activeRow(1)),'-CUR',num2str(activeRow(2)),'FREQ',num2str(activeRow(4)),'PULSE',num2str(activeRow(5)),HSS,'\ROI.mat');
% % %         end
% % %         
% % %         % Load ROIs
% % %         roiData = load(fileNameROIH);
% % %         regions = roiData.regions;
% % %         midActive = roiData.midActive;
% % %         
% % %         regionsVer2 = regions;
% % %         cents = cell(numel(regions),1);
% % %         for i = 1:numel(regions)
% % %             % load each region
% % %             curPixels = regions(i).PixelList;
% % %             curVol = false(size(midActive));
% % %             shiftedRoiVol = false(size(midActive));
% % %             for m = 1:size(curPixels,1)
% % %                 curVol(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
% % %             end
% % %             curVol = curVol(:,:,1:sharedDepth); % Trim volume
% % %             
% % %             zRange = curPixels(:,3);
% % %             zRange(zRange>sharedDepth)=[];
% % %             
% % %             % Shift Volume
% % %             for d = min(zRange):(max(zRange)-vOFF(cd))
% % %                 shiftedRoiVol(:,:,d) = spatial_interp_patchwarp(curVol(:,:,d), ptForms{cd}, 'euclidean', 1:512, 1:512);            
% % %             end
% % %             shiftedRoiVol = shiftedRoiVol(:,:,1:sharedDepth); % Trim volume
% % %             
% % %             shiftedRoiVol(overlay3D==0)=0;
% % %             curR = regionprops(shiftedRoiVol,'PixelList','Centroid');
% % %             if(numel(curR)>0)
% % %                 regionsVer2(i).PixelList = curR(1).PixelList;
% % %                 cents{i} =  curR(1).Centroid;
% % %             else
% % %                 regionsVer2(i).PixelList = NaN;
% % %                 cents{i} = NaN;
% % %             end
% % %         end
% % %         [regionsVer2.Centroid] = cents{:};
% % %         
% % %         sharedROIs{cd1} = regionsVer2;
% % %     end
% % %     
% % %     % Add all centroids to arrays
% % %     allCents = cell(numel(trueSourceFolders),1);
% % %     for i = 1:numel(trueSourceFolders)
% % %         curRegions = sharedROIs{i};
% % %         curCents = zeros(numel(curRegions),3);
% % %         for j = 1:numel(curRegions)
% % %             curCents(j,:) = curRegions(j).Centroid;
% % %         end
% % %         allCents{i} = curCents;
% % %     end
% % %     
% % %     
% % %     
% % %     
% % %     
% % % 
% % %     % Compare ROIs between all sessions to find those that overlap
% % %     sharedROIMapping = cell(numel(trueSourceFolders),1);
% % %     for cd = 1:numel(trueSourceFolders)
% % %         % compare ROI to other sessions and find nearest neighbors.  If
% % %         % there is overlap then valid shared region
% % %         curCents = allCents{cd};
% % %         curRegions = sharedROIs{cd};
% % % 
% % %         % Array of overlapping ROIs for current session
% % %         mappedRegions = zeros(numel(curRegions),numel(trueSourceFolders));
% % %         
% % %         for cd2 = 1:numel(trueSourceFolders)
% % %             if(cd ~= cd2)
% % %         
% % %                 targetCents = allCents{cd2};
% % %                 targetRegions = sharedROIs{cd2};
% % %                 targetCents(:,3) = targetCents(:,3)*1000;
% % %                 curCents(:,3) = curCents(:,3)*1000;
% % %                 idx = knnsearch(targetCents,curCents);
% % %                 
% % %                 for i = 1:numel(curRegions)
% % %                     if(~isnan(curRegions(i).PixelList))
% % %                         % Check overlap
% % %                         vol1 = false(512,512,sharedDepth);
% % %                         curPixels = curRegions(i).PixelList;
% % %                         for m = 1:size(curPixels,1)
% % %                             vol1(curPixels(m,2),curPixels(m,1))=1;
% % %                         end
% % % 
% % %                         vol2 = false(512,512,sharedDepth);
% % %                         curPixels = targetRegions(idx(i)).PixelList;
% % %                         for m = 1:size(curPixels,1)
% % %                             vol2(curPixels(m,2),curPixels(m,1))=1;
% % %                         end
% % % 
% % %                         overlap = vol1;
% % %                         overlap(vol2==0)=0;
% % % 
% % %                         if(sum(overlap,'all')>0)
% % %                             mappedRegions(i,cd2) = idx(i);
% % %                         end
% % %                     end
% % %                 end
% % %             end
% % %         end
% % %         sharedROIMapping{cd} = mappedRegions;
% % %     end
% % %      
% % %     
    
    
   
    
    
    
    % Initially ony focus on HIT Scenarios ... expand to MISS later
    % HSS = 'MISS';
    sharedROIs = cell(numel(trueSourceFolders),1);
    sharedROIMapping = cell(numel(trueSourceFolders),1);
    pixelInds = cell(numel(trueSourceFolders),1);
    for cd = 1:numel(trueSourceFolders) 

        if(mod(cd,2)==1)
            HSS = 'HIT';
        else
            HSS = 'MISS';
        end
        % Select file for dataset in current session - all share the same
        % core regions/raster base
        ind = find(chanCurrPermsOrig(:,3)==cd,1);
        activeRow = chanCurrPermsOrig(ind,:);
        curpath = trueSourceFolders{cd};
        fileNameROIH = strcat(curpath,'\Ch',num2str(activeRow(1)),'-CUR',num2str(activeRow(2)),HSS,'\ROI.mat'); 
        if(~isfile(fileNameROIH)) % check if older filename style is used or new format
            fileNameROIH = strcat(curpath,'\Ch',num2str(activeRow(1)),'-CUR',num2str(activeRow(2)),'FREQ',num2str(activeRow(4)),'PULSE',num2str(activeRow(5)),HSS,'\ROI.mat');
        end
        
        % Load ROIs
        roiData = load(fileNameROIH);
        regions = roiData.regions;
        regionsVer2 = regions;
        midActive = roiData.midActive;
        curPxInds = zeros(numel(regions),2); % construct it and expand to it....
        
        for i = 1:numel(regions)
            % load each region
            curPixels = regions(i).PixelList;
            curVol = false(size(midActive));
            shiftedRoiVol = false(size(midActive));
            for m = 1:size(curPixels,1)
                curVol(curPixels(m,2),curPixels(m,1),curPixels(m,3))=1;
            end
            curVol = curVol(:,:,1:sharedDepth); % Trim volume
            
            zRange = curPixels(:,3);
            zRange(zRange>sharedDepth)=[];
            
%             Shift Volume
            for d = min(zRange):(max(zRange)-vOFF(cd))
                shiftedRoiVol(:,:,d) = spatial_interp_patchwarp(curVol(:,:,d), ptForms{cd}, 'euclidean', 1:512, 1:512);            
            end
            shiftedRoiVol = shiftedRoiVol(:,:,1:sharedDepth); % Trim volume
            
            shiftedRoiVol(overlay3D==0)=0;
            curR = regionprops(shiftedRoiVol,'PixelList');
            if(numel(curR)>0)
                regionsVer2(i).PixelList = curR(1).PixelList;
            else
                regionsVer2(i).PixelList = NaN;
            end
            flatROI = reshape(shiftedRoiVol,1,[]);
            pxInds = find(flatROI);
            if(numel(pxInds)>1)
                curPxInds(i,1:numel(pxInds)) = pxInds;
            end
        end
        sharedROIs{cd} = regionsVer2;
        pixelInds{cd} = curPxInds;
    end
    


 %% Compare ROIs between all sessions to find those that overlap
 
    for cd = 1:numel(trueSourceFolders) 
        curRegions = sharedROIs{cd};
        curPixelInds = pixelInds{cd};
        
        % Array of overlapping ROIs
        mappedRegions = zeros(numel(curRegions),numel(trueSourceFolders));
        mappedRegionsAlternatives = cell(numel(curRegions),numel(trueSourceFolders));

        for i = 1:numel(curRegions)
            if(~isnan(curRegions(i).PixelList))
                
                % Load current pixel listing
                activeInds = curPixelInds(i,:);
                inds = find(activeInds);
                

                for cd2 = cd+1:numel(trueSourceFolders) 
                    if(cd2 ~= cd) % check between different sessions
                        searchPixelInds = pixelInds{cd2};
                        
                        matchingInd = zeros(size(searchPixelInds));
                        % perform search against all other regions
                        for h = 1:numel(inds)
                            matchingInd(searchPixelInds == activeInds(h))=1;
                        end
                        searchArr = sum(matchingInd,2);
                        if(sum(searchArr)>0)
                            % Region with maximal overlap is aligned
                            [~,mr]= max(searchArr);
                            mappedRegions(i,cd2) = mr(1); % catch for multiple maxes
                            mappedRegionsAlternatives{i,cd2} = find(searchArr); % all possible overlaps
                        end
                    end
                end
            end
        end
        
        
        % For each compared session, ensure that all mappings are unique
        % with no repeates.  If there are repeates find alternatives
        for t = 1:numel(trueSourceFolders) 
            fInds = find(mappedRegions(:,t));
            allN = numel(fInds);
            unN = numel(unique(mappedRegions(fInds,t)));
            
            if(allN>unN)
                % There are repeate instance in the directories find unique
                % alternatives
                
                % Find repeated entries
                uReg = unique(mappedRegions(fInds,t));
                for k = 1:numel(uReg)
                    invalid = uReg;
                    invalid(k)=[];
                    curReps = find(mappedRegions(:,t)==uReg(k));
                    if(numel(curReps)>1) % multiple instances- hence repeat
                        % find if any of the indicies had alternatives that
                        % are not currently in use
                        
                        possInds = zeros(numel(curReps),1);
                        possCnts = zeros(numel(curReps),1);
                        for cr = 1:numel(curReps)
                            curOps = mappedRegionsAlternatives{curReps(cr),t};
                            toRemove = ismember(curOps,invalid); % Remove indecies that are already in use elsewhere
                            curOps(toRemove==1)=[];
                            possInds(cr,1:numel(curOps)) = curOps;
                            possCnts(cr) = numel(curOps);
                        end
                        [~,ord2Fill] = sort(possCnts); % determine optimal order to fill regions
                        
                        for o = 1:numel(ord2Fill)
                            % Save first necissary region
                            CI = ord2Fill(o);
                            newInd = find(possInds(CI,:));
                            if(~isempty(newInd))
                                % valid Index available for mapping
                                mappedRegions(curReps(CI),t) = possInds(CI,newInd(1));

                                % update following possible regions to prevent
                                % mapping to same location
                                possInds(possInds == possInds(CI,1))= 0;
                            else
                                % no valid Index, clear located region
                                mappedRegions(curReps(CI),t) = 0;
                            end
                        end
                    end
                end
            end
        end
        
        
        sharedROIMapping{cd} = mappedRegions;
    end
    
    
    for cd1 = 1:numel(trueSourceFolders) 
        mappedRegionsDest = sharedROIMapping{cd1};
        for cd2 = 1:cd1-1
            mappedRegionsSource = sharedROIMapping{cd2} ;
            
            %construct reverse mapping of shared ROIs
            curMapping = mappedRegionsSource(:,cd1);
            mapInds = find(curMapping);
            
            for m = 1: numel(mapInds)
                mappedRegionsDest(curMapping(mapInds(m)),cd2) = mapInds(m);
            end
        end
        sharedROIMapping{cd1} = mappedRegionsDest;
    end
    
    

    %% Analyze the number of active neurons over time
    figure()
    hold on
    legText = cell(size(uChanCurrs,1),1);
    lgN = 0;
    uCurr = unique(uChanCurrs(:,2));
    uCurr(uCurr<3)=[];
    uCurr(uCurr>6)=[];
    uChan = unique(uChanCurrs(:,1));
    groupedCurrChansConsistency = NaN(numel(daysTrained),numel(uChan),numel(uCurr));

    currColors = {'r','m','b','y','g','k','#D95319','#7E2F8E','#77AC30','#4DBEEE','#EDB120','#7E2F8E','#FF00FF','#0072BD'}; % Hopefully these sets of colors and symbols will be enough
    chanSyms = {'o','+','square','x','diamond','*','^','<','>','|'};
    for cInd = 1:numel(uCurr)
        curCurr = uCurr(cInd); % process each current
        groupedCurrs = NaN(numel(daysTrained),size(uChanCurrs,1));

        for p = 1:size(uChanCurrs,1)
            if(uChanCurrs(p,2)==curCurr) % if valid current include in average

                chanInd = find(uChan==uChanCurrs(p,1));
                dateInds = correlationDates{p};
                if(~isempty(dateInds))
                    if(numel(dateInds)>2)
                        % Load stimulation scenario and associated Rasters
                        rasters = allRasters{p};

                        % load each associated shared ROI mapping
                        conAveraging = zeros(numel(rasters),numel(rasters));
                        for r = 1:numel(rasters)
                            curRast = rasters{r};
                            sharedMapping =  sharedROIMapping{dateInds(r)};
                            sharedInds = sharedMapping(curRast==1,:);
                            sharedMask = false(size(sharedInds));
                            sharedMask(sharedInds>0)=1;
                            timeSums = sum(sharedMask,1);
                            timeSums = timeSums(dateInds);
                            conVals = timeSums./size(sharedMask,1);
                            conVals(r) = NaN;
                            conAveraging(r,:) = conVals;
                        end
                        concistencyCnts = mean(conAveraging,1,'omitnan');
                        curDays = daysTrained(dateInds);

                        p1 = plot(curDays,concistencyCnts,'--'); % plot the number of neurons active for this channel and current
                        p1.Marker = chanSyms{chanInd};
                        p1.Color = currColors{cInd};
                        groupedCurrs(dateInds,p)=concistencyCnts;
                        lgN = lgN + 1;
                        legText{lgN} = strcat('CH',num2str(uChanCurrs(p,1)),'--CURR',num2str(uChanCurrs(p,2)));
                        groupedCurrChansConsistency(dateInds,chanInd,cInd) = concistencyCnts;
                    end
                end
            end
        end    
        % plot grouped neural activation for this current
        grpTrace = mean(groupedCurrs,2,'omitnan');
        curTrace2 = grpTrace;
        curTrace2(isnan(curTrace2))=0;
        ctX = find(curTrace2);
        ctY = curTrace2(ctX);
        if(~isempty(ctX))
            p2 = plot(daysTrained(ctX),ctY,'linewidth',2);
            p2.Color = currColors{cInd};
            lgN = lgN + 1;
            legText{lgN} = strcat('Grouped CURR',num2str(curCurr));
        end
    end
    toKeep = false(numel(legText),1);
    for k = 1:numel(legText)
        if(~isempty(legText{k}))
            toKeep(k)=1;
        end
    end
    legText = legText(toKeep);
    legend(legText)
    title('Neural Activation Concistency Vs Time')
    ylabel('Concistency of Activation')
    xlabel('Days of Training')




        
%         
%         
%         % Overall neural activation with respect to threshold
%         figure()
%         hold on
%         for cInd = 1:numel(uCurr)
%             curCurr = uCurr(cInd); % process each current
%             groupedCurrsThresh = NaN(numel(uThresh),size(uChanCurrs,1));
%             for p = 1:size(uChanCurrs,1)
%                 if(uChanCurrs(p,2)==curCurr) % if valid current include in average
% 
%                     chanInd = find(uChan==uChanCurrs(p,1));
%                     dateInds = correlationDates{p};
%                     curThresholds = daysThresholds(chanInd,dateInds);
%                     threshInds = zeros(numel(curThresholds),1);
%                     
%                     
%                     
%                     for c = 1:numel(curThresholds)
%                         threshInds(c) = find(uThresh == curThresholds(c));
%                     end
%                     utInds = unique(threshInds);
%                     neuCnt = correlationMatriciesNeuronCount{p};
%                     neuCnts2 = zeros(numel(utInds),1);
%                     neuCnts2Num = zeros(numel(utInds),1);
%                     indsInds = zeros(numel(threshInds),1); % I've gone too deep with indexing, send help
%                     for c = 1:numel(curThresholds)
%                         indsInds(c) = find(utInds == threshInds(c));
%                     end
% 
%                     for c = 1:numel(curThresholds)
%                         if(neuCnts2(indsInds(c))==0)
%                             neuCnts2(indsInds(c)) = neuCnt(c);
%                             neuCnts2Num(indsInds(c)) = 1;
%                         else
%                             neuCnts2(indsInds(c)) = neuCnts2(indsInds(c)) + neuCnt(c);
%                             neuCnts2Num(indsInds(c)) = neuCnts2Num(indsInds(c))+1;
%                         end
%                     end
%                     neuCnts2 = neuCnts2./neuCnts2Num; % average for multiple sessions at same threshold
%                     
%                     
%                     
%                     % Load stimulation scenario and associated Rasters
%                     rasters = allRasters{p};
% 
%                     % load each associated shared ROI mapping
%                     conAveraging = zeros(numel(rasters),numel(rasters));
%                     for r = 1:numel(rasters)
%                         curRast = rasters{r};
%                         sharedMapping =  dDiv{dateInds(r)};
%                         sharedInds = sharedMapping(curRast==1,:);
%                         sharedMask = false(size(sharedInds));
%                         sharedMask(sharedInds>0)=1;
%                         timeSums = sum(sharedMask,1);
%                         timeSums = timeSums(dateInds);
%                         conVals = timeSums./size(sharedMask,1);
%                         conVals(r) = NaN;
%                         conAveraging(r,:) = conVals;
%                     end
%                     concistencyCnts = mean(conAveraging,1,'omitnan');
%                     curDays = daysTrained(dateInds);
% 
%                     
%                     
%                     
% 
%                     p1 = plot(uThresh(utInds),concistencyCnts,'--'); % plot the number of neurons active for this channel and current
%                     p1.Marker = chanSyms{chanInd};
%                     p1.Color = currColors{cInd};
%                     groupedCurrsThresh(utInds,p)=concistencyCnts;
%                     lgN = lgN + 1;
%                     legText{lgN} = strcat('CH',num2str(uChanCurrs(p,1)),'--CURR',num2str(uChanCurrs(p,2)));
%                 end
%             end    
%             % plot grouped neural activation for this current
%             grpTrace = mean(groupedCurrsThresh,2,'omitnan');
%             curTrace2 = grpTrace;
%             curTrace2(isnan(curTrace2))=0;
%             ctX = find(curTrace2);
%             ctY = curTrace2(ctX);
%             p2 = plot(uThresh(ctX),ctY,'linewidth',2);
%             p2.Color = currColors{cInd};
%             lgN = lgN + 1;
%             legText{lgN} = strcat('Grouped CURR',num2str(curCurr));
%         end
%         legend(legText)
%         title(strcat('Consistency Vs Threshold', {'  '},HSS))
%         ylabel('Consistency')
%         xlabel('Threshold (\muA)')
% 
% 



        
        % Save data for analysis
        save(strcat('AlternativeConcistencyDataALL_ver7',HSS,num2str(icms),'.mat'),'uChanCurrs','weekLim','allPopulation',...
            'allDensity','uCurr','daysTrained','enclosingRadius','popDates','populationDistance','populationHistogram',...
            'populationDistances','groupedCurrChans','groupedCurrChansThresh','allRadiiThresh','allDensityThresh',...
            'daysThresholds','groupedCurrChansThreshBaseMean','groupedCurrChansThreshBaseMin',...
            'groupedCurrChansThreshBaseMax','groupedCurrChansThreshBaseSTD','groupedCurrChansThreshBaseMeanUpper10',...
            'groupedCurrChansBaseMean','groupedCurrChansBaseMin','groupedCurrChansBaseMax','groupedCurrChansBaseSTD',...
            'groupedCurrChansBaseMeanUpper10','groupedCurrChansNorm','neighborComps','correlationDates',...
            'correlationMatricies','correlationMatricies3D','correlationMatriciesBase','correlationMatricies3DBase',...
            'diffMatricies','diffMatriciesSTD','diffMatriciesVals','diffMatricies3D','correlationMatriciesROI',...
            'correlationMatriciesROI3D','correlationMatriciesNeuronCount','correlationMatriciesNeuronCountNorm',...
            'sharedROIMapping','groupedCurrChansConsistency','allRasters','regions');
%     end
end  
