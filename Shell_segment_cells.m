% Shell_segment_cells.m
% This script takes in a *.tif file and segments the cell boundaries.
%
% Thien-Khoi N. Phung (August 15, 2022)

addpath('lib')

%% Example of apical cell segmentation
% Image file path
apicalfile = 'test_images/U13_72_C2_4_ACB.tif';

% Pixel resolution of the image (hardcoded)
pr = 0.586; % um/pixel

% Segmentation
% Load ACB file
ACBimage = imread(apicalfile);

% Segment Apical Cell Boundaries (function assumes 20x mag)
apicalbody = apicalseg(ACBimage,pr,'visualize',true);


%% Example of basal cell segmentation
% Image file path
basalfile = 'test_images/U13_72_C2_4_BCB.tif';

% Pixel resolution of the image (hardcoded)
pr = 0.586; % um/pixel

% Segmentation
% Load ACB file
BCBimage = imread(basalfile);

% Segment Apical Cell Boundaries (function assumes 20x mag)
basalbody = basalseg(BCBimage,pr,'visualize',true);


%% Quanity & Visualize cell morphology
% Calculate apical cell metrics
    areas        = cell2mat(struct2cell(regionprops(apicalbody,'Area'))');
    mals         = cell2mat(struct2cell(regionprops(apicalbody,'MajorAxisLength'))');
    mils         = cell2mat(struct2cell(regionprops(apicalbody,'MinorAxisLength'))');
    aspectratio  = mals./mils;

% Map cell area
    ACBlabel      = bwlabel(apicalbody);
    [lia,locb]    = ismember(ACBlabel,1:numel(aspectratio));
    ACBlabel(lia) = areas(locb(lia));
    
    figure('WindowStyle','docked','NumberTitle','off','name',...
           'Apical Area')
    axt = axes('Units', 'normalized', 'Position', [0 0 1 1]);
    imshow(imadjust(ACBimage),[])
    hold on
    axtoo = axes('Units', 'normalized', 'Position', [0 0 1 1]);
    imagesc(ACBlabel,'AlphaData',imerode(apicalbody,1))
    colormap(axtoo,'spring')
    cb=colorbar;
    cb.Position = cb.Position + 1e-10;
    caxis([500 2000])
    axis equal tight off
    linkaxes([axt axtoo])

% Map cell aspect ratio
    ACBlabel      = bwlabel(apicalbody);
    [lia,locb]    = ismember(ACBlabel,1:numel(aspectratio));
    ACBlabel(lia) = aspectratio(locb(lia));
    
    figure('WindowStyle','docked','NumberTitle','off','name',...
           'Apical Aspect Ratio')
    axt = axes('Units', 'normalized', 'Position', [0 0 1 1]);
    imshow(imadjust(ACBimage),[])
    hold on
    axtoo = axes('Units', 'normalized', 'Position', [0 0 1 1]);
    imagesc(ACBlabel,'AlphaData',imerode(apicalbody,1))
    colormap(axtoo,'winter')
    cb=colorbar;
    cb.Position = cb.Position + 1e-10;
    caxis([1.5 3.0])
    axis equal tight off
    linkaxes([axt axtoo])
