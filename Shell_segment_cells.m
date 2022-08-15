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