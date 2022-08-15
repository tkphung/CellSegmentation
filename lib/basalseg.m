function [BCBMASK,varargout] = basalseg(BCBIMG,pr,varargin)
%basalseg(BCBIMG,varargin): Segmentation of basal cell bodies developed
%specifically for 20x OR 63x magnification Maximum Intensity Projection of
%Phalloidin-stained Z-stack. Borders for cells are dilated by a disk of
%0.75 um in radius.
% 
% [BCBMASK,varargout] = basalseg(BCBIMG,varargin)
%   INPUTS:
%       BCBIMG- Basal Cell Image file (phalloidin stained, 20x mag)
%       pr-     pixel resolution (microns/pixel)
%       varargin: 'visualize' followed by true or false (default)
%
%   OUTPUTS:
%       BCBMASK- binary mask of cell bodies (body = 1)
%       varargout: {1} all borders
%
% Created by Thien-Khoi N. Phung (April 17, 2020)

% Deal with VARARGIN
visual_flag = false;
if ~isempty(varargin)
    for jz = 1:2:numel(varargin)
        switch varargin{jz}
            case 'visualize'
                visual_flag = varargin{jz+1};
            otherwise
                error('BCB ERROR: Check your varargins.')
        end
    end
end % varargin


% Remove image noise
    wfilt = 2;%round(0.5/pr);
    CBnorm = wiener2(BCBIMG,[wfilt wfilt]);

% Find seeds
    CBsmooth = imgaussfilt(CBnorm,2.32/pr); % sigma determined from visual optimization
    seeds    = imregionalmin(CBsmooth); % use min if cell cyto is black
                                        % use max if cell cyto is white

% Watershed Segmentation
    % Smooth image and impose seeds as minimum value in image
    CBsmooth = imgaussfilt(CBnorm,0.37/pr); % sigma determined from visual optimization
    CBimpmin = imimposemin(CBsmooth,seeds);
    % Watershed algorithm
    CBwshed  = watershed(CBimpmin);

% Clear cells touching the border
    CBwshedbord  = imclearborder(CBwshed);

% Dilate cell boundaries (based on pixel resolution)
    % Cell outlines
    celloutl = CBwshedbord==0;
    celloutl = ~imfill(~celloutl,'holes');
    
    % Identify objects
    bcellid = bwlabel(~celloutl);
    
% Filter cells by size
    % Calculate areas of regions
    areas = cell2mat(struct2cell(regionprops(logical(bcellid),'Area')));
    % Remove based on mean
    msegareas = mean(areas);
    % Exclusion window
    exclcell = bsxfun(@or,areas<0.20*msegareas,...
                          areas>2.5*msegareas);
    % Exclude cells
    bcellid(ismember(bcellid,find(exclcell))) = uint16(0);

% OUTPUT: Basal Cell Body Mask
    BCBMASK      = bcellid>0; % Cell Body Mask
    varargout{1} = CBwshed;   % All borders logical mask

% (OPTION) Visualization of segmentation
if visual_flag
    I6 = labeloverlay(imadjust(BCBIMG),BCBMASK,'Transparency',0.50);
    figure('WindowStyle','docked','NumberTitle','off','name','BCB Body')
    axt = axes('Units', 'normalized', 'Position', [0 0 1 1]);
    imshow(I6,'Parent',axt)
    
    I7 = labeloverlay(imadjust(BCBIMG),~BCBMASK,'Transparency',0.50);
    figure('WindowStyle','docked','NumberTitle','off','name','BCB Border')
    axt = axes('Units', 'normalized', 'Position', [0 0 1 1]);
    imshow(I7,'Parent',axt)
end % visual_flag

end % basalseg()