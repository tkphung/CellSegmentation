function ACBMASK = apicalseg(ACBIMG,pr,varargin)
%apicalseg(ACBIMG,varargin): Segmentation of apical cell bodies developed
%specifically for 20x magnification Maximum Intensity Projection of
%Phalloidin-stained Z-stack
% 
% ACBMASK = apicalseg(ACBIMG,varargin)
%   INPUTS:
%       ACBIMG- Apical Cell Image file (phalloidin stained, 20x mag)
%       pr-     pixel resolution (microns/pixel)
%       varargin: 'visualize' followed by true or false (default)
%
%   OUTPUTS:
%       ACBMASK- binary mask of cell bodies (body = 1)
%
% References used:
% Marker-controlled watershed segmentation
% https://www.mathworks.com/help/images/marker-controlled-watershed-segmentation.html
% 
% Created by Thien-Khoi N. Phung (April 16, 2020)

% Deal with VARARGIN
visual_flag = false;
if ~isempty(varargin)
    for jz = 1:2:numel(varargin)
        switch varargin{jz}
            case 'visualize'
                visual_flag = varargin{jz+1};
            otherwise
                error('ACB ERROR: Check your varargins.')
        end
    end
end % varargin


% Opening-Closing by Reconstruction (removing small blemishes without
% affecting the overall shapes of the objects)
    se      = strel('disk',round(1.465/pr));
    Ie      = imerode(ACBIMG,se);
    Iobr    = imreconstruct(Ie,ACBIMG);
    Iobrd   = imdilate(Iobr,se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);

% Regional Minima to obtain foreground markers
    fgm = imregionalmin(Iobrcbr);
% Impose minima
    Iimpmin = imimposemin(ACBIMG,fgm);

% Watershed
    L = watershed(Iimpmin);
    borders = L==0;
% Connected Components to eliminate small rings (keep only the largest
% network- which should be the apical cell borders)
    CC = bwconncomp(borders);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = max(numPixels);
    borders = zeros(size(borders));
    borders(CC.PixelIdxList{idx}) = 1;

% Generate "quality" metrics for borders
% endpoint-to-endpoint length & pixel intensity
    % Count neighbors for each border (using 4 neighbor defn)
        four = zeros(3); four(2:2:end) = 1;
        cnbrs = conv2(borders,four,'same');
    % cnbrs>2 means a corner, Remove corners
        bordersegs = borders;
        bordersegs(cnbrs>2) = 0;
    % bwlabel the line segments
        labsegs = bwlabel(bordersegs,4);
    % Segment lengths calculated as number of pixels (area)
        seglengths = cell2mat(struct2cell(regionprops(labsegs,'area')));
    % All endpoints of line segments as linear indices 
        cnbrs       = conv2(bordersegs,four,'same');
        [erow,ecol] = find(cnbrs==1);
        eptslist    = sub2ind(size(cnbrs),erow,ecol);
    % line segment pixels
        pxlist = struct2cell(regionprops(labsegs,'PixelIdxList'));
    % Find end points in each line segment & calculate distance
        epdist = ones(size(pxlist));
        pxint  = zeros(size(pxlist));
        for jz = 1:numel(pxlist)
            
            % Calculate distance
            if numel(pxlist{jz})>1
                % Find endpoints in the segment
                endpts = ismember(pxlist{jz},eptslist);
                % Convert from linear indices to [row,col]
                [epx,epy] = ind2sub(size(cnbrs),pxlist{jz}(endpts));
                % Calculate distance
                epdist(jz) = sqrt(diff(epx)^2 + diff(epy)^2);
            end

            % mean intensity
            pxint(jz) = mean(ACBIMG(pxlist{jz}));
        end % working through line segments

% Calculate tortuosity of segments
    % (length of curve [number of pixels])/(endpoint to endpoint length)
        tort    = seglengths(:)./epdist(:);
    % Assign Tortuosity and Pixel Intensity to image map
        tortimage       = labsegs;
        pxintimage      = labsegs;
        seglengthsimage = labsegs;
        epdistimage     = labsegs;
        [rlog,ridx]     = ismember(labsegs,1:numel(tort));
        tortimage(rlog)       = tort(ridx(rlog));
        pxintimage(rlog)      = pxint(ridx(rlog));
        seglengthsimage(rlog) = seglengths(ridx(rlog));
        epdistimage(rlog)     = epdist(ridx(rlog));

% Define bad borders
    % Thresholds
    tortthresh   = 1.5;                         % Tortuosity
    pxthresh     = mean(pxint) - std(pxint);    % Pixel Intensity
    pxthreshlow  = mean(pxint) - 2*std(pxint);    
    seglengthlow = mean(seglengths);            % Segment Length
    
%     % Version 1
%     % (Tortuosity>torthresh AND MeanPixelIntensit<pxthresh)
%     % OR (0<MeanPixelIntensity<pxthreshlow)
%     % badborders  = (tortimage>tortthresh & pxintimage<pxthresh) | (pxintimage<pxthreshlow & pxintimage>0);

    % Version 2
    % Borders of interest for exclusion
    % (Tortuosity>torthresh AND MeanPixelIntensit<pxthresh AND SegLength>seglengthlow)
    % OR (0<MeanPixelIntensity<pxthreshlow)
    badborders  = (tortimage>tortthresh & pxintimage<pxthresh & seglengthsimage>seglengthlow); %| (pxintimage<pxthreshlow & pxintimage>0);

% Remove bad borders
    keepborders = logical(borders - badborders);
% Now remove any free floating branches
    kbord = bwmorph(keepborders,'spur',inf);
% And remove any small segmentation areas
    kbord = bwareaopen(kbord,round(2.1462/(pr^2)));

% Clear cells touching the border & label cells
    acells  = imclearborder(~kbord);
    acellid = bwlabel(acells);

% Filter cells by size
    % Calculate areas of regions
    areas = cell2mat(struct2cell(regionprops(logical(acellid),'Area')));
    % Remove based on mean
    msegareas = mean(areas);
    % Exclusion window
    exclcell = bsxfun(@or,areas<0.05*msegareas,...
                          areas>3.0*msegareas);
    % Exclude cells
    acellid(ismember(acellid,find(exclcell))) = uint16(0);
    
% OUTPUT: Apical cell bodies logical mask
    ACBMASK = acellid>0;
    
% (OPTION) Visualization of segmentation
if visual_flag
    I6 = labeloverlay(imadjust(ACBIMG),ACBMASK,'Transparency',0.50);
    figure('WindowStyle','docked','NumberTitle','off','name','ACB Body')
    axt = axes('Units', 'normalized', 'Position', [0 0 1 1]);
    imshow(I6,'Parent',axt)
    
    I7 = labeloverlay(imadjust(ACBIMG),~ACBMASK,'Transparency',0.50);
    figure('WindowStyle','docked','NumberTitle','off','name','ACB Border')
    axt = axes('Units', 'normalized', 'Position', [0 0 1 1]);
    imshow(I7,'Parent',axt)
end % visual_flag

end % apicalseg()
