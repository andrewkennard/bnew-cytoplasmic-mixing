function [avgbr, avgbrseg] = getImgBrightness(img,segcont,varargin)
% get average brightness of whole image and segmented contour

% use mean or median brightness for average
avgtype = 'mean';
dodisplay = 0;

for vc = 1:2:length(varargin)
    switch(varargin{vc})
        case('avgtype')
            avgtype = varargin{vc+1};
        case('dodisplay')
            dodisplay = varargin{vc+1};
    end
end

mask = poly2mask(segcont(:,1),segcont(:,2),size(img,1),size(img,2));
if (dodisplay>1)
    imshow(mask,[])
    drawnow
end

if (strcmp(avgtype,'mean'))
    meanint = mean(img(mask(:)));
    avgbrseg = meanint;
    avgbr = mean(img(:));
elseif (strcmp(avgtype,'median'))
    meanint = median(double(img(mask(:))));
    avgbrseg = meanint;
    avgbr = median(double(img(:)));
else
    error('bad value for avgtype: ', avgtype)
end

end