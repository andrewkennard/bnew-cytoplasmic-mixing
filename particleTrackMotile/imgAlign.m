function [offset,ccval,max_c,mean_c] = imgAlign(img1,img2,varargin)
% align a template image img1 to a larger image img2
% based on a rectangle (rect) subset of img1
% using cross-correlation
% default rectangles span entire image
% optionally, supply a mask to apply to cropped img1 region first
% mask must match rect1 in size
% WARNING: does not work quite right if rectangles are non-integer

rect1 = [1,1,size(img1,2),size(img1,1)]; % rectangle to crop img1
rect2 = [1,1,size(img2,2),size(img2,1)]; % rectangle to crop img2
usemask = 0; % apply mask before aligning
dodisplay = 1; % display aligned images
subpixel = 0; % do subpixel alignment
subpixscl = 0.01; % resolution for subpixel alignment

for vc = 1:2:length(varargin)
    switch (varargin{vc})
        case('mask')
            mask = varargin{vc+1};
            usemask = 1;
        case('rect1')
            rect1 = varargin{vc+1};
        case('rect2')
            rect2 = varargin{vc+1};
        case('dodisplay')
            dodisplay = varargin{vc+1}; 
        case('subpixel')
            subpixel = varargin{vc+1};
        case('subpixscl')
            subpixscl = varargin{vc+1};    
    end
end

if (usemask & (size(mask,1)~=rect1(3)+1 | size(mask,2)~=rect1(4)+1))
    error('mask size must match rect1 size')
end

% subtract out background
m1 = mode(double(img1(:)));
m2 = mode(double(img2(:)));
img1 = img1-m1;
img2 = img2-m2;

img1C = imcrop(img1,rect1);
img2C = imcrop(img2,rect2);

if (usemask)
    img1C = double(img1C).*mask;
end

%figure;
%imshow(img1C,[],'initialmagnification','fit')
%figure;
%imshow(img2C,[],'initialmagnification','fit')

%% do the cross-correlation

    ccval = normxcorr2(img1C,img2C);    

%figure, surf(ccval), shading flat

% if (normarea)
%     %% normalize by overlap area
%     w = size(img1C,2); h = size(img1C,1);
%     if (isequal(size(img1C),size(img2C)))
%         xshift = -(w-1):(w-1);
%         yshift = -(h-1):(h-1);
%         [XS,YS] = meshgrid(xshift,yshift);
%         overlaparea = (w-abs(XS)).*(h-abs(YS));
%         ccval = ccval./overlaparea;
%         ccval = ccval./ccval(h,w);
%     else
%         error('Alignment normalized by overlap area not implemented for unequal img sizes')
%     end
% end
% %% check the peak
% goodind = find(overlaparea>prod(size(img1C))/2);
% [max_c, imax] = max(abs(ccval(goodind)));
% [ypeak, xpeak] = ind2sub(size(ccval),goodind(imax(1)))
%% check the peak
[max_c, imax] = max(abs(ccval(:)));
[ypeak, xpeak] = ind2sub(size(ccval),imax(1));

mean_c = mean(abs(ccval(:)));
%max_c/mean_c


% get the overall offset
corr_offset = [(xpeak-size(img1C,2))
               (ypeak-size(img1C,1))];
           
if (subpixel)
    % do some hackish subpixel alignment
    xspan = 10;
    yspan = 10;
    rect = [xpeak-xspan,ypeak-yspan,2*xspan,2*yspan];
    rect = expandRect(rect,1,[1,size(ccval,2)],[1,size(ccval,1)]);
    
    I = imcrop(ccval,rect);
    xvals=1:rect(3)+1;
    yvals = 1:rect(4)+1;
    [X,Y] = meshgrid(xvals,yvals);
    
    xvalsI=1:subpixscl:rect(3)+1;
    yvalsI = 1:subpixscl:rect(4)+1;
    [XI,YI] = meshgrid(xvalsI,yvalsI);
    
    II = interp2(X,Y,I,XI,YI,'spline');
    
    [max_cI, imaxI] = max(abs(II(:)));
    [ypeakI, xpeakI] = ind2sub(size(II),imaxI(1));
    
    % adjust
    xpeakA = rect(1)+xvalsI(xpeakI)-1;
    ypeakA = rect(2) +yvalsI(ypeakI)-1;  
    
    corr_offset = corr_offset + [xpeakA ypeakA]' - [xpeak ypeak]';
end
           
           
rect_offset = [rect2(1)-rect1(1) 
                rect2(2)-rect1(2)];

offset = corr_offset + rect_offset;
xoffset = offset(1);
yoffset = offset(2);

%% display images one on top of another
if (dodisplay)
xbegin = round(xoffset+1);
xend   = round(xoffset+ size(img1,2));
ybegin = round(yoffset+1);
yend   = round(yoffset+size(img1,1));

xmin = min(xbegin,1);
ymin = min(ybegin,1);
xmax = max(xend,size(img2,2));
ymax = max(yend,size(img2,1));

% pad image2 to enable overlap with image 1
img2pad = uint8(zeros(ymax-ymin+1,xmax-xmin+1));

% where in padded image does img2 start
x0 = -xmin+2;
y0 = -ymin+2;

img2pad(y0:y0+size(img2,1)-1,x0:x0+size(img2,2)-1) = img2;

%imshow(img2pad,[],'initialmagnification','fit')
%%
recovered = uint8(zeros(size(img2pad)));

ybegin0 = ybegin+y0-1;
yend0 = yend+y0-1;
xbegin0 = xbegin +x0-1;
xend0 = xend+x0-1;

recovered(ybegin0:yend0,xbegin0:xend0,:) = img1;
%figure, imshow(recovered, 'initialmagnification','fit')

%%
%imshowpair(img2pad,recovered)
img12 = imfuse(img2pad,recovered);
imshow(img12,[],'initialmagnification','fit')
drawnow;
end

offset = [xoffset, yoffset];
