function [jumps,cellmv,opt] = checkStageJumps(cellmv,options)
%% check for stage jumps by aligning images
% use largest connected component after thresholding for image alignment

% default options
opt = struct();
% which frames to check for jump (before that frame)
opt.checkforjump = ones(1,length(cellmv));
opt.checkforjump(1) = 0;
% display alignment
opt.dodisplay = 0;
% thresholding scale (multiply by the automatically set threshold)
opt.threshscl=0.4; 
% max out pixels above this percentile
opt.adjmax = 1;
opt.adjmin=0.001;
% minimal area of cell (in px). throw out anything smaller than this
opt.minarea = 1e3;
% dilation around thresholded area in pixels
opt.dilate = 50;
% maximum fractional jump in avg image brightness; if jump is too big,
% then assume phase is being turned on/of and DO NOT align images
% THIS WILL BREAK IF THERE IS A STAGE JUMP DURING PHASE ON/OFF
opt.brightjump = 0.5;
% maximum offset (in px) when alignining consecutive frames, to assume no stage
% jump
opt.maxalignoff = 2;

% use the jumps to set cell positions
opt.setCellPos = 1;

% more printing
opt.verbose=0;

% calculate segmentation via thresholding, as opposed to using existing
% segcont
opt.getsegcont = 1;

% input options
opt = copyStruct(options,opt);

%%
jumps = zeros(length(cellmv),2);

if (opt.verbose)
    disp(sprintf('Checking stage jump for %d frames',nnz(opt.checkforjump)))
end

checkframes = fliplr(find(opt.checkforjump));

if (~opt.getsegcont & isempty(cellmv(checkframes(1)).segcont))
    error('Cannot do alignment if first segmented contour is not defined')
end
    
currect = [];
%%
for fc = checkframes
    % align rectangle around segmented contour on frame fc
    % against full image of frame fc-1
    % if segmented contour is undefined (and no getsegcont), use the
    % rectangle based on the previous alignment, shifted according to jump
    
    if (fc==1); continue;end
    
    if (opt.verbose)
        disp(sprintf('Checking for frame jump before frame %d,',fc))
    end
    imgprev = imread(cellmv(fc-1).fname, cellmv(fc-1).fpage);
    img = imread(cellmv(fc).fname, cellmv(fc).fpage);
%%
    if (opt.getsegcont)
        %%
        m1 = mean(imgprev(:));
        m2 = mean(img(:));
        
        
        img = imadjust(img,stretchlim(img,[0.1,1]));
        imgprev = imadjust(imgprev,stretchlim(img,[0.1,1]));
        %%
        if(abs(m2/m1-1)>opt.brightjump)
            % brightness jumps too much to make alignment meaningful, assume no stage jump
            jumps(fc,:) = [0,0];
            disp(sprintf('Too big a brightness jump on frame %d. Cant align. Assume no stage jump. %f %',fc,m1,m2))
            continue;
        end
        
        %imshow(imgprev,[],'initialmagnification','fit')
        
        %% thresholding
        % normalize image
        I2 = double(img);
        I2 = I2/max(I2(:));
        I2 = imadjust(I2,stretchlim(img,[opt.adjmin,opt.adjmax]));
        %h = fspecial('gaussian', 4, 4)
        %I3 = imfilter(I2,h,'replicate');
        %imshow(I3,[]);
        %
        level = graythresh(I2);
        bw = im2bw(I2,level*opt.threshscl);
        bw = bwareaopen(bw, opt.minarea);
        if (opt.dodisplay>1)
            imshow(bw,[]);
        end
        %% dilate using distance to nearest pixel
        bw2 = bwdist(bw) <= opt.dilate;
        %imshow(bw2)
        
        % pick out connected components
        concomp = bwconncomp(bw2, 8);
        
        % largest connected component
        areas = cellfun(@(x) length(x), concomp.PixelIdxList);
        [maxa,maxind] = max(areas);
        bwcc = zeros(size(bw2));
        bwcc(concomp.PixelIdxList{maxind}) = 1;
        
        if (opt.dodisplay>=2)
            imshow(bwcc,'initialmagnification','fit')
        end
        
        %% draw rectangle around largest connected component
        xvals = 1:size(bwcc,2);
        yvals = 1:size(bwcc,1);
        [X,Y] = meshgrid(xvals,yvals);
        
        inds = find(bwcc(:));
        xmin = min(X(inds));
        xmax = max(X(inds));
        ymin = min(Y(inds));
        ymax = max(Y(inds));
        rect = [xmin ymin xmax-xmin ymax-ymin];
        currect = rect;
    else            
        if (~isempty(cellmv(fc).segcont))
            % align using the largest segmented area 
            for sc = 1:length(cellmv(fc).segcont)
                tmp = polygeom(cellmv(fc).segcont{sc}(:,1),cellmv(fc).segcont{sc}(:,2));
                areas(sc) = tmp(1);
            end
            [a,b] = max(areas);
            % get rectangle from saved segmentation contour
            xmin = min(cellmv(fc).segcont{b}(:,1));
            xmax = max(cellmv(fc).segcont{b}(:,1));
            ymin = min(cellmv(fc).segcont{b}(:,2));
            ymax = max(cellmv(fc).segcont{b}(:,2));
            rect = [xmin ymin xmax-xmin ymax-ymin];
            currect = rect;
        else
            rect = currect;
        end
    end
    
    if (~isequal(size(rect),[1,4]))
        error('bad rectangle')
    end
    %%
    if (opt.dodisplay>=2)
        imshow(img,[],'initialmagnification','fit')
        hold all
        rectangle('Position',rect,'EdgeColor','y')
        hold off
    end
    
    %%
    [offset,ccval,max_c,mean_c] = imgAlign(imgprev,img,'rect1',rect,'dodisplay',opt.dodisplay);  
    if (opt.dodisplay)
        hold all
        text(100,100,sprintf('Frame %d',fc),'Color','r','FontSize',20)
        hold off
        drawnow;
    end
        
    %% if offset is less than cutoff, assume no stage jump
    if norm(offset)<opt.maxalignoff
        jumps(fc,:) = [0,0];        
    else
        disp(sprintf('Stage jump detected before frame %d. Offset: %f %f',fc,offset))
        jumps(fc,:) = offset;        
    end
    currect(1:2) = currect(1:2)-jumps(fc,:);
end

%% use jumps to set cell positions
mincheck = min(find(opt.checkforjump));
if (opt.setCellPos)
    for fc = mincheck:length(cellmv)        
        cellmv(fc).pos = cellmv(fc-1).pos + jumps(fc,:);
    end    
end
end