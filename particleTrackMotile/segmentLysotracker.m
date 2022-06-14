function [cellmvout,opt] = segmentLysotracker(cellmv,framelist,savefile,options)

% get APPROXIMATE cell segmentation based on thresholding lysotracker
% fluorescence
% this is mostly just for tracking COM, not serious measurements of cell
% shape
% segmentation based on thresholding

%% default options
opt = struct();

% dilation to pick out individual cells
opt.bigdilate = 30;
% dilation for smoothing
opt.smalldilate=3;
% line dilation at end to get rid of small invaginations
opt.linedil = 50;
opt.linerfrac=0.9; % fraction of line to erode
% show images as we go
opt.dodisplay = 1;
% how often to save to file (after so many frames)
opt.saveevery = 10;

% number of degrees of freedom for cubic spline smoothing
opt.ndeg = 50; 

% connect disconnected components using convex hull
opt.useconvexhull = 0;
% when calculating the alpha shape around the image, this is the radius of
% the ball to roll around
% in principle can do this after connecting components with convex hull
% (though somewhat redundant)
% if alpharad set to 0 or negative then do not get the alpha shape and just
% use explicit boundary
opt.alpharad = 50
% do line dilation and erosion at the end to join concave regions
% NOTE: this is probably unnecessary if using alpha shape
opt.dolinedilate = 0;

% scale for thresholding relative to the automatic level
opt.threshscl=0.5;
% max out all pixels above this percentile
opt.adjmax=0.99;
% and min out all pixels below this
opt.adjmin=0.01;
% run a bandpass filter to clean up 
% noise in the image before thresholding
% if bpassw is empty then do not run the filter
% otherwise, pair of numbers giving min and max wavelent
opt.bpassw = [];

% which cells to segment; -1 is all of them; numbers above cellcount are
% ignored
opt.whichcells = -1;

% minimal area of cell (in px). throw out anything smaller than this at the
% point of thresholding
opt.minarea = 1e4;
opt.minarea2 = 1e4;

% input options
% only copy over the ones appropriate for this function
inputopt = fieldnames(options);
for c = 1:length(inputopt)
    s = inputopt(c); s=s{1};
    if (isfield(opt,s))
        opt.(s) = options.(s);
    end
end


nframe = length(cellmv);

% variables from options
whichcells = opt.whichcells;
dodisplay = opt.dodisplay;
saveevery= opt.saveevery;

disp([opt.adjmax, opt.threshscl])

%%
framect = 0;
for fc = framelist
    disp(sprintf('Segmenting cells in frame %d',fc));
    
    framect = framect+1;    
    %% 
    img = imread(cellmv(fc).fname,cellmv(fc).fpage);    
    
    %I2 = imadjust(img);
    %% increase contrast by maxing out top % of pixels
    imgn=double(img);
    % normalize image so that stretchlim works properly
    imgn = imgn/max(imgn(:)); 
    I1 = imadjust(imgn,stretchlim(imgn,[opt.adjmin,opt.adjmax]),[0,1]);
    %I2=img;
    %
    if (~isempty(opt.bpassw))
        padn=floor(opt.bpassw(2)/2);
        I2 = padarray(I1,[padn,padn]);
        I2 = bpass(I2,opt.bpassw(1),opt.bpassw(2));
        crect = [padn+1,padn+1,size(img,2)-1,size(img,1)-1];
        I2 = imcrop(I2,crect);      
    else
        I2 = I1;
    end
    if (dodisplay>=2)
        imshow(I2,[],'initialmagnification','fit');
    end
    %%
    [level, EM] = graythresh(I2);
    bw = im2bw(I2,level*opt.threshscl);
    if (dodisplay>=2)
    imshow(bw,[],'initialmagnification','fit')
    end
    %%
    bw = bwareaopen(bw, opt.minarea);
    
    if (dodisplay>=2)
        imshow(bw,[],'initialmagnification','fit')
    end

    %% dilate using distance to nearest nonzero pixel
    % lump together everything within a certain distance cutoff pixels together
    bw2 = bwdist(bw) <= opt.bigdilate;
    %bw2 = bwareaopen(bw2, opt.minarea2);
    
    if (dodisplay>=2)
    imshow(bw2,'initialmagnification','fit')
    end
    
    %% pick out connected components
    concomp = bwconncomp(bw2, 8);
    
    %%
    % slightly dilated version of image
    bwdil = bwdist(bw)<opt.smalldilate;
   
    if (dodisplay>=2)
    subplot(1,2,1)
    imshow(bw)
    subplot(1,2,2)
    imshow(bwdil)
    end
    
    %% go through each individual cell
    cellmv(fc).segcont = {};
    if (whichcells<0)
        celllist = 1:concomp.NumObjects;
    else
        celllist = whichcells(whichcells<=concomp.NumObjects);
    end
    %%
    cellct = 0;
    for cc = celllist
        %%
        imgcell = zeros(size(img));
        inds = concomp.PixelIdxList{cc};
        % keep only those parts of cluster that are in slightly dilated image
        i2 = find(bwdil(inds));
        
        % throw out components that are too small
        if (length(i2)<opt.minarea2)
            continue
        end
        cellct = cellct+1;        
        
        inds = inds(i2);

        imgcell(inds) = 1;
        if (dodisplay>=2)
            imshow(imgcell,'initialmagnification','fit')
        end
       
         
        %% break cell up into connected components and get their boundaries
        [bounds,labcell] = bwboundaries(imgcell,'noholes');
        %labcell = bwlabel(imgcell1);
        % make boundaries clockwise
        for bc = 1:length(bounds)
            [X,Y] = poly2ccw(bounds{bc}(:,2),bounds{bc}(:,1));
            bounds{bc} = [X,Y];
        end
        %%
        if (opt.useconvexhull)
            % define boundaries based on convex hull that spans
            % disconnected components, then replace with true boundaries
            % within connected components
            %% get the convex hull
            imsize = size(img);
            [Y,X] = ind2sub(imsize,inds);
            
            ch = convhull(X,Y);
            
            chbound = [X(ch) Y(ch)];
            
            if (dodisplay>=2)
                imshow(imgcell,'initialmagnification','fit')
                hold all
                plot(chbound(:,1),chbound(:,2),'r-')
                hold off
            end
        
               
            %% go around the convex hull and find segments that span different connected
            % components
            chind = sub2ind(imsize,chbound(:,2),chbound(:,1));
            
            dch = diff(labcell(chind)); % differences in connected components
            
            dchind = find(dch);
            
            if (dodisplay>=2)
                imshow(imgcell,'initialmagnification','fit')
                hold all
                plot(chbound(:,1),chbound(:,2),'r-')
                plot(chbound(dchind,1),chbound(dchind,2),'c*')
                plot(chbound(dchind+1,1),chbound(dchind+1,2),'y*')
                hold off
            end
        
            %% for each convex hull span within a constant component, replace the hull points with actual boundary
            if (isempty(dchind) )
                fullbound=[];
                for bc = 1:length(bounds)
                    fullbound = [fullbound; bounds{bc}];
                end
                [X,Y] = poly2cw(fullbound(:,1),fullbound(:,2));
                fullbound=[X Y];
                if (dodisplay>=2)
                    hold all
                    %plot(chbound(indstart(bc):indend(bc),1),chbound(indstart(bc):indend(bc),2),'g-')
                    plot(fullbound(:,1),fullbound(:,2),'c-')
                    hold off
                end
            else
                indstart = [1,dchind'+1];
                indend = [dchind',size(chbound,1)];
                
                fullbound = [];
                for bc = 1:length(indstart)
                    % which component is involved for this segment
                    regnum = labcell(chind(indstart(bc)));
                    
                    % find nearest points on boundary
                    diffs = bsxfun(@minus, bounds{regnum}',chbound(indstart(bc),:)');
                    [a,ind1] = min(sum(diffs.^2,1));
                    diffs = bsxfun(@minus, bounds{regnum}',chbound(indend(bc),:)');
                    [a,ind2] = min(sum(diffs.^2,1));
                    
                    if (ind2<ind1)
                        brend = size(bounds{regnum},1);
                        allind = [ind1:brend,1:ind2];
                    else
                        allind = ind1:ind2;
                    end
                    
                    if (dodisplay>=2)
                        hold all
                        %plot(chbound(indstart(bc):indend(bc),1),chbound(indstart(bc):indend(bc),2),'g-')
                        plot(bounds{regnum}(allind,1),bounds{regnum}(allind,2),'c-')
                        hold off
                    end
                    
                    if (bc<length(indstart))
                        fullbound = [fullbound;bounds{regnum}(allind,:); chbound(indend(bc)+1,:); chbound(indstart(bc+1),:)];
                    else
                        fullbound = [fullbound;bounds{regnum}(allind,:)];
                    end
                end
            end
        end
        
        if (opt.alpharad>0)
            % fill holes in mask
            imgcell = imfill(imgcell);
            
            % define boundaries based on alpha shape of thresholded region
            % get pixels within cell region
            [X,Y] = meshgrid(1:size(imgcell,2),1:size(imgcell,1));
            X = X(find(imgcell(:))); Y = Y(find(imgcell(:)));
            if (dodisplay>=3)
                imshow(imgcell,[],'initialmagnification','fit');
                hold all
                plot(X(:),Y(:),'r.')
                hold off
            end
            %%
            [V,S] = alphavol([X Y],opt.alpharad,1);
            %
            fullbound = [X(S.bnd(:,1)) Y(S.bnd(:,1))];
            %fullbound(:,1) = X(S.bnd(:,1),:);     
            %fullbound(:,2) = Y(S.bnd(:,1),:);
             if (dodisplay>=3)
                imshow(imgcell,[],'initialmagnification','fit');
                hold all
                plot(fullbound(:,1),fullbound(:,2),'r-')
                hold off
             end     
        end
           
        if (~opt.useconvexhull & opt.alpharad<=0)
            %% set basic boundary if using neither convex hull or alpha
            % shape
            fullbound=[];
            for bc = 1:length(bounds)
                fullbound = [fullbound; bounds{bc}];
            end
            [X,Y] = poly2cw(fullbound(:,1),fullbound(:,2));
            fullbound=[X Y];
            if (dodisplay>=2)
                hold all
                %plot(chbound(indstart(bc):indend(bc),1),chbound(indstart(bc):indend(bc),2),'g-')
                plot(fullbound(:,1),fullbound(:,2),'c-')
                hold off
            end
        end
        
%         %% remove redundant points
%         dfb = diff(fullbound);
%         ndiff = sum(dfb.^2,2);
%         goodind = [find(ndiff>0);length(fullbound)];
%         fullbound = fullbound(goodind,:);

        %% view final image
        if (dodisplay>=2)
            imshow(img,[],'initialmagnification','fit')
            hold all
            plot(fullbound(:,1),fullbound(:,2),'m-')
            hold off
        end

        %% interpolate boundary by arclength
        param = arclenparam(fullbound');
        paramnew = linspace(param(1),param(end),round(param(end)/10));
        paramnew = paramnew(1:end-1);
        
        interpbound = interp1(param,fullbound,paramnew,'linear');
        interpbound = round(interpbound);
        if (dodisplay>=2)
            imshow(img,[],'initialmagnification','fit')
            hold all
            plot(interpbound(:,1),interpbound(:,2),'m.-')
            hold off
        end

        if (opt.dolinedilate)
            %% line dilation of boundary
            imgbound = zeros(size(img));
            inds = sub2ind(size(img),interpbound(:,2),interpbound(:,1));
            imgbound(inds) = 1;
            %imshow(imgbound,[],'initialmagnification','fit')
            
            
            se90 = strel('line', opt.linedil, 90);
            se0 = strel('line', opt.linedil, 0);
            imgbounddil = imdilate(imgbound, [se90 se0]);
            
            if (dodisplay>=2)
                imshow(imgbounddil,'initialmagnification','fit')
            end
            
            %% line erosion
            se90 = strel('line', opt.linedil*opt.linerfrac, 90);
            se0 = strel('line', opt.linedil*opt.linerfrac, 0);
            imgbounder = imerode(imgbounddil, [se90 se0]);
            if (dodisplay>=2)
                imshow(imgbounder,'initialmagnification','fit')
            end
            
            % one more dilation to fill in small gaps
            se90 = strel('line', 2, 90);
            se0 = strel('line', 2, 0);
            imgbounder = imdilate(imgbounder, [se90 se0]);
            
            finbound = bwboundaries(imgbounder,'noholes');
            finbound = fliplr(finbound{1});
            %         %% fill interior holes
            imgf = imfill(imgbounder, 'holes');
            if (dodisplay>=2)
                imshow(imgf,'initialmagnification','fit')
            end
        else
            finbound = fullbound;
        end
        %% smooth boundary with cubic splines
        param = arclenparam(finbound');
        paramnew = linspace(param(1),param(end),opt.ndeg+1);
        paramnew = paramnew(1:end-1);
        
        smoothbound = interp1(param,finbound,paramnew,'spline');
        smoothbound = [smoothbound; smoothbound(1,:)];
        
        if (dodisplay>1)
            imshow(img,[],'initialmagnification','fit')
            hold all
            %plot(fullbound(:,1),fullbound(:,2),'m-')
            plot(smoothbound(:,1),smoothbound(:,2),'g-')
            hold off
            title(sprintf('Frame %d, Cell %d', fc,cellct))
            drawnow
        end
        
        cellmv(fc).segcont{cellct} = smoothbound;                             
    end
    
    if (dodisplay)
        cmap = cool(length(cellmv(fc).segcont));        
        imshow(I1,[],'initialmagnification','fit')
        for cellct = 1:length(cellmv(fc).segcont)    
            smoothbound = cellmv(fc).segcont{cellct};
            hold all            
            plot(smoothbound(:,1),smoothbound(:,2),'Color',cmap(cellct,:))
            hold off
            title(sprintf('Frame %d, Cell %d', fc,cellct))
            drawnow
        end
    end
    
    if (~isempty(savefile) & mod(framect,saveevery)==0)
        save(savefile)
    end
end

cellmvout = cellmv;