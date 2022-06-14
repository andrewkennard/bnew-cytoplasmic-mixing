% class definition for a cell object, containing data on particle
% positions, segmentation, etc.

classdef CellObj < handle
%     properties (SetAccess = immutable)
%         % cell name
%         Name = 'cell';  
%     end
    
    properties (SetAccess = private)
        % cell contour spline for interpolation
        SpCont = [];
        % set spline for interpolating contour center of mass
        SpCOM = [];
    end
    
    properties                               
        % Cell Name
        Name = 'cell';
        
        % Full names of movie files associated with this cell
        MovieFiles = {};
        
        % image size
        ImgSize = [];
        
        % directory where (input) data files for this cell are stored
        DirName = '';
        % filename originally used to load data into cell object
        LoadFile='';  
        
        % number of frames
        NFrame = 0;
        % structure array containing data on individual frames
        % fields are:
        % fname: movie file name 
        % fpage: page within the movie file
        % eltime: elapsed time directly from metadata
        % time: time in sec associated with each frame
        % bad: marker of bad frames
        % pos: stage position approximated by image alignment
        % rect: cropping rectangle containing this cell within the frame
        % meanbr: mean brightness
        % segcont: segmented cell contour
        % COM: cell contour center of mass
        % ptdata = nx6 array of data on found (not tracked) points
        % columns are: x,y,mass,rg,eccentricity
        % pts: single structure with fields:
        % xy = nx2 array of particle coordinates within the frame
        % pi = corresponding indices in Particles array
        % pic = corresponding indices within each individual particle track
        % np = number of particles in this cell, in this frame
        Frames = [];         
        
        % number of particles
        NParticle = 0;
        % structure array containing data on individual particles
        % fields are:
        % xy = x,y coordinates over multiple frames
        % fi = corresponding frame indices
        % fic = index of each particle in the tracks array for each frame
        % mass = integrated feature mass
        % rg2 = square of radius of gyration
        % ecc = feature eccentricity
        % trackopt = options used for linking tracks
        Particles = [];     
        % options used for linking tracks
        TrackOpt = struct();
        
        % Wavelet analysis object
        Wavelet = [];
        
        % avg signal to noise
        sig2noise = 0;                
    end
    
    methods
        function CL = CellObj(name)
            % create a cell object with the given name
            CL.Name = name;
            
            % set up empty Frames array
            CL.Frames = struct('fname',{},'fpage',{},'eltime',{},'time',{},'bad',{},'pos',{},'rect',{},'meanbr',{},'segcont',{},'COM',{},'tracks',[]);
            % set up empty Particles array
            CL.clearParticles;   
                       
        end        
        
        function setSpCont(CL)
            % create a cubic spline object for interpolating cell contours
            % over time
            % uses data in the Frames.segcont array
            % sets the SpCont field
                                   
            valframes =find(arrayfun(@(x) ~isempty(x.segcont),CL.Frames));
            
            if (isempty(valframes))
                error('Cannot set contour splines because segcont is not set')
            end
            
            times = [CL.Frames.time];
            valtimes = times(valframes);
            
            % segmentation contour spline, in relative coordinates            
            segcont = cat(3,CL.Frames(valframes).segcont);  
            for sc = 1:size(segcont,3)
                segcont(:,:,sc) = bsxfun(@minus, segcont(:,:,sc), CL.Frames(valframes(sc)).pos);
            end
            CL.SpCont = spline(valtimes,segcont);
            
            % center of mass spline, in relative coordinates
            COM = cat(1,CL.Frames(valframes).COM); 
            COM = COM-vertcat(CL.Frames(valframes).pos);
            CL.SpCOM = spline(valtimes,COM');            
        end
        
        function setRect(CL,expf)
            % for each frame in the cell, set the rectangle approximately
            % containing the cell contour
            % for frames that do not have the cell contour defined, use
            % interpolation
            % expand by some number of pixels in both directions
            % expansion border is given by factor expf-1 (expf>1 for positive border) times
            % sqrt(cell area)           
            
            % make contour spline if not already done
            if (isempty(CL.SpCont))
                CL.setSpCont();
            end
            
            cutx = [1,CL.ImgSize(2)]; 
            cuty = [1,CL.ImgSize(1)];
            
            for fc = 1:length(CL.Frames)
                %%
                frame = CL.Frames(fc);
                
                % interpolate contour
                cellcont = fnval(CL.SpCont,frame.time);
                % absolute coordinates
                cellcont = bsxfun(@plus,cellcont,frame.pos);
                
                % rectangle containing contour
                xmin = min(cellcont(:,1));
                xmax = max(cellcont(:,1));
                ymin = min(cellcont(:,2));
                ymax = max(cellcont(:,2));
                rect = [xmin,ymin,xmax-xmin+1,ymax-ymin+1];
                
                % pad around the contour
                tmp = polygeom(cellcont(:,1),cellcont(:,2));
                padn = sqrt(tmp(1))*(expf-1);
                xmin = max(cutx(1),xmin-padn);
                ymin = max(cuty(1),ymin-padn);
                xmax = min(cutx(2),xmax+padn);
                ymax = min(cuty(2),ymax+padn);
                exprect = [xmin,ymin,xmax-xmin+1,ymax-ymin+1];
                %%
                %exprect = expandRect(rect,expf,cutx,cuty);
                CL.Frames(fc).rect = exprect;
            end                        
        end
        
        function [img] = viewFrame(CL,fc,varargin)
            % crop to cell rectangle
            docrop=1;
            % show particles
            particles=1;
            % show segcont
            segcont = 1;
            % contrast limits in fractions of full range
            contrast = [0,1];
            % show results from PIV analysis of flow
            PIV = 0;
            pivscl=10; % scaling of arrows in showing piv
            % show rectangle
            showrect = 0;
            % rectangle to crop to (default is the rectangle associated
            % with the frame)
            rect = [];
            % initial magnification
            initialmag = 'fit';
            % show cell frame of reference position
            showcellpos = 0;
            % show COM
            showCOM = 0;
            % arguments to imshow
            imshowargs = {'initialmagnification','fit'};
            % show particle number index in pts array
            particlenum = 0;
            for vc = 1:2:length(varargin)
                switch(varargin{vc})
                    case('docrop')
                        docrop=varargin{vc+1};
                    case('particles')
                        particles=varargin{vc+1};
                    case('segcont')
                        segcont = varargin{vc+1};
                    case('contrast')
                        contrast = varargin{vc+1};
                    case('rect')
                        rect = varargin{vc+1};
                    case('showrect')
                        showrect = varargin{vc+1};
                    case('PIV')
                        PIV = varargin{vc+1};
                    case('pivscl')
                        pivscl = varargin{vc+1};
                    case('initialmag')
                        initialmag = varargin{vc+1};
                    case('showcellpos')
                        showcellpos = varargin{vc+1};
                    case('showCOM')
                        showCOM = varargin{vc+1};
                    case('imshowargs')
                        imshowargs = varargin{vc+1};
                    case('particlenum')
                        particlenum = varargin{vc+1};
                end
            end
            
            % make contour spline if not already done
            if ((segcont && isempty(CL.SpCont)) || (showCOM && isempty(CL.SpCOM)))
                CL.setSpCont();
            end
            
            %% view one specific frame of the cell
            
            frame = CL.Frames(fc);
            if (isempty(rect))
                rect = frame.rect;
            end
            img = imread(frame.fname,frame.fpage);
                        
            brightrange = [prctile(img(:),contrast(1)*100),prctile(img(:),contrast(2)*100)];
            
            if (segcont)
                segcont = fnval(CL.SpCont,frame.time);
                segcont = bsxfun(@plus, segcont, frame.pos);
            end
            
            if (docrop)
                imgC = imcrop(img,rect);
                if (~isequal(contrast,[0,1]))
                    brightrange = [prctile(imgC(:),contrast(1)*100),prctile(imgC(:),contrast(2)*100)];
                else
                    brightrange = [];
                end
                imshow(imgC,brightrange,'initialmagnification',initialmag); 
                if (segcont)
                    segcont = bsxfun(@minus,segcont,rect(1:2));
                end
            else
                if (~isequal(contrast,[0,1]))
                    brightrange = [prctile(img(:),contrast(1)*100),prctile(img(:),contrast(2)*100)];
                else
                    brightrange = [];
                end
                if (isempty(initialmag))
                     imshow(img,brightrange);
                else
                    imshow(img,brightrange,'initialmagnification',initialmag);
                end
                %segcont = frame.segcont;
            end
            
            hold all   
            if (segcont)
                plot(segcont(:,1),segcont(:,2),'c')
            end
            if (particles & ~isempty(frame.pts.xy)>0)      
                if (docrop)
                    plot(frame.pts.xy(:,1)-rect(1)+1,frame.pts.xy(:,2)-rect(2)+1,'r.')
                else
                    plot(frame.pts.xy(:,1),frame.pts.xy(:,2),'r.')
                end
                if (particlenum)
                    if (docrop)
                    error('particlenum with docrop not set up yet')
                    else
                        hold all
                        for pc = 1:size(frame.pts.xy,1);
                            text(frame.pts.xy(pc,1),frame.pts.xy(pc,2),sprintf('%d',pc),'Color','y')
                        end
                        hold off
                    end
                end
            end            
            hold off            
            
            if (showrect)
                hold all
                rectangle('Position',rect,'EdgeColor','y')
                hold off
            end
            
            hold all
            text(10,10,sprintf('Frame %d',fc),'Color','y')
            hold off
            
            if (PIV)
                if (~isfield(CL.Frames(fc),'PIV') || isempty(CL.Frames(fc).PIV))
                    warning('No PIV results found')
                else
                    PIV = frame.PIV;
                    if (docrop)
                        PIV(:,1) = PIV(:,1)-rect(1)+1;
                        PIV(:,2) = PIV(:,2)-rect(2)+1;
                    end
                    hold all
                    quiver(PIV(:,1),PIV(:,2),PIV(:,3)*pivscl,PIV(:,4)*pivscl,0,'g')
                    hold off
                end
                
                
            end
            
            if (showcellpos)
                if (docrop)
                    pos = CL.Frames(fc).cellpos-rect(1:2)+1;
                else
                    pos = CL.Frames(fc).cellpos;
                end
                hold all
                plot(pos(1),pos(2),'g*')
                hold off
            end
            
            if (showCOM)
                COM = fnval(CL.SpCOM,CL.Frames(fc).time);
                COM = COM+CL.Frames(fc).pos';
                if (docrop)
                    COM = COM'-rect(1:2)+1;
                end
                hold all
                plot(COM(1),COM(2),'m*')
                hold off
            end
        end
        
        function M = viewTrack(CL,tc,varargin)
            % view a specific track over time
            
            % crop to cell rectangle
            docrop = 0;
            % show cell frame of reference
            showcellpos = 0;
            % which frame indeces to show
            % default is all of them
            frameind = [];
            for vc = 1:2:length(varargin)
                switch(varargin{vc})
                    case('docrop')
                        docrop = varargin{vc+1};
                    case('showcellpos')
                        showcellpos = 1;
                    case('frameind')
                        frameind = [];
                end
            end
            
            %%
            track0 = CL.Particles(tc).xy;           
            
            if (isempty(frameind))
                frameind = 1:size(track0,1);
            end
            ct = 0;
            for cc = frameind%1:size(track0,1)
                fc = CL.Particles(tc).fi(cc);
                CL.viewFrame(fc,'particles',0,'docrop',docrop,'showcellpos',showcellpos);
                rect = CL.Frames(fc).rect;
                
                if (docrop)
                    trackrel = bsxfun(@plus,track0,-rect(1:2)+1);
                else
                    trackrel = track0;
                end
                
                trackrel = bsxfun(@plus,trackrel,CL.Frames(fc).pos);
                
                hold all
                plot(trackrel(1:cc,1),trackrel(1:cc,2),'r-')
                plot(trackrel(cc,1),trackrel(cc,2),'c.','MarkerSize',10)
                hold off
                               
                drawnow
                
                if (nargout>0)
                    ct=ct+1;
                    M(ct) = getframe(gcf);
                end
            end
        end
        
        function [COM,M] = viewCOMmovie(CL,varargin)
            % crop to cell rectangle
            docrop=1;
            % which frames to show in movie
            frames = 1:length(CL.Frames);
            % show cell frame of reference fiduciary point rather than
            % segmented COM
            cellpos = 0;
            cellframe=0;
            % show segcont
            segcont = 1;
            % save avi movie to file
            moviefile = '';
            %magnification
            initialmag = '';
            
            for vc = 1:2:length(varargin)
                switch(varargin{vc})
                    case('docrop')
                        docrop=varargin{vc+1};   
                    case('frames')
                        frames = varargin{vc+1};
                    case('cellpos')
                        cellpos = varargin{vc+1};
                    case('cellframe')
                        cellframe = varargin{vc+1};
                    case('segcont')
                        segcont = varargin{vc+1};
                    case('moviefile')
                        moviefile = varargin{vc+1};
                    case('initialmag')
                        initialmag = varargin{vc+1};
                end
            end
            
            
            goodf = find(arrayfun(@(x) ~isempty(x.pos),CL.Frames))
            frames = frames((frames<=max(goodf)));
            
            % open movie file for writing
            if (~isempty(moviefile))
                figpos = [150 150 fliplr(CL.ImgSize)];
                set(gcf,'Position',figpos);    
               % M(step) = getframe(gcf);
            end
    
            %    writerObj = VideoWriter(moviefile);
            %    framerate = 1/mean(diff([CL.Frames.time]))
            %    writerObj.FrameRate = framerate;
            %    open(writerObj);
            %end
            
            % make contour spline if not already done
            if (isempty(CL.SpCont))
                CL.setSpCont();
            end
            
            %%
            tlist = [CL.Frames(frames).time];
                
            
                rect = vertcat(CL.Frames(frames).rect);
            
            
            if (cellpos || cellframe)
                % cell frame of reference positions
                cellpos = vertcat(CL.Frames(frames).cellpos);
            end
            
            if (cellframe) % set rectangle in cell's frame of reference
                maxwidth = max(round(rect(:,3)));
                if (mod(maxwidth,2)==1); maxwidth = maxwidth+1; end
                maxheight = max(round(rect(:,4)));                
                if (mod(maxheight,2)==1); maxheight = maxheight+1; end
                wh = ones(size(cellpos));                
                wh(:,1) = maxwidth; wh(:,2) = maxheight;                
                rect = [cellpos(:,1)-maxwidth/2, cellpos(:,2)-maxheight/2, wh];   
                
                badind= find(rect(:,1)<0)
                rect(badind,3) = rect(badind,3)+rect(badind,1)-1;
                rect(badind,1) = 1;
                
                
                badind = find(rect(:,2)<0)
                rect(badind,4) = rect(badind,4)+rect(badind,2)-1;
                rect(badind,2) = 1;                                
            end
            
            if (cellpos)                
                COM = cellpos;
                 if (docrop)
                    COM = COM-rect(:,1:2);
                end
            else
                % relative COM positions
                COM = fnval(CL.SpCOM,tlist)';
                pos = vertcat(CL.Frames(frames).pos);
                if (docrop)
                    COM = COM+pos-rect(:,1:2);
                else
                    COM = COM+pos ;
                end
            end
            
            
            %%
            clear M
            for ct = 1:length(frames)
                fc = frames(ct);
                % show each frame                   
                CL.viewFrame(fc,'particles',0,'docrop',docrop, 'rect',rect(ct,:),'segcont',segcont,'initialmag',initialmag);
                hold all                
                plot(COM(ct,1),COM(ct,2),'r*')
                hold off
                drawnow    
                
                if (~isempty(moviefile))                    
                    %set(gcf,'Position',[1 1 rect(ct,3:4)])
                    %set(gca,'Position',[0 0 0.95 0.95])
                    %fr = getframe(gcf);
                    %writeVideo(writerObj,fr);
                    M(ct) = getframe(gcf);
                end
            end
           
            %if (~isempty(moviefile))                
               % close(writerObj);
            %end
        end
        
        function clearParticles(CL,varargin)
            
           clearptdata = 1;
           for vc = 1:2:length(varargin)
               switch(varargin{vc})
                   case('clearptdata')
                       clearptdata = varargin{vc+1};
               end               
           end
           
           %% clear particle lists and tracks for this cell
            CL.NParticle = 0;
            CL.Particles = struct('xy',{},'fi',{},'fic',{});            
            for fc = 1:length(CL.Frames)
                CL.Frames(fc).pts = struct('xy',[],'pi',[],'pic',[],'np',0);
                if (clearptdata)
                % data on all points found (including those not in tracks)                
                CL.Frames(fc).ptdata = [];
                end
            end            
        end
        
        function tracklist = getTracklist(CL,varargin)
            % tracks relative to cell COM
            % otherwise, given relative to stage position
            relCOM = 0;
            relCellPos = 0;
            % frames that are excluded from tracks (and tracks broken up)
            badframes = [];
            % split tracks at every skipped frame
            splitskip = 0;
            % linearly interpolate through missing frames
            interptracks = 0;
            % minimal track length to keep
            mintracklen=2;
            % usegood
            % 0: only particles/frames not marked good
            % 1: only those marked good
            % ignore good marking
            usegood = 2;
            % which particle indices to include
            particleind = 1:length(CL.Particles);
            
            % extend tracks by the given number of points on each side
            extendtracks = 0;           
            
            for vc = 1:2:length(varargin)
                switch (varargin{vc})
                    case('relCOM')
                        relCOM = varargin{vc+1};
                    case('relCellPos')
                        relCellPos = varargin{vc+1}; 
                    case('badframes')
                        badframes = varargin{vc+1};
                    case('splitskip')
                        splitskip = varargin{vc+1};
                    case('interptracks')
                        interptracks = varargin{vc+1};
                    case('mintracklen')
                        mintracklen = varargin{vc+1};
                    case('usegood') 
                        usegood = varargin{vc+1};
                    case('particleind')
                        particleind = varargin{vc+1};       
                    case('extendtracks')
                        extendtracks = varargin{vc+1};                   
                end
            end
            
            % COM gives the COM of cell frame position (relative to stage position)   
            if (relCOM)
                % use tracks relative to cell COM
                COMpos = fnval(CL.SpCOM,[CL.Frames.time])';
                %minframe = min(arrayfun(@(x) x.xy(1,6),CL.Particles));
                %maxframe = min(arrayfun(@(x) x.xy(1,6),CL.Particles));                
            elseif (relCellPos)
                COMpos = vertcat(CL.Frames.cellpos)-vertcat(CL.Frames.pos);
            end
            
            %% get tracklist dictionary for cell
            tracklist = {};
            ct=0;
            for tc = particleind%1:CL.NParticle
                %if (length(CL.Particles(tc).xy)>2*max(nvals))
                    ct =ct+1;
                    if (relCOM || relCellPos)
                        tracklist{ct} = CL.Particles(tc).xy - COMpos(CL.Particles(tc).fi,:);
                    else
                        tracklist{ct} = CL.Particles(tc).xy;
                    end
                    tracklist{ct}(:,6) = CL.Particles(tc).fi;
                    tracklist{ct}(:,8) = tc;
                    if (usegood==0)
                        goodind = find(~CL.Particles(tc).good);
                    elseif (usegood==1)
                        goodind = find(CL.Particles(tc).good);
                    else
                        goodind = 1:size(CL.Particles(tc).xy,1);
                    end
                    tracklist{ct} = tracklist{ct}(goodind,:);
                %end
            end            
            
            tracklens = cellfun(@(x) size(x,1), tracklist);
            ind = find(tracklens>=mintracklen);
            tracklist = tracklist(ind);
            
            tracklist = breakupTracks(tracklist,badframes);
            
            if (splitskip) % break up tracks at missing frames
                tracklist = splitTrackSkip(tracklist);
            end
            
            tracklens = cellfun(@(x) size(x,1), tracklist);
            ind = find(tracklens>=mintracklen);
            tracklist = tracklist(ind);
            
            if (interptracks) % linearly interpolate missing frames
                for tc = 1:length(tracklist)
                    tracklist{tc} = interpolateTrack(tracklist{tc});
                end
            end
            
            % extend tracks
            if (extendtracks>0)
                for tc = 1:length(tracklist)
                    track = tracklist{tc};                   
                    
                    %% do a linear interpolation to get edge
                    % velocities                    
                    framewant1 = track(1,6)-(1:extendtracks);
                    framewant2 = track(end,6)+(1:extendtracks);
                    framewant = [framewant1 framewant2];
                    extratrack = interp1(track(:,6),track(:,1:2),framewant,'linear','extrap');
                    newtrack = zeros(2*extendtracks,size(track,2));
                    newtrack(:,1:2) = extratrack;
                    newtrack(:,6) = framewant;
                    newtrack(:,8) = track(1,8);
                    newtrack = [flipud(newtrack(1:extendtracks,:)); track; newtrack(extendtracks+1:end,:)];
                    
                    tracklist{tc} = newtrack;                   
                end
            end
        end
        
        function [tracklist,BN] = runBNEWAnalysis(CL,nvals,varargin)
            % options for fitting coefficients
            fitopt = struct();
            % file storing rescaling coefficients
            % if empty, then recalculate from scratch
            coeffile = '';
            % tracks relative to cell COM
            relCOM = 0;
            relCellPos =0;
            % only keep portions of tracks within these frame limits
            framelim = [];
            %
            timeavg=1;
            
            for vc = 1:2:length(varargin)
                switch varargin{vc}
                    case('fitopt')
                        fitopt = varargin{vc+1};                            
                    case('coeffile')
                        coeffile = varargin{vc+1};
                    case('relCOM')
                        relCOM = varargin{vc+1};
                    case('relCellPos')
                        relCellPos = varargin{vc+1};
                    case('framelim')
                        framelim = varargin{vc+1};
                    case('timeavg')
                        timeavg = varargin{vc+1};                    
                end
            end                       
            
            tracklist = CL.getTracklist(varargin{:});   
            
            
            ntrack = length(tracklist);
            
            BN = BNEWobj(nvals,varargin{:});
            
            % get wavelet coefficients
            BN = BN.getCoefficients();
            
            % run wavelet analysis on
            tic
            BN = BN.analyzeTracks(tracklist,'timeavg',timeavg);
            toc
            
            %% rescale to collapse everything to universal line, using A^n_k and B^n_k functions
            if (isempty(coeffile))
                BN = BN.rescaleData('recalc',1);
            else
                BN = BN.loadWaveletCoeff(coeffile);
                BN = BN.rescaleData('recalc',0);
            end

            [BN,Dfit] = BN.fitDcoeff(fitopt);
            CL.Wavelet = BN;
% 
%             
%             %% run wavelet analysis on all track data for this cell
%             % using half-spans nvals              
%             
%             % make a wavelet analysis object
%             WL = WaveletAnalysisObj(nvals,varargin{:});
% 
%             % get wavelet coefficients
%             WL = WL.getCoefficients();
%             
%             % run wavelet analysis
%             WL = WL.analyzeTracks(tracklist,'timeavg',timeavg);
% 
%             %% rescale everything
%             tic
%             if (isempty(coeffile))
%                 WL = WL.rescaleData('recalc',1);
%             else
%                 WL = WL.loadWaveletCoeff(coeffile)
%                 WL = WL.rescaleData('recalc',0);
%             end
%             toc
%             
%             %% fit coefficients
%             WL = WL.fitDcoeff(fitopt);
%             
%             % attach wavelet object to cell
%             CL.Wavelet=WL;
        end
        
        function [rmsvel, interpCOM,tlist] = getCOMvel(CL,del,varargin)
            % get approximate center of mass velocities 
            % using direct difference of spline-interpolated COM (in
            % absolute coords)
            % del = time difference to use for velocities
            % WARNING: interpolation over stage shifts may be subpar
            % should recode eventually to do interpolation over relative
            % coordinates
            
            % maximal jump allowed for COM
            maxjump = inf;
            
            for vc = 1:2:length(varargin)
                switch (varargin{vc})
                    case('maxjump')
                        maxjump = varargin{vc+1};
                end
            end
            
            % make contour spline if not already done
            if (isempty(CL.SpCOM))
                CL.setSpCont();
            end
            
            valframes = find(arrayfun(@(x) ~isempty(x.segcont), CL.Frames));
            
            
            tlist = CL.Frames(valframes(1)).time:del:CL.Frames(valframes(end)).time;
            interpCOM = fnval(CL.SpCOM,tlist)';            
            
            diffs = interpCOM(2:end,:)-interpCOM(1:end-1,:);
            ndiffs = sum(diffs.^2,2);
            ndiffs = ndiffs(ndiffs<maxjump^2);
            
            rmsvel = sqrt(mean(ndiffs));
            rmsvel = rmsvel./del;
        end
        
        function partopt = findParticles(CL,options)
            % find particles for each frame of cell
            % assumes alignment to take care of stage jump has already been
            % done
            % Also get average brightness of whole frame and segmented cell
            
            % default options
            opt = struct();
            % expand segmentation contour by a given distance
            % and search only particles within the expanded contour 
            % if empty then search for all particles within preset cell rectangle
            opt.expandcont = [];
            % which frames to search in
            opt.frames = 1:length(CL.Frames);
            % amount of displaying (0 is no figures drawn, >1 does
            % intermediate steps)
            opt.dodisplay = 0;            
            
            % particle finding parameters
            % upper limit width (diameter) of features, used in bandpass filter
            opt.w = 8;
            % intensity ratio cutoff
            % only keep features with a brightness per unit area that's at least this
            % fraction of the maximum (with some top bright outliers cut out to
            % calculate maximum)
            opt.intratscl = 0.5;
            % percentile cutoff for overall feature mass
            opt.ppmass = 80;
            
            % maximal eccentricity
            opt.maxecc = 0.3;
            
            % invert image to find dark particles
            opt.invert = 0;
            % limits for contrast stretching
            % percentiles to min out and max out
            opt.stretchlim=[0,1];

            % save every so many frames
            opt.saveevery = 20;
            % file to save to (empty means do not save)
            opt.savefile = '';
            % save cell to file even if empty
            opt.saveemptycell = 0;            
            
            % brightness scaling for displaying image
            opt.showscale = [];
           
            % type of average to use in brightness calculations
            opt.bravgtype = 'median';
            
            % crop to integer pixel values
            opt.intpixcrop = 1;
            
            % copy over input options
            opt = copyStruct(options, opt);
            
            partopt = opt;
            
            %% assumes spline and cell rectangle have been set             
            if (isempty(CL.SpCont))
                error('Contour spline must be preset to run particle search')                
            end
            if (~isfield(CL.Frames,'rect'))
                error('Cell rectangle must be preset to run particle search')
            end
            
            % relative segmentation contours for all frames
            segcontrel = fnval(CL.SpCont,[CL.Frames.time]);
            
            ct = 0;
            for fc = opt.frames                
                ct = ct+1;
                %%
                frame = CL.Frames(fc);
                img = imread(frame.fname, frame.fpage);
                img0 = img;
                if (opt.invert)
                    img = max(img(:))-img;
                end
                
                % absolute contour
                segcont = bsxfun(@plus,segcontrel(:,:,fc),frame.pos);
                
                [CL.Frames(fc).avgbr,CL.Frames(fc).avgbrseg] = getImgBrightness(img,segcont,'bravgtype',opt.bravgtype);
                
                % expand contour and set rectangle for search
                if (isempty(opt.expandcont))
                    rect = frame.rect;
                else
                    %%
                    % convert to counterclockwise
                    % when visualized will look clockwise since image has y axis flipped
                    [x,y] = poly2ccw(segcont(:,1),segcont(:,2));                    
                    contccw = [x y];                              
                    contexp = expandPolygon(contccw,opt.expandcont);
                    clen = cellfun(@(x) length(x),contexp);
                    [mlen,maxind] = max(clen);
                    contexp = contexp{maxind};
                    contexp = [contexp;contexp(1,:)];
                    
                    minx = max(min(contexp(:,1)),1);
                    maxx = min(max(contexp(:,1)),size(img,2));
                    miny = max(min(contexp(:,2)),1);
                    maxy = min(max(contexp(:,2)),size(img,1));
                    rect = [minx miny maxx-minx maxy-miny];
                    if (opt.intpixcrop)
                        rect = round(rect);                    
                    end
                end
               
                %%
                if (opt.dodisplay>1)
                    %%
                    imshow(img,[],'initialmagnification','fit')
                    hold all
                    plot(segcont(:,1),segcont(:,2),'c')                   
                    
                    if (~isempty(opt.expandcont))
                        plot(contexp(:,1),contexp(:,2),'m')
                    end
                    
                    rectangle('Position',rect,'EdgeColor','y')
                    text(100,100,sprintf('Frame %d',fc),'Color','y')
                    hold off
                    drawnow
                end
                
                %%
                % crop image
                imgC = imcrop(img,rect);
                imgC = double(imgC)/max(double(imgC(:)));
                
                %% adjust image contrast
                imgA = imadjust(imgC,stretchlim(imgC,opt.stretchlim));

                % bandpass filter
                imgclean = bpass(imgA,1,opt.w);    
                maxclean = max(imgclean(:));
                
                % parameters for feature finding
                Imax= 0.1;      
                massmin = prctile(imgclean(:),opt.ppmass)*(2*opt.w)^2;
                massmax = inf;
                dodisplay=0;
                % allow for a bright outlier
                pp = 1-(2*opt.w+1)^2/prod(size(imgclean));
                maxclean = prctile(imgclean(:),pp*100);
                intrat = opt.intratscl*maxclean;
        
                [goodpts,tmp,pts] = detectParticles(imgclean,opt.w,'Imax',Imax,'massmin',massmin,'massmax',massmax,...
                    'maxecc',opt.maxecc,'intrat',intrat);
                
                if (opt.dodisplay>1)
                    imshow(imgC,[])
                    hold all
                    plot(goodpts(:,1),goodpts(:,2),'b.')
                    hold off
                    drawnow
                end
                
                % adjust back to absolute coords
                fnum = ones(size(goodpts,1),1)*(fc);
                if (~isempty(goodpts)) 
                    goodpts(:,1:2) = bsxfun(@plus,goodpts(:,1:2),rect(1:2)-1);                    
                end
                
                %% keep only those points within the expanded contour
                if (~isempty(opt.expandcont) & ~isempty(goodpts))
                    IN = inpolygon(goodpts(:,1),goodpts(:,2),contexp(:,1),contexp(:,2));
                    goodpts = goodpts(IN,:);
                end                                
                
                if (opt.dodisplay & ~isempty(goodpts))
                    imgshow = imadjust(double(img0)/max(double(img0(:))),stretchlim(opt.stretchlim));
                    imshow(imgshow,opt.showscale,'initialmagnification','fit')
                    hold all
                    plot(segcont(:,1),segcont(:,2),'c')                   
                    
                    if (~isempty(opt.expandcont))
                        plot(contexp(:,1),contexp(:,2),'m')
                    end
                    
                    rectangle('Position',rect,'EdgeColor','y')
                    text(100,100,sprintf('Frame %d',fc),'Color','y')
                    plot(goodpts(:,1),goodpts(:,2),'b.')
                    hold off
                    drawnow
                end
                
                %% attach points to frame
                CL.Frames(fc).ptdata = goodpts;
                
                if (~isempty(opt.savefile) & mod(ct,opt.saveevery)==0)
                    save(opt.savefile,'CL','partopt');
                end
                
                % options used to find particles in this frame
                CL.Frames(fc).partopt = partopt;
            end
        end
        
        function [gaussopt] = gaussOptParticles(CL,options)
            % optimize gaussian positions of all particles
            
            % which frames to run
            opt.frames = 1:CL.NFrame;
            
            % expected particle diameter
            opt.w = 8;
            %opt.span = opt.w;
            
            opt.dodisplay = 1;
            
            opt.filtertype = 'bpass'
            
            opt.saveevery = 100;
            opt.savefile = 'tmp';
            
            opt.expandcont = 50;
            %%
            if (nargin>1)
                % input options
                opt = copyStruct(options,opt);
            end                
            
            dodisplay = opt.dodisplay;
            gaussopt = opt;
            
            for fc = opt.frames
                fr = CL.Frames(fc);
                img0 = imread(CL.Frames(fc).fname, CL.Frames(fc).fpage);
                
                segcont = bsxfun(@plus,fnval(CL.SpCont,fr.time),fr.pos);
                if (isempty(opt.expandcont))
                    rect = fr.rect;
                else
                    %%
                    % convert to counterclockwise
                    % when visualized will look clockwise since image has y axis flipped
                    [x,y] = poly2ccw(segcont(:,1),segcont(:,2));                    
                    contccw = [x y];                              
                    contexp = expandPolygon(contccw,opt.expandcont);
                    clen = cellfun(@(x) length(x),contexp);
                    [mlen,maxind] = max(clen);
                    contexp = contexp{maxind};
                    contexp = [contexp;contexp(1,:)];
                    
                    minx = max(min(contexp(:,1)),1);
                    maxx = min(max(contexp(:,1)),size(img0,2));
                    miny = max(min(contexp(:,2)),1);
                    maxy = min(max(contexp(:,2)),size(img0,1));
                    rect = round([minx miny maxx-minx maxy-miny]);                   
                 end
                
                img = imcrop(img0,rect);
                
                img = double(img)/max(double(img(:)));
                
                
                
                medval = median(img(:));
                if strcmp(opt.filtertype,'bpass')
                    % bandpass filter the image
                    imgclean = bpass(img,1,opt.w);
                else
                    imgclean = img;
                end
                newptdata = zeros(size(fr.ptdata,1),5);
                %%
                for pc = 1:size(fr.ptdata,1)
                    %%
                    pt = fr.ptdata(pc,1:2)-rect(1:2)+1;
                    
                    if (pt(1)<opt.w*2 | pt(1)>size(img,2)-opt.w*2 ...
                            | pt(2)<opt.w*2 | pt(2)>size(img,1)-opt.w*2)
                        newptdata(pc,:) = NaN;
                        continue
                    end
                    c0 = [1,(opt.w/2)^2,pt,medval];
                    opt2 = opt;
                    opt2.dodisplay = 0;
                    [cfit, fval,I] = gaussianFit(imgclean,c0,opt2);
                    newptdata(pc,:) = cfit;
                end
                
                newptdata(:,3:4) = bsxfun(@plus,newptdata(:,3:4),rect(1:2)-1);
                
                %% get rid of points where gaussian is too wide
                goodind = find(newptdata(:,2)<(opt.w/2)^2*2);
                goodptdata = fr.ptdata(goodind,:);
                goodptdata(:,1:2) = newptdata(goodind,3:4);
                goodptdata(:,6:8) = newptdata(goodind,[1,2,5]);
                %%
                if (opt.dodisplay)
                    imshow(img0,[]);
                    hold all
                    plot(fr.ptdata(:,1),fr.ptdata(:,2),'r.')
                    plot(goodptdata(:,1),goodptdata(:,2),'b.')
                    title(sprintf('Frame %d of %d',fc, CL.NFrame))
                    drawnow
                end
                
                %%
                CL.Frames(fc).gaussptdata = goodptdata;
                
                %% 
                if (mod(fc,opt.saveevery)==0)
                    save(opt.savefile,'CL','gaussopt');
                end
            end
            
        end
        
        function [tracklist,opt] = linkTracks(CL,options)
            % link up particles found on each frame into tracks
            
            %% default options
            opt = struct();
            % allowed to skip so many frames in linking
            opt.memory = 3;
            % minimal length of tracks
            opt.goodenough=40;
            % maximal particle displacement in pixels
            opt.maxdisp=6;
            
            % use gaussian optimized points
            opt.usegausspts = 0;
            
            if (nargin>1)
                % input options
                opt = copyStruct(options,opt);
            end                       
            
            %% create allpts list of all particles from all frames
            % convert to relative coords (relative to stage shift)
            
            allpts = [];
            for fc = 1:CL.NFrame;
                if (opt.usegausspts)
                    if isempty(CL.Frames(fc).gaussptdata)
                        continue
                    end
                    pts = CL.Frames(fc).gaussptdata(:,1:5);
                else
                    pts = CL.Frames(fc).ptdata;
                end
                if (isempty(pts) | sum(pts(3:5))==0) % put placeholder particles in empty frames
                    adjpts= [1, zeros(1,4), fc];
                else
                    adjpts = [pts,fc*ones(size(pts,1),1)];
                    adjpts(:,1) = pts(:,1)-CL.Frames(fc).pos(1);
                    adjpts(:,2) = pts(:,2)-CL.Frames(fc).pos(2);
                end
                allpts = [allpts; adjpts];
            end

            %% get particle tracks by linking particles
            dim = 2;
            
            % pad with last column of zeros
            npt = length(allpts);
            xyzs = [allpts,zeros(npt,1)];
            
            % shift positions to avoid negative coordinates since this breaks the
            % tracking code
            xshift=min(xyzs(:,1));
            yshift = min(xyzs(:,2));
            xyzs(:,1) = xyzs(:,1)-xshift+1;
            xyzs(:,2) = xyzs(:,2)-yshift+1;
            
            [lub] = trackmem(xyzs,opt.maxdisp,dim,opt.goodenough,opt.memory);
            lub(:,1) = lub(:,1)+xshift-1;
            lub(:,2) = lub(:,2)+yshift-1;
            
            % get rid of place holder particles
            lub = lub(lub(:,1)~=1 | lub(:,2)~=0,:);
            
            if (isempty(lub)) 
                disp('Found 0 tracks')
                tracklist={};
                return
            end
            
            % number of particle traces found
            disp(sprintf('Found %d tracks.', lub(end,end)))
            
            %% create array of individual tracks
            tracklist = {};
            [b,ind,n] = unique(lub(:,8),'first');
            for pc = 1:lub(end,8)-1
                tracklist{pc} = lub(ind(pc):ind(pc+1)-1,:);
            end
            tracklist{lub(end,8)} = lub(ind(end):end,:);
            
            %% Set up Particles array for cell
            CL.NParticle = length(tracklist);
            
            for tc = 1:length(tracklist)
                % track data for each point
                track = tracklist{tc};
                CL.Particles(tc).xy = track(:,1:2);        
                CL.Particles(tc).mass = track(:,3);
                CL.Particles(tc).rg2 = track(:,4);
                CL.Particles(tc).ecc = track(:,5);
                CL.Particles(tc).fi = track(:,6);    
    
                % append to the pts field of each frame (absolute coords)
                for pc = 1:size(track,1)
                    fc = track(pc,6);
                    CL.Frames(fc).pts.np = CL.Frames(fc).pts.np +1;
                    nn=CL.Frames(fc).pts.np;
                    trackabs = track(pc,1:2)+CL.Frames(fc).pos;
                    CL.Frames(fc).pts.xy(nn,:) = trackabs;
                    % index inside Particles list
                    CL.Frames(fc).pts.pi(nn) = tc;
                    % index within the particle track
                    CL.Frames(fc).pts.pic(nn) = pc;
                    CL.Particles(tc).fic(pc) = nn;
                end
            end
            
            CL.TrackOpt = opt;
        end
        
        function [belowthresh,frackeep,opt] = thresholdParticleTracks(CL,options)
            % belowthresh contains trajectories where velocities always
            % stay below the threshold
            % frackeep contains the fraction of total trajectory length
            % kept
            
            %% default options
            opt = struct();            
            % break trajectories at stage jump frames?
            opt.breakjump = 1;
            % wavelet span to use in getting smoothed velocity for
            % thresholding
            opt.threshnsmooth = 20;
            % threshold trajectories with smoothed velocities below so many standard
            % deviations
            opt.threshsig = 2;
            % minimal thresholded track length to keep
            opt.mintracklen=40;
            % trajectories relative to cell position
            opt.relCellPos = 0;
            % trajectories relative to contour COM
            opt.relCOM = 0;
            
            if (nargin>1)
                % input options
                opt = copyStruct(options,opt);                
            end  
            
            % ------------
            if opt.breakjump
                % break trajectories at stage jump frames
                pos = vertcat(CL.Frames.pos);
                dpos = diff(pos); ndpos = dpos(:,1).^2+dpos(:,2).^2;
                jumpind = find(ndpos>0)';
                badframes = jumpind;
            else
                badframes = [];
            end
            
            % get trajectories split at badframes
            tracklist = CL.getTracklist('mintracklen',2*opt.threshnsmooth+2,'relCellPos',opt.relCellPos,'splitskip',1,'badframes',jumpind,'relCOM',opt.relCOM);
            if numel(tracklist)>0
                % find the threshold for velocities to be above
                % threshsig*standard deviation
                [thresh,vv,smvels,meanvels,totmeanvel] = getThreshold(tracklist,opt.threshnsmooth,'flowweight',opt.threshsig,'noiseweight',0);
                
                % break tracks whenever velocities go above threshold
                [belowthresh,abovethresh] = thresholdTracks(tracklist,smvels,opt.threshnsmooth,thresh,'meanshift',vv,'minabovelen',1);
                
                lenorig = length(vertcat(tracklist{:}));
                lenthresh = length(vertcat(belowthresh{:}));
                frackeep = lenthresh/lenorig;
                
                % keep trajectories above minimal length
                tracklen = cellfun(@(x) size(x,1),belowthresh);
                goodtrackind = find(tracklen>opt.mintracklen);
                belowthresh = belowthresh(goodtrackind);
            else
                warning("No items in tracklist!");
                belowthresh = [];
                frackeep = 0;
            end
        end
        
        function opt = setCellFrameRef(CL,options)
            %% default options
            opt = struct();
            % factor by which to expand 2nd image rectangle before aligning
            opt.exprect2=1.2;
            % which frames to run
            opt.framelist = 1:CL.NFrame;
            % clear all values of cell position before running
            opt.clearcellpos= 0;
            
            % show aligned images
            opt.dodisplay=1;
            
            % do (hackish) subpixel alignment? To what resolution?
            opt.subpixel = 1;
            opt.subpixelscl = 0.01;
            
            if (nargin>1)
                % input options
                opt = copyStruct(options,opt);                
            end  
            
            if (~isfield(CL.Frames(1),'cellpos'))
                CL.Frames(1).cellpos = [];
            end
            
            %% clear
            if (opt.clearcellpos)
                CL.Frames = rmfield(CL.Frames,'cellpos');
                CL.Frames(1).cellpos = [];
            end
            
            
            %% set the cell frame of reference by overall cross-correlation
            
            fcprev = opt.framelist(1);
            if (~isfield(CL.Frames(fcprev),'cellpos'))
                CL.Frames(fcprev).cellpos = [0,0];
            end
            if (isempty(CL.Frames(fcprev).cellpos))
                CL.Frames(fcprev).cellpos = [0,0];
            end
            
            imgprev = imread(CL.Frames(fcprev).fname,CL.Frames(fcprev).fpage);
                        
            for fcc=2:length(opt.framelist)
                fcprev = opt.framelist(fcc-1);
                fc = opt.framelist(fcc);
                
                rect1 = round(CL.Frames(fcprev).rect);
                
                img = imread(CL.Frames(fc).fname,CL.Frames(fc).fpage);                                
                rect2 = CL.Frames(fc).rect;
                rect2 = expandRect(rect2,opt.exprect2,[1,CL.ImgSize(2)],[1,CL.ImgSize(1)]);
                rect2 = round(rect2);
                
                if (rect2(3)<rect1(3) || rect2(4)<rect1(4))
                    % if rect2 is too small, use whole image
                    rect2 = [1 1 CL.ImgSize(2) CL.ImgSize(1)]
                end
                
                % subtract out background                
                [offset,ccval,max_c,mean_c] = imgAlign(imgprev,img,'rect1',rect1,'rect2',rect2,'dodisplay',opt.dodisplay,'subpixel',opt.subpixel,'subpixelscl',opt.subpixelscl);
                
                if (opt.dodisplay)
                    hold all
                    text(100,100,sprintf('Frame %d: %f %f',fc,offset),'Color','y','FontSize',16)
                    hold off
                    drawnow
                end
                
                CL.Frames(fc).cellpos = CL.Frames(fcprev).cellpos + offset;
                imgprev = img;
            end
        end
        
        function interpCellFrameRef(CL,refframe)
            % interpolate cell frame for reference (cellpos) field to all
            % frames. Set such that first frame matches to COM
            % interpolate relative to stage position
            % ref frame is the reference frame where cell position = COM
            if (nargin<2)
                refframe = 1;
            end
            
            framelist = find(arrayfun(@(x) ~isempty(x.cellpos),CL.Frames));
            % cell position relative to stage position
            cellpts =  vertcat(CL.Frames(framelist).cellpos) - vertcat(CL.Frames(framelist).pos);
            COM = fnval(CL.Frames(refframe).time,CL.SpCOM);
            cellpts = bsxfun(@plus, cellpts,COM'-cellpts(1,:));           
            
            cellpostimes = [CL.Frames(framelist).time];
            cellposSp = spline(cellpostimes,cellpts');
            times = [CL.Frames.time];
            allcellpts = fnval(cellposSp,times)';
            valframes = find(arrayfun(@(x) ~isempty(x.pos),CL.Frames))
            for fc = valframes
                CL.Frames(fc).cellpos = allcellpts(fc,:)+CL.Frames(fc).pos;
            end
        end
        
        function CL2 = copyCell(CL,varargin)
            % copy over all fields to a new cell
            % optionally, specify 'Name' argument to give the new cell a
            % different name
            
            Name = CL.Name;
            for vc = 1:2:length(varargin)
                switch varargin{vc}
                    case('Name')
                        Name = varargin{vc+1};
                end
            end
            
            CL2 = CellObj(Name);
            CL2 = copyStruct(CL,CL2,'addnew',1,'exclude',{'Name','SpCont','SpCOM'});
            
            if (~isempty(CL.SpCont))
                CL2.setSpCont();
            end
        end
        
        function s2ncell = getSignal2Noise(CL,options)
            % get average signal to noise for the cell
            
            %% default options
            opt = struct();
            % which frames to average over
            opt.frames = 1:100:CL.NFrame;
            
            % show images
            opt.dodisplay=0;
            
            % width for particle squares
            opt.w = 8;
            
            % clear signal to noise data for all frames
            opt.clears2n = 0;
                        
            if (nargin>1)
                % input options
                opt = copyStruct(options,opt);
            end
                        
            %% clear s2n data
            if (opt.clears2n)
                for fc = 1:CL.NFrame
                    CL.Frames(fc).s2n = 0;
                end
            end
            
            %% get signal to noise for each frame
            s2nf = zeros(1,length(opt.frames));
            
            for fcc = 1:length(opt.frames)
                fc = opt.frames(fcc);
                %%
                frame = CL.Frames(fc);
                
                if (isempty(frame.pts.xy))
                    continue
                end
                img = imread(frame.fname,frame.fpage);
                imgnopt = img;
                
                segcont = fnval(CL.SpCont,frame.time);
                % absolute position
                segcont = bsxfun(@plus, segcont,frame.pos);
                %% avg intensity per pixel in square of size w around points
                intval = 0;
                ct = 0;
                %allintval = [];
                for i = -opt.w/2:opt.w/2
                    for j = -opt.w/2:opt.w/2
                        %%
                        x = round(frame.pts.xy(:,1)+i);
                        y = round(frame.pts.xy(:,2)+j);
                        ok = find(x<=size(img,2) & y <=size(img,1) & x>0 & y>0);
                        x = x(ok); y = y(ok);
                        ind = sub2ind(size(img),y,x);
                        %allintval = [allintval;img(ind)];
                        intval = intval + sum(img(ind));
                        ct = ct + length(ind);
                        % image with points cut out
                        imgnopt(ind) = 0;
                    end
                end
                intval = intval/ct;
                %intval = median(double(allintval))
                if (opt.dodisplay)
                    imshow(imgnopt,[])
                end
                %% avg intensity per pixel within image with the points removed
                % keeping only within the segmentation contour
                BW = poly2mask(segcont(:,1),segcont(:,2),size(img,1),size(img,2));
                imgnopt = double(imgnopt).*BW;
                
                if (opt.dodisplay)
                    imshow(imgnopt,[])
                end
                
                % background intensity
                bgval = sum(imgnopt(:))/sum(BW(:));
                
                s2nf(fcc) = intval/bgval;
                CL.Frames(fc).s2n = s2nf(fcc);
                %[fc s2nf(fc)]
            end
            
            % average over all the frames
            s2ncell = mean(s2nf(find(s2nf>0)));
            
            CL.sig2noise = s2ncell;
        end
        
        function getAvgBrightness(CL,varargin)
            % calculate the average brightness in the frame at large and in
            % the cell
            % save to avgbr and avgbrseg fields of the Frames structure

            % for each frame, calculate overall median brightness and median brightness within cell
            for fc = 1:CL.NFrame
                if (mod(fc,20)==0); disp([fc CL.NFrame]); end
                frame = CL.Frames(fc);
                img = imread(frame.fname,frame.fpage);
                %rect1 = frame.rect;
                %img = imcrop(img,rect1);
                
                segcont = fnval(CL.SpCont,frame.time);
                %segcont = bsxfun(@plus,segcont, frame.pos-rect1(1:2)+1);
                segcont = bsxfun(@plus,segcont, frame.pos+1);
                [CL.Frames(fc).meanbr, CL.Frames(fc).meanbrseg] = getImgBrightness(img,segcont,varargin{:});
            end
    
        end
        
        function [outliers,goodind] = getPhaseFrames(CL,varargin)
            % get phase frames as brightness outliers
            
            pct = 50; 
            multiplier = 0.2;
            dodisplay = 0;
            
            for vc = 1:2:length(varargin)
                switch(varargin{vc})
                    case('pct')
                        pct = varargin{vc+1};
                    case('multiplier')
                        multiplier = varargin{vc+1};
                    case('dodisplay')
                        dodisplay=varargin{vc+1};
                end
            end
            %% robust smoothing
            framelist = 1:CL.NFrame;
            avgbr = [CL.Frames.avgbr];
            smbr = smooth(framelist,avgbr,100,'rlowess');
            %plot(segframelist,avgbr,'.-',segframelist,smbr)
            
            %% cutoff on residuals to find phase frames
            res = avgbr-smbr';
            %hist(res)
            outliers = find(res>prctile(abs(res),pct)*multiplier);
            goodind = find(res<prctile(abs(res),pct)*multiplier);
            if (dodisplay)
                plot(framelist,avgbr,'.-',framelist(goodind),avgbr(goodind),'.-',framelist,smbr)
            end
        end
         
        function [M,opt] = viewTrackMovie(CL,savefile,options)
            % view movie of tracks superimposed on original images
            
            % default options
            opt = struct();
            
            % additional options to pass to getTrackList;
            opt.tracklistopts = {};
            
            % which frames to show
            % default is all the ones where tracks exist
            opt.frames = [];
            
            % additional buffer for overall movie (on each side)
            opt.buffer = 100;
            
            % brightness range (absolute)
            opt.imrange = [];
            
            % how much to expand cropping rectangle beyond cell outline
            opt.rectexpand = 1.1;
            
            % print a time label on the movie
            opt.printlabel = 1;
            
            % minimal frame at which to start tracks            
            opt.minframe = 1;
            
            % shift time label
            opt.timeshift = 0;
            
            % how often to save
            opt.saveevery = 100;
            
            % add track number
            opt.labeltracks = 0;
            
            % skip certain track numbers
            opt.badtracks = [];
            
            % specific tracklist to use
            % default, pull via getTracklist
            opt.tracklist = [];
            
            % still show tracks that finished before current frame            
            opt.showfinished = 1;
            
            % color tracks based on value in certain column
            % default: use track number instead
            opt.colorcol = []
            
            % colormap
            opt.cmap = [];
            
            % plot point at end of track
            opt.endpt = 1;
            
            % crop image around cell
            opt.docrop = 1;
            % shift to account for stage shift
            opt.shiftpos = 1;
            
            % copy over input options
            if (nargin>2)
                opt = copyStruct(options, opt);
            end

            
            %% get dimensions of overall movie
            
            if (isempty(opt.tracklist))
                tracklist = CL.getTracklist(opt.tracklistopts{:});
            else
                tracklist = opt.tracklist;
            end
            
            allpts = vertcat(tracklist{:});
            minfr = min(allpts(:,6));
            maxfr = max(allpts(:,6));
            minx = min(allpts(:,1));
            maxx = max(allpts(:,1));
            miny = min(allpts(:,2));
            maxy = max(allpts(:,2));
                       
            minx = minx-opt.buffer; 
            miny = miny-opt.buffer; 
            maxx = maxx+opt.buffer; 
            maxy = maxy+opt.buffer;            
            
            %% run the movie
            if isempty(opt.frames)
                frameind = minfr:maxfr;
            else
                frameind = opt.frames;
            end

            if (isempty(opt.cmap))  
                cmap = lines(length(tracklist));
            else
                cmap = opt.cmap;
            end
            
            % starting and ending frame of each track
            firstframes = cellfun(@(x) x(1,6),tracklist);
            lastframes = cellfun(@(x) x(end,6),tracklist);
            
            % original image size
            img = imread(CL.Frames(frameind(1)).fname, CL.Frames(frameind(1)).fpage);
            cutx = [1,size(img,2)];
            cuty = [1,size(img,1)];
            
            ct = 0;
            for cc = frameind
                img = imread(CL.Frames(cc).fname, CL.Frames(cc).fpage);
    
                % all tracks that have shown up by this frame
                if (opt.showfinished)
                    trackind = find(firstframes<=cc);      
                else
                    trackind = find(firstframes<=cc & lastframes>=cc);
                end
                
                % crop around cell outline
                if (opt.docrop)
                    segcont = fnval(CL.SpCont,CL.Frames(cc).time);
                    croprect = [min(segcont(:,1)),min(segcont(:,2)),range(segcont(:,1)),range(segcont(:,2))];
                    croprect(1:2) = croprect(1:2) + CL.Frames(cc).pos;
                    croprect = expandRect(croprect,opt.rectexpand,cutx,cuty);
                    
                    % crop around cell
                    imgclean = imcrop(img,croprect); 
                else
                    imgclean = img;
                    minx = 1; miny = 1;
                    maxx = size(img,2); maxy = size(img,1);
                    croprect = [1 1 size(img,1) size(img,2)]
                end
              
                if (opt.shiftpos)
                    % translate image, rounding to nearest pixel
                    imgshift = zeros(round(maxy-miny),round(maxx-minx));
                    vec = round(-CL.Frames(cc).pos+(croprect(1:2)-[minx miny]));
                    imgshift(vec(2)+1:vec(2)+size(imgclean,1),vec(1)+1:vec(1)+size(imgclean,2)) = imgclean;
                else
                    imgshift=imgclean;
                end
                
                % display image
                imgscl = double(imgshift);
                imgscl = imgscl/max(imgscl(:));
                imshow(imgscl,opt.imrange)
       
                % plot the tracks
                hold all
                for tc = trackind
                    track = tracklist{tc};
                    
                    if (ismember(track(end,8),opt.badtracks)); continue; end
                    
                    goodind = find(track(:,6)>=opt.minframe);
                    if (isempty(goodind)); continue; end
                    track = track(goodind,:);
                    
                    % find endpoint for track
                    if (track(end,6)<cc)
                        endind = size(track,1);
                    else
                        endind = find(track(:,6)>=cc); endind = endind(1);
                    end
%                     if (track(endind,6)>cc)
%                         endind = endind-1;
%                     end
                    %%
                    if (opt.shiftpos)
                        trackshift = bsxfun(@minus,track(1:endind,1:2),[minx-1,miny-1]);
                    else
                        trackshift = bsxfun(@minus,track(1:endind,1:2),-CL.Frames(cc).pos+1);
                    end
                    if (isempty(opt.colorcol))
                        plot(trackshift(:,1),trackshift(:,2),'-','Color',cmap(tc,:));
                        if (opt.endpt); plot(trackshift(end,1),trackshift(end,2),'.','Color',cmap(tc,:)); end
                    else
                        plot(trackshift(:,1),trackshift(:,2),'-','Color',cmap(track(1,opt.colorcol),:));
                        if (opt.endpt); plot(trackshift(end,1),trackshift(end,2),'.','Color',cmap(track(1,opt.colorcol),:)); end
                    end
                    if (opt.labeltracks)
                        text(trackshift(end,1),trackshift(end,2),sprintf('%d',track(end,8)),'Color',cmap(tc,:))
                    end
                end
    
                if (opt.printlabel)
                   lbl = sprintf('%d sec', round(CL.Frames(cc).time+opt.timeshift))
                   text(size(imgshift,2)-200,size(imgshift,1)-60,lbl,'Color','y','FontSize',50)
                end
                
                hold off
                set(gca,'Visible','off','Position',[0.05 0.05 0.9 0.9],'XTick',[],'YTick',[])
                drawnow
                
                ct = ct+1;
                M(ct) = getframe;
                
                if (mod(ct,opt.saveevery)==0)
                    save(savefile,'M')
                end
            end
            
             save(savefile,'M')
        end
        
       function makeShortMovie(CL,tiffile,varargin)
            %% make short, downsampled movie of cell motion for easy viewing

            frameskip = 10; % sample every so many frames
            expandcont = 50; % how far out to expand beyond cell contour
            resizescl = 0.3; % scaling factor for resizing image
            dodisplay = 0; % display frames as they're drawn
            
            for vc = 1:2:length(varargin)
                switch varargin{vc}
                    case('frameskip')
                        frameskip = varargin{vc+1};
                    case('expandcont')
                        expandcont = varargin{vc+1};
                    case('resizescl')
                        resizescl = varargin{vc+1};
                    case('dodisplay')
                        dodisplay = varargin{vc+1};
                end
            end
            
            filename = tiffile;    
    
            framelist = frameskip:frameskip:CL.NFrame;
            for fcc = 1:length(framelist)
                fc=framelist(fcc);
                
                img = imread(CL.Frames(fc).fname, CL.Frames(fc).fpage);                
                
                %% expanded mask around cell
                segcontrel = fnval(CL.SpCont,CL.Frames(fc).time);
                segcont = bsxfun(@plus, segcontrel, CL.Frames(fc).pos);
                [x,y] = poly2ccw(segcont(:,1),segcont(:,2));
                contccw = [x y];
                contexp = expandPolygon(contccw,expandcont);
                clen = cellfun(@(x) length(x),contexp);
                [mlen,maxind] = max(clen);
                contexp = contexp{maxind};
                contexp = [contexp;contexp(1,:)];
                mask = poly2mask(contexp(:,1),contexp(:,2),size(img,1),size(img,2));
                
                img = img.*uint16(mask);
                if (resizescl<1)
                    img = imresize(img,resizescl);
                end
                
                if (dodisplay)
                    imshow(img,[])
                    drawnow
                end
                %%
                if (fcc==1)
                    imwrite(img, filename,'Compression','none');
                else
                    imwrite(img, filename, 'WriteMode', 'append','Compression','none');
                end
            end
        end
    end
end


