function [Cells,cellmv,tracklist] = cellObjFromData(matfile,varargin)
% create cell objects based on precomputed data file

rectexpf = 2; % expansion factor for defining cropping rectangle around each cell
% output directory to save cell objects
% if empty, do not save
outputdir = '';
% also save cells that are empty of particles
saveemptycells = 1;
% load information on particle tracks and sort into appropriate cells
loadparticles=1;
memory=20;
% minimal number of frames to keep a cell
goodenough = 2;

%%
for vc = 1:2:length(varargin)
    switch (varargin{vc})
        case('rectexpf')
            rectexpf = varargin{vc+1};
        case('outputdir')
            outputdir=varargin{vc+1};
        case('loadparticles')
            loadparticles = varargin{vc+1};
        case('saveemptycells')
            saveemptycells = varargin{vc+1};
        case('goodenough')
            goodenough = varargin{vc+1};   
        case('memory')
            memory= varargin{vc+1};
    end
end

%%
%matfile = '~/data/Caleb/lysophase20141113/lysoctrl2_segcont.mat';
%matfile = '~/data/Caleb/lysophase20140821/lysoph6/lysoph6_segcont.mat';
%matfile = '~/data/Rikki/lysotracker20141112/lysorikki2_segcont';
%% load data file
if (~exist(matfile,'file'))
    error('File %s does not exist.', matfile)
end

load(matfile)

%% get image size
img = imread(cellmv(1).fname,cellmv(1).fpage);
imsize  = size(img);
%% identify individual cells based on movie segmentation results

display('Linking up cell COMs')
% link up individual cell tracks within the movie
[cellind, COMlist, cellmv] = linkCells(cellmv,'memory',memory,varargin{:});

%% create cell objects
ncell = length(COMlist);

Cells = [];
for cc = 1:ncell
    name = sprintf('%s_%d', prefix,cc);
    CL = CellObj(name);    
    
    Cells = [Cells,CL];
end

%% populate variables for each cell (excluding particle stuff)
display('Setting up frame information for cells...')
for cc = 1:ncell
    CL = Cells(cc);
    
    CL.ImgSize = imsize;
    CL.DirName = dirname;
    CL.LoadFile = matfile;
    
    % data on each frame
    CL.NFrame = length(cellmv);
    CL.Frames = cellmv;
    CL.clearParticles;
    
    %% keep only the appropriate contour and COM    
    valframes = find(~isnan(cellind(cc,:)));
    badframes = find(isnan(cellind(cc,:)));
    
    inds = cellind(cc,:);
    
    for fc = valframes
        CL.Frames(fc).segcont = CL.Frames(fc).segcont{inds(fc)};
        CL.Frames(fc).COM = CL.Frames(fc).COM(inds(fc),:);
    end
    for fc = badframes
        CL.Frames(fc).segcont=[];
        CL.Frames(fc).COM = [];
    end   
    
    if isfield(CL.Frames,'tracks')
        CL.Frames = rmfield(CL.Frames,{'tracks'});
    end
    if (isfield(CL.Frames,'fc'))
        CL.Frames = rmfield(CL.Frames,'fc');
    end
    
    %% trim frames at start and end; throw out those outside first and last segmented frame
    minfr = min(valframes);
    maxfr = max(valframes);
    CL.Frames = CL.Frames(minfr:maxfr);
    CL.NFrame = length(CL.Frames);
    
    %% set spline of cell contours for interpolation over time
    CL.setSpCont();
    
    % set cropping rectangle for each frame
    CL.setRect(rectexpf);  
    
   % valframesave{cc} = valframes;
end
%%
if (loadparticles)
%% sort particle tracks into specific cells
disp('Sorting particles into cells...')
for tc = 1:length(tracklist)
    track = tracklist{tc};
    
    % check which cell this belongs to
    valframes = find(arrayfun(@(x) ~isempty(x.COM),cellmv(track(:,6))));
        
    if (isempty(valframes))
        disp(sprintf('Particle %d does not go over any segmented frames',tc))
        continue
    end
    
    ci = zeros(1,length(valframes));
    for ct = 1:length(valframes)
        fc = valframes(ct);        
        fn = track(fc,6);
        trackabs = track(fc,1:2)+cellmv(fn).pos; % position in absolute coords
        % check if its within a segmented contour
        for ccell = 1:length(cellmv(fn).segcont)
            segcont = cellmv(fn).segcont{ccell};
            IN = inpolygon(trackabs(1),trackabs(2),segcont(:,1),segcont(:,2));
            if (IN)
                ci(ct) = cellmv(fn).cellid(ccell);
                break
            end
        end
        if (~IN) % points not inside any segmented contour
            diffs = bsxfun(@minus, cellmv(fn).COM,trackabs);
            % find nearest cell
            [~,b] = min(abs(sum(diffs.^2,2)));
            ci(ct) = cellmv(fn).cellid(b);
            
            % check if within cell rectangle
            rect = Cells(ci(ct)).Frames(fn).rect;
            if (trackabs(1)<rect(1) || trackabs(1)>rect(1)+rect(3) || ...
                    trackabs(2)<rect(2) || trackabs(2)>rect(2)+rect(4))
                %outside of rectangle
                ci(ct) = NaN;
            end
        end
    end
    
    if (any(isnan(ci)))
        warning('Particle %d does not belong to any cell at some point in the trajectory',tc)
        continue
    end
    if (any(abs(ci-ci(1))>0))
        warning('Particle %d spans multiple cells, or partially exits cell rectangle',tc)
        continue
    end
             
    % append to cell's list of particles
    CL = Cells(ci(1));
    CL.NParticle = CL.NParticle+1;
    np = CL.NParticle;
    %CL.Particles(np) = FeatureObj(track(:,1:2));
    
    CL.Particles(np).xy = track(:,1:2);        
    CL.Particles(np).fi = track(:,6);    
    
    % append to the pts field of each frame
    for pc = 1:size(track,1)
        fc = track(pc,6);
        CL.Frames(fc).pts.np = CL.Frames(fc).pts.np +1;        
        nn=CL.Frames(fc).pts.np;
        trackabs = track(pc,1:2)+CL.Frames(fc).pos;
        CL.Frames(fc).pts.xy(nn,:) = trackabs;
        % index inside Particles list
        CL.Frames(fc).pts.pi(nn) = np;
        % index within the particle track
        CL.Frames(fc).pts.pic(nn) = pc;        
        CL.Particles(np).fic(pc) = nn;
    end
end
end

% %% trim frames from each cell object
% % cut out those outside of first:last segmented frame
% for cc = 1:length(Cells)
%     CL = Cells(cc);
%     minf = valframesave{cc}(1);
%     maxf = valframesave{cc}(end);
%         
%     CL.Frames = CL.Frames(minf:maxf);
%     CL.NFrame = length(CL.Frames);
%     
%     for pp = 1:CL.NParticle
%         CL.Particles(pp).fi = CL.Particles(pp).fi -minf + 1;  
%     end
% end


%% save data in output files

if (~isempty(outputdir))    
    for cc = 1:length(Cells)
        CL = Cells(cc);
        if (saveemptycells || CL.NParticle > 0)
            savefile = [outputdir,Cells(cc).Name,'.mat'];
            disp(sprintf('Saving cell %d to %s', cc, savefile));
            save(savefile,'CL')
        end
    end
end
end
