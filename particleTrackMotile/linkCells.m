function [cellind, COMlist, cellmv] = linkCells(cellmv,varargin)
% starting from a cell movie with segmented contours, link up the contours
% to keep track of individual cells
% returns array cellind where cellind(i,j,n) = index of i-th cell in j-th
% frame (index within segcont array)
% NaN if cell does not exist in that frame
% COMlist is a cell list of center of mass tracks for individual cells
% (relative to cell pos field)
% cellmv is an updated cellmv structure with a COM field tracking the COMs
% (in absolute coordinates)
% indexing within this array is same as for segcont field
% also contains field cellid giving the cell to which each contour is
% assigned (NaN if not assigned)

% maximum displacement in pixels
% negative means use sqrt(avg. area);
maxdisp = -1;
% allow memory over so many frames
% SOME KIND OF BUG HERE. WHY DOES THIS FAIL FOR INTEGER VALUES?
memory = 2;
%memory = 2*mean(diff(COMpts(:,6)));
% tracks must be at least this long
goodenough = 2;

for vc = 1:2:length(varargin)
    switch(varargin{vc})
        case('maxdisp')
            maxdisp = varargin{vc+1};
        case('memory')
            memory = varargin{vc+1};
        case('goodenough')
            goodenough = varargin{vc+1};
    end
end


%% valid frames are those with both segcont and pos defined
valframes = find(arrayfun(@(x) ~isempty(x.segcont) & ~isempty(x.pos), cellmv));
%arrayfun(@(x) length(x.segcont),cellmv(valframes))
%% track COM for each cell in each frame
COMpts = []; % keep track of all COMs
areas = []; % keep track of cell size
for fc = valframes
    ncell = length(cellmv(fc).segcont); % number of segmented cells in this frame
    cellmv(fc).COM = zeros(ncell,2);
    for cc = 1:ncell
        tmp = polygeom(cellmv(fc).segcont{cc}(:,1),cellmv(fc).segcont{cc}(:,2));
        cellmv(fc).COM(cc,:) = tmp(2:3);
        COMpts = [COMpts; [tmp(2:3)-cellmv(fc).pos,zeros(1,3),fc]];
        areas = [areas;tmp(1)];
    end        
end

%% get cell tracks by linking COMs in the same way as particles (from Kilfoil code)

%% get particle tracks by linking particles
dim = 2;
if (maxdisp<0)
    maxdisp = round(sqrt(mean(areas)));
end


% pad with last column of zeros
npt = length(COMpts);
xyzs = COMpts(:,[1:2,6]);


% older stuff trying to deal with kilfoil code
%xyzs = [COMpts,zeros(npt,1)];
% add phantom particles in each frame
%nf = length(valframes);
%xyzs = [xyzs;[ones(nf,2)+rand(nf,2),zeros(nf,3),valframes',zeros(nf,1)]];
% sort by time
%xyzs = sortrows(xyzs,6);
% adjust time vector
%xyzs(:,7) = xyzs(:,6); % save actual frame number for sorting later
%mdt = max(1,round(mean(diff(COMpts(:,6)))));
%xyzs(:,6) = round((xyzs(:,6)-min(xyzs(:,6)))/mdt)+1;


% shift positions to avoid negative coordinates since this breaks the
% tracking code
xshift=min(xyzs(:,1))-1;
yshift = min(xyzs(:,2))-1;
xyzs(:,1) = xyzs(:,1)-xshift;
xyzs(:,2) = xyzs(:,2)-yshift;
%%
%[lub] = trackmem(xyzs,maxdisp,dim,goodenough,memory);
% use dufresne code
%maxdisp=800;
lub = track_dufresne(xyzs,maxdisp,struct('mem',memory,'good',goodenough,'dim',2,'quiet',0));
lub(:,1) = lub(:,1)+xshift;
lub(:,2) = lub(:,2)+yshift;

% HACK: add last point to previous track (some bug in dufresne code)
% this should be okk since real tracks that are only one timepoint long
% will be thrown out within the tracking routine;
if lub(end,4) ~= lub(end-1,4)
    lub(end,4) = lub(end-1,4);
end
%% convert back to format where 6th column has frame number; 8th column has track number
lub = [lub(:,1:2),zeros(size(lub,1),3),lub(:,3),zeros(size(lub,1),1),lub(:,4)];

% get back actual frame numbers
%lub(:,6) = lub(:,7);
%lub(:,7) = 0;
%lub = lub((lub(:,1)>2 | lub(:,2)>2),:);

% number of particle traces found
disp(sprintf('Found %d tracks.', lub(end,end)))

%% create array of individual tracks
COMlist = {};
[b,ind,n] = unique(lub(:,8),'first');
for pc = 1:lub(end,8)-1
    inds = ind(pc):ind(pc+1)-1;
    COMlist{pc} = lub(inds,:);
end
inds = ind(end):size(lub,1);
COMlist{lub(end,8)} = lub(inds,:);

%% get indices of each cell for the segmentation contours / COMs in each frame
cellind = NaN*ones(length(COMlist),length(cellmv));

for fc= 1:length(cellmv)
    cellmv(fc).cellid = NaN*ones(1,length(cellmv(fc).segcont));
end

for cc = 1:length(COMlist)
    for tc = 1:size(COMlist{cc},1)
        fc = COMlist{cc}(tc,6);
        % find nearest COM in this frame
        diffs = bsxfun(@minus, cellmv(fc).COM,cellmv(fc).pos+COMlist{cc}(tc,1:2));
        [~,b] = min(abs(sum(diffs.^2,2)));
        cellind(cc,fc) = b;
        cellmv(fc).cellid(b) = cc;
    end
end