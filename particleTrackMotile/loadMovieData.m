function cellmv = loadMovieData(dirname,files,metafile,varargin)

% go through a movie consisting of sequential multipage tiff files
% create a struct array containing relevant data about the movie
% each element of the array is a frame
% If metafile is empty sets times to 0.05sec intervals by default
% fields in the array are:
% fc = overall frame number (starting from 1)
% fname = full file name (including path)
% fpage = index in the original file (starting from 1)
% time = time for that frame, in sec; starts at 0

% example input for debugging
%dirname = '/data/lenafabr/Caleb/lysotracker20140716/20140716 - Zyla Poly-L-lysine and LysoTracker/M1_100x_1.0x_PolyLysine_LysoTracker_CellG_2 (mCherry - G)/';
%files = {'M1_100x_1.0x_PolyLysine_LysoTracker_CellG_2_MMStack.ome.tif'};
%for c = 1:2
%    files{c+1} = sprintf('M1_100x_1.0x_PolyLysine_LysoTracker_CellG_2_MMStack_%d.ome.tif',c);
%end

%metafile = 'M1_100x_1.0x_PolyLysine_LysoTracker_CellG_2_MMStack_metadata.txt';

% load every nth frame
framestep = 1;
for vc = 1:2:length(varargin)
    switch(varargin{vc})
        case('framestep')
            framestep = varargin{vc+1};
    end
end
%%
disp(sprintf('Number of files: %d', length(files)))
%% create struct array and
% find total number of frames

cellmv = struct();

for c = 1:length(files)
    fname = [dirname files{c}];
    info = imfinfo(fname);
    
    nframe(c) = length(info);    
    
    inds = 1:framestep:nframe(c);    
    
    if (c>1)
        lastind = cellmv(end).fc;
        newmv = struct('fc',num2cell(inds+lastind),'fname',fname,'fpage',num2cell(inds));
        cellmv = [cellmv, newmv];
    else
        cellmv = struct('fc',num2cell(inds),'fname',fname,'fpage',num2cell(inds));
    end
    
    disp(sprintf('File %s has %d frames.',fname, nframe(c)))
end
nftot = sum(nframe);

disp(sprintf('Total Number of frames: %d', nftot))

%% read metadata (with Bioformats) to find time for each frame
fname = [dirname files{1}];
[time, eltime] = getTimeBetweenFrames(fname);
for c = 1:length(cellmv)
    cellmv(c).time = time(c); %in seconds
    cellmv(c).eltime = eltime(c); %in seconds
end
%%
for c = 1:nftot
    cellmv(c).bad = 0;
end
