% set up paths
restoredefaultpath
addpath('./cellobj')
addpath(genpath('./tools'))
addpath(genpath('C:\Users\pseudopod\Documents\Code\bfmatlab'))
% path to the BNEW code
% CHANGE THIS to the location of your BNEW repository
BNEWRelDir = '../BNEW-master/BNEW-master/';
addpath(BNEWRelDir)


%% Some general parameter values

% path to the directory containing the tif image stacks
storagedirpath = 'F:\DATA 01\AndrewK\exampleMovies\';
% glob for the name of the directory containing each movie
% %d will be replaced with the movie number
storagedirname = 'mitored_example';

% directory name in which to save results arising from processing the
% movies
dirname = './example/results/';
% directory name in which to save individual cell objects
celldirname = './example/results/';

% prefix used for naming files associated with each movie
% integer movie number is inserted for %d
prefixglob = 'Ng_example%d';

% which movie numbers will we be processing?
% can process multiple ones automatically if multiple movies taken
% during the same imaging session
movnums = 1; 

% segment cell outline every so many frames
segevery = 5;

% how much displaying to do during segmentation, particle tracking, etc.
dodisplay = 1;

% threshold out trajectories of unusually fast-moving particles
dothresh = 1;
% break trajectories when velocity goes above the given number of stdevs
threshsig = 2;
% wavelet span to use in calculating smoothed velocities for thresholding
threshnsmooth = 20;

% minimal track length to keep for doing BNEW analysis
mintracklen = 40;

% wavelet spans to use for BNEW analysis
% recommended for HL60s: 2-17 if doing trajectory thresholding, 2-13
% otherwise
bnewnvals = 2:17;

% file containing precalculated BNEW coefficients
% if file not available or using something other than svg3 wavelets, set to
% empty
bnewcoeffile = fullfile(BNEWRelDir,'svg3waveletcoeff.mat');

% --------------------
%% The cell segmentation is intervention-heavy
% It is recommended that the parameters are adjusted manually by trying to
% segment the first and last frame and making sure it looks good

% some default parameter values 
threshscl = 0.5*ones(1,length(movnums));
startframe = ones(1,length(movnums));
adjmin = 0.01*ones(1,length(movnums));
adjmax = 0.99*ones(1,length(movnums));

%% Use GUI to adjust initial parameters for segmentation 
% resulting segmentation options will be saved into file example/results/example15_seginit.mat
ind = 1;
% output file for initial segmentation calculations
outputfile = [dirname prefixglob '_seginit.mat']
tuneSegmentationParams('outputfile',outputfile,'mnum',movnums(1),'threshscl',threshscl(ind),...
    'startframe',startframe(ind),'adjmin',adjmin(ind),'adjmax',adjmax(ind),...
    'storagedirpath',storagedirpath,'storagedirname',storagedirname)
% will load segmentation parameters from file
loadsegparams = ones(1,length(movnums));

%% ---------------MAIN LOOP FOR MULTIPLE MOVIES ----------
% once the segmentation parameters have been adjusted and saved, 
% should be able to let this run on autopilot

% Process each of the desired imaging movies
for mc = 1:length(movnums);
    %%
    % mnum is the current movie number
    mnum = movnums(mc);
    
    % prefix that will be used to name each of the cells in this movie
    prefix = sprintf(prefixglob,mnum)
    
    if (loadsegparams(mc))
        % load movie data and segmentation parameters from file
        load([dirname prefix '_seginit.mat'],'cellmv')
    else
        % extract movie data
        
        % glob (directory name, with * wildcards) for the directory
        % containing tif stacks for this specific movie
        dirglob = fullfile(storagedirpath, sprintf(storagedirname,mnum))
        
        % get list of tif stacks for this movie
        storagedirname = dir(dirglob);
        storagedirname = fullfile(storagedirname(1).folder, filesep)
        files = dir(fullfile(storagedirname, '*.ome.tif'));
        files = {files.name};
        
        % get metadata file
        metafile = dir([storagedirname '*metadata.xml']);
        metafile = metafile.name;
        disp(sprintf('Processing data for Cell %d: %d',mc,mnum))
        
        % read in movie data
        disp('Reading in data from file...')
        cellmv = loadMovieData(storagedirname,files,metafile);
    end
    %% set frames for cell segmentation
    
    % set this to first and last for adjusting parameters
    % segframelist = [startframe(mc), length(cellmv)];
    
    % otherwise segment every segevery frames, including 1st and last
    if (mod(startframe(mc),segevery)==0)
        segframelist=[startframe(mc):segevery:length(cellmv),length(cellmv)];
    else
        tmp = ceil(startframe(mc)/segevery)*segevery;
        segframelist = [startframe(mc), tmp:segevery:length(cellmv),length(cellmv)];
    end
    
    % file for saving segmentation results
    segsavefile = [dirname prefix '_segcont.mat'];
    
    %% Run the approximate cell segmentation
    % brjumpcutoff: set to negative means don't worry about intermittent
    % phase frames
    % clearsegcont: 0 means don't redo segmentation of frames that have
    % already been segmented
    if (loadsegparams(mc))
        % load segmentation parameters from file
        load([dirname prefix '_seginit.mat'],'segoptions')
        segoptions.dodisplay = dodisplay;
        segoptions.clearsegcont = 0;
        segoptions.brjumpcutoff = -1;
        segoptions.segframelist = segframelist;
    else
        segoptions = struct('dodisplay',dodisplay,'adjmin',adjmin(mc),'adjmax',adjmax(mc),'threshscl',threshscl(mc),...
            'clearsegcont',0,'brjumpcutoff',-1,'segframelist',segframelist);
    end
    [cellmv, segopt] = runSegmentation(cellmv,segsavefile,segoptions);
    
    %% align movie frames based on segmented cells (deals with stage jumps)
    % use largest cell contour to specify area for stage jump alignment
    alignoptions = struct('dodisplay',dodisplay,'getsegcont',0);
    [cellmv,jumps,alignopt] = cellAlignDriver(cellmv,alignoptions);        
    
    %% save current results
    savefile = [dirname prefix '_results.mat'];
    save(savefile, 'cellmv','segopt','alignopt','segsavefile','savefile',...
        'dirname','prefix','storagedirpath')
    
    %% convert to individual cell objects
    disp(sprintf('Converting data to cell files...'))
    [Cells,cellmv] = cellObjFromData(savefile,'outputdir','','loadparticles',0,'maxdisp',80,'memory',2);
    disp(sprintf('%d Cells Found', length(Cells)))
    
    % filenames in which to save data for individual cells
    for ccell = 1:length(Cells)
        cellfiles{ccell} = [celldirname sprintf('%s.mat',Cells(ccell).Name)];
    end
    
    % view each of the cell objects
    if (dodisplay)
        for sc = 1:length(Cells)
            subplot(1,length(Cells),sc)
            Cells(sc).viewFrame(1,'docrop',0);
            xlabel(sprintf('Cell %d',sc))
        end
    end
    
    %% find particles for each cell    
    for ccell = 1:length(Cells)
        %%
        CL = Cells(ccell);
        partoptions = struct('expandcont',50,'dodisplay',dodisplay,'savefile',cellfiles{ccell});
        partopt = CL.findParticles(partoptions);        
        save(cellfiles{ccell},'CL','partopt')        
    end        
    
    %% link up particles into tracks for each cell
    for ccell = 1:length(Cells)        
        CL = Cells(ccell);
        [tracklist,trackopt] = CL.linkTracks();
        save(cellfiles{ccell},'CL','trackopt')        
    end
    
    %% Run BNEW analysis    
    for ccell = 1:length(Cells)
        CL = Cells(ccell);
        
        if (dothresh)
            % Threshold trajectories to get rid of unusually fast-moving particles
            threshoptions = struct('mintracklen',mintracklen,'threshsig',threshsig,'threshnsmooth',threshnsmooth);
            [tracklist, frackeep] = CL.thresholdParticleTracks(threshoptions);
        else       
            % break trajectories at stage jump frames
            pos = vertcat(CL.Frames.pos);
            dpos = diff(pos); ndpos = dpos(:,1).^2+dpos(:,2).^2;
            jumpind = find(ndpos>0)';                           
            tracklist = CL.getTracklist('mintracklen',mintracklen,'splitskip',1,'badframes',jumpind);            
            frackeep = 1;
        end
        % Now do the BNEW analysis
        BN = BNEWobj(bnewnvals);
        BN = BN.getCoefficients();
        if numel(tracklist) > 0
            BN = BN.analyzeTracks(tracklist);       
            if (isempty(bnewcoeffile))
                BN = BN.rescaleData('recalc',1);
            else
                BN = BN.loadWaveletCoeff(bnewcoeffile);
                BN = BN.rescaleData('recalc',0);
            end        
            % fit diffusion coefficient and scaling
            % average timestep
            dt = mean(diff([CL.Frames.time]));
            fitopt = struct('kmax',0.74,'fitd0',[3,0.1,1.2],'fitfunc',@(D,t) (4*D(1)*(dt*t).^D(3)+4*D(2)^2));
            [BN,Dfit,~,allMSE(ccell)] = BN.fitDcoeff(fitopt);
            %
            BN.plotRescaled('kmax',0.74,'doplot',@plot);
            
            % BN.Dfit gives the following, in order
            % 1) diffusion coefficient D (units of px^2/sec)
            % 2) localization error (units of px)
            % 3) scaling coefficient alpha
            % frackeep is the fraction of total trajectory length that survived
            % the thresholding away of unusually fast trajectories
            
            [ccell BN.Dfit frackeep]
        else
            warning("No tracks in this cell!")
        end
        CL.Wavelet = BN;
        save(cellfiles{ccell},'CL')        
    end
    
end

% ----------------
%% reload cell objects in case you need to restart part-way
% load all cell files with this prefix
prefix = 'example15'
% directory to load from
dirname = './example/results/';
filenames = {};
for ccell = 1:3
    filenames{ccell} = sprintf('%s_%d.mat',prefix,ccell);
    cellfiles{ccell} = [dirname filenames{ccell}];
end
Cells = loadCellObj(dirname,filenames);
% ----------------