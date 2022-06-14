classdef BNEWobj
    % class defining BNEW analysis results and methods
    properties
        % spans of wavelet
        Nvals = [];
        % Coefficients of position points (for each span n)
        Ws = {};
        % Coefficients of velocities (for each span)
        Whats = {};
        % type of wavelet (svg=Savitzky-Golay,haar,mean)
        WaveType='svg';
        % degree of velocity smoothing for svg wavelet
        WaveDeg=3;
        % number of tracks analyzed
        NTrack = {};
        % corrected mean squared displacements (for each span)
        MSD = {};
        % times (as frame numbers) associated with MSDs
        Tbin = {};
        % number of data points in each adjusted MSD calculation
        Cnt = {};
        % standard error of each adjusted MSD calculation
        Sterr = {};
        % smoothed velocities for each track, for each span
        SmVel = {}
        
        % rescaled times, MSD
        Tscaled = {};
        MSDscaled = {};
        
        % two functions used for rescaling
        % for Brownian motion with constant drift,  
        % MSD = 4*D*A + 4*epsilon^2*B;
        % listed for each n
        Afunc = {};
        Bfunc = {};
        
        % function for fitting diffusion coefficient/scaling
        FitFunc=[];
        % Diffusion coefficient, localization error, scaling
        Dfit = [];
        % covariance matrix for fitting; this is not really useful since
        % cannot assume independent errors in data points
        CovFit = [];
        
        % fit was done for log-log data
        FitLog = 0;
    end
    
    methods
        function BN = BNEWobj(nvals,varargin)
            %Initialize the BNEW object with a particular set of spans
            BN.Nvals = nvals;
            
            % default arguments
            BN.WaveType='svg';
            BN.WaveDeg=3;
            
            for vc = 1:2:length(varargin)
                switch (varargin{vc})

                    case('wavetype')
                        BN.WaveType=varargin{vc+1};
                    case('wavedeg')
                        BN.WaveDeg=varargin{vc+1};
                    case('Afunc')
                        BN.Afunc = varargin{vc+1};
                    case('Bfunc')
                        BN.Bfunc = varargin{vc+1};
                end
            end
        end
        
        function BN = getCoefficients(BN)
            % get the coefficients defining the wavelet
            
            for cc = 1:length(BN.Nvals)
                if (strcmp(BN.WaveType,'svg')) % Savitzky-Golay wavelet
                    [ws, what] = svgwavelet(BN.Nvals(cc),BN.WaveDeg);
                elseif (strcmp(BN.WaveType,'haar')) % Haar wavelet
                    [ws,what] = haarwavelet(BN.Nvals(cc));
                elseif (strcmp(BN.WaveType,'mean')) % sliding mean wavelet
                    [ws,what] = meanwavelet(BN.Nvals(cc));
                else
                    error('Wavetype %s is invalid. Must be svg, haar, or mean')
                end
                BN.Ws{cc} = ws;
                BN.Whats{cc} = what;
            end
            
        end
        
        function BN = analyzeTracks(BN,tracklist,varargin)
            % run BNEW analysis on tracklist, using preset coefficients to get
            % corrected MSD for each span   
            % Each element of tracklist is an Nx2 particle trajectory,
            % where N is the number of time points. The nmax datapoints at
            % the start and end of the track are not used in calculating
            % the corrected MSD (though they are used in the wavelet
            % smoothing)
            % trajectories of length < 2*nmax are not included in the
            % analysis
            % time points are assumed to be equispaced

            [BN.MSD, BN.Cnt, BN.Sterr,BN.SmVel] = BNEWanalyzeTracks(tracklist,BN.Nvals,BN.Ws,varargin{:});
            
            % time points
            for nc = 1:length(BN.Nvals)
                BN.Tbin{nc} = 1:length(BN.MSD{nc});
            end
        end
        
        function plotMSD(WL,varargin)
            % plot corrected MSD results post-wavelet analysis
            
            % optional arguments
            % type of plot to use (default is loglog)
            doplot = @loglog;
            
            % color map for the plotting
            cmap = lines(length(WL.Nvals));
            % type of line to use for the plotting
            dottype = '.-';
            % additional arguments for the plotting
            plotargs = {};
            
            for vc = 1:2:length(varargin)
                switch (varargin{vc})
                    case('cmap')
                       cmap=varargin{vc+1};     
                    case('doplot')
                        doplot = varargin{vc+1};
                    case('dottype')
                        dottype = varargin{vc+1};
                    case('plotargs')
                        plotargs = varargin{vc+1};
                end
            end                        
            
            for nc = 1:length(WL.Nvals)
                tbinall = WL.Tbin{nc};
                MSDtot = WL.MSD{nc};    
      
                doplot(tbinall,MSDtot,dottype,'Color',cmap(nc,:),plotargs{:});
              
                hold all   
            end
            hold off
        end
        
        function BN = rescaleData(BN,varargin)
            % rescale corrected MSD results
            
            % recalculate rescaling coefficients rather than using saved
            % ones
            recalc = 1;
            
            for vc = 1:2:length(varargin)
                switch(varargin{vc})
                    case('recalc')
                        recalc = varargin{vc+1};
                end                
            end                       
                    
            if (recalc)
                for nc = 1:length(BN.Nvals)                
                    nn = BN.Nvals(nc);    
                    [BN.Tscaled{nc},BN.MSDscaled{nc},A,B] = rescaleMSD(BN.Tbin{nc},BN.MSD{nc},BN.Whats{nc},nn);        
                    BN.Afunc{nc}(BN.Tbin{nc}) = A;
                    BN.Bfunc{nc}(BN.Tbin{nc}) = B;
                end
            else
                % use precalculated rescaling coefficients                
                for nc = 1:length(BN.Nvals)
                   BN.Tscaled{nc} = BN.Afunc{nc}(BN.Tbin{nc})./BN.Bfunc{nc}(BN.Tbin{nc});
                   BN.MSDscaled{nc} = BN.MSD{nc}./BN.Bfunc{nc}(BN.Tbin{nc});
                end
            end                
        end
        
        function [allxvals, allyvals] = plotRescaled(BN,varargin)
            % plot rescaled data
            
            % default coloring is lines spectrum
            % plotting color matrix (different for each n)
            cmat = [];
            
            % additional plotting arguments
            % (eg: constant color, line width)
            plotargs={};
            
            % plot fitted line if defined
            plotfit = 1;
            % types of dots/lines to use for plotting
            dottype = '.-';
            % which spans to save into allyvals and allxvals and to plot
            nvalind = 1:length(BN.Nvals);
            % how far out in k to plot (up to floor(n*kmax))
            kmax = 1;
            
            % type of plot (default is linear; can specify loglog, etc)
            doplot = @plot;
            
            % additional arguments for plotting the fitted line
            fitplotargs = {};
            
            
            % scaling factor for MSD
            msdscl = 1;
            % scaling factor for rescaled time
            xscl = 1;
            
            for vc = 1:2:length(varargin)
                switch(varargin{vc})
                    case('cmat')
                        cmat = varargin{vc+1};
                    case('plotfit')
                        plotfit = varargin{vc+1};
                    case('doplot')
                        doplot = varargin{vc+1};
                    case('dottype')
                        dottype = varargin{vc+1};
                    case('plotargs')
                        plotargs = varargin{vc+1};
                    case('nvalind')
                        nvalind = varargin{vc+1};
                    case('kmax')
                        kmax = varargin{vc+1};
                    case('fitplotargs')
                        fitplotargs=varargin{vc+1};
                    case('msdscl')
                        msdscl = varargin{vc+1};
                    case('xscl')
                        xscl = varargin{vc+1};
                end
            end
            
            % set up color matrices
            if (isempty(cmat))
                cmat = lines(length(nvalind));
            end
            
            allxvals = [];
            allyvals = [];
            
            mint = inf; maxt = -inf;
            
            ct = 0;
            for nc = nvalind
                
                nn = BN.Nvals(nc);
                
                if isempty(BN.Tscaled{nc})
                    % no data, no plotting
                    continue
                end
                
                xval = BN.Tscaled{nc};
                
                MSDtot = BN.MSDscaled{nc};
                %sterrtot = BN.Sterr{nc};
                %cnttot= BN.Cnt{nc};
                
                if (kmax>0)
                    nend = min(floor(nn*kmax), numel(xval));
                else
                    [~,nend] = max(xval);
                end
                mint = min([mint,xval(1:nend)]);
                maxt = max([maxt,xval(1:nend)]);
                
                doplot(xval(1:nend)*xscl,MSDtot(1:nend)*msdscl,dottype,'Color',cmat(nc,:),plotargs{:})
                drawnow
                hold all
                
                % if (ismember(nc,nvalind))
                allxvals = [allxvals, xval(1:nend)];
                allyvals = [allyvals, MSDtot(1:nend)];
                % end
                ct = ct+1;
                labels{ct} = sprintf('n=%d',nn);
            end
            hold off
            
            
            if (~isempty(BN.Dfit) & plotfit)
                tlist = linspace(mint,maxt);
                
                if (BN.FitLog)
                    fitvals = exp(BN.FitFunc(BN.Dfit,log(tlist)));
                else
                    fitvals = BN.FitFunc(BN.Dfit,tlist);
                end
                
                hold all
                plot(tlist*xscl,fitvals*msdscl,'k',fitplotargs{:})
                hold off
                
                ct = ct+1;
                labels{ct} = 'fit';
                
%                 % write out fit parameters
%                 tmin = min(tlist*xscl);
%                 mmax = max(fitvals*msdscl);
%                 mmin = min(fitvals*msdscl);
%                 
%                 hold all
%                 text(tmin,mmax-(mmax-mmin)*0.1,sprintf('D=%0.2f\n eps=%0.2f \n alpha=%0.2f',...
%                     BN.Dfit(1)*msdscl/xscl,BN.Dfit(2)*sqrt(msdscl),BN.Dfit(3)))
%                 hold off;
            end
            
            xlabel('rescaled time')
            ylabel('rescaled MSD')
            legend(labels)
        end
        
        function [BN,Dfit,CovFit,R2] = fitDcoeff(BN,options)
            % fit coefficients to rescaled results of BNEW analysis
            
            % outputs:
            % Dfit = [D,epsilon,alpha]
            % CovFit = covariance matrix of fitted coefficients (doesn't
            % mean much since assumes independent data point noise)
            % R2 is sum of residuals squared
            
            % default options:
            opt = struct();
            
            % indices of which spans to use for fitting
            opt.nvalind = 1:length(BN.Nvals);
            % fit log-log function
            opt.fitlog = 0;
            % kmax: use k values 1:floor(kmax*n) from each span for fitting
            % kmax<0 means use up to the maximum in that
            opt.kmax = 1;
            

            % weights for the data points
            % weight = 0: unweighted
            % weight =1: use data counts
            % weight = 2: use standard error
            opt.weight = 0; 
            % alternately, specify explicitly a list of weights
            opt.weightlist=[];
            
            % supply function to use for fitting
            opt.fitfunc=[];
            
            % starting point for fitting
            opt.fitd0 = [0.1 0.1 1];           
            
            % input options
            if (nargin>1)
                opt = copyStruct(options,opt);
            end
            
            
            % concatenate all data
            alldatatbin = [];
            alldataMSD = [];
            alldataweights = [];
            for nc = opt.nvalind
                if isempty(BN.Tscaled{nc})
                    continue
                end
                xval = BN.Tscaled{nc};
                MSDtot = BN.MSDscaled{nc};
                sterrtot = BN.Sterr{nc};
                cnttot= BN.Cnt{nc};
                
                if (opt.kmax<0)
                    [~,nend] = max(xval);
                else
                    %To fix out-of-range error. Not sure what effect this
                    %will have?
                    warning("BNEW.fitDcoeff: Using end-of-range value for nend")
                    nend = min(floor(BN.Nvals(nc)*opt.kmax), numel(xval));
                end
                
                alldatatbin = [alldatatbin xval(1:nend)];
                alldataMSD = [alldataMSD MSDtot(1:nend)];
                if (opt.weight==1)
                    alldataweights = [alldataweights sqrt(cnttot(1:nend))];
                elseif (opt.weight==2)
                    alldataweights = [alldataweights 1./sterrtot(1:end)];
                end
            end
            
            BN.FitLog = opt.fitlog;
            
            ind = 1:length(alldatatbin);
            if (opt.weight)
                weights = {'Weights',alldataweights(ind).^2};
            else
                weights={};
            end
            
            if (~isempty(opt.weightlist))
                weights = {'Weights',opt.weightlist};
            end
            
            if (~isempty(opt.fitfunc))
                BN.FitFunc = opt.fitfunc;
            end           
            
            if (opt.fitlog)
                if (isempty(opt.fitfunc)); BN.FitFunc = @(D,lt) log(4*D(1)*exp(lt).^D(3)+4*D(2)^2); end
                [Dfit,R,J,CovFit,MSE,ErrorModelInfo] = nlinfit(log(alldatatbin(ind)),log(alldataMSD(ind)),BN.FitFunc,opt.fitd0,weights{:})
                
            else
                if (isempty(opt.fitfunc)); BN.FitFunc = @(D,t) 4*D(1)*t.^D(3)+4*D(2)^2; end
                [Dfit,R,J,CovFit,MSE,ErrorModelInfo] = nlinfit(alldatatbin(ind),alldataMSD(ind),BN.FitFunc,opt.fitd0,weights{:});
            end                        
            
            BN.Dfit = Dfit;
            BN.CovFit = CovFit;
            
            R2 = R'*R;
        end
        
        function BN = loadWaveletCoeff(BN,coeffile)
            % load in wavelet coefficients and rescaling functions from file
            % data files are made using savewaveletcoeff.m           
            
            load(coeffile,'wavetype','wavedeg','Afunc','Bfunc','ws','what')
            
            % check that wavelet type and degree match up
            if (~strcmp(BN.WaveType,wavetype))
                error('Wavelet coefficient file does not have the right wavelet type: %s %s', WL.WaveType,wavetype)
            end
            
            if (strcmp(wavetype,'svg') || strcmp(wavetype,'poly'))
                if (BN.WaveDeg ~= wavedeg)
                    error('Wavelet coefficient file does not have the right degree: %d %d', WL.WaveDeg, wavedeg)
                end
            end   
            
            for scc = 1:length(BN.Nvals)
                nn = BN.Nvals(scc);
                BN.Afunc{scc} = Afunc(nn,:);
                BN.Bfunc{scc} = Bfunc(nn,:);
                
                BN.Ws{scc} = ws{nn};
                BN.Whats{scc} = what{nn};
            end
            
        end
    end
end