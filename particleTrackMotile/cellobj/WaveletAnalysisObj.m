classdef WaveletAnalysisObj
    % class defining wavelet analysis results and methods
    properties 
        % spans of wavelet
        Nvals = [];
        % Coefficients of position points
        Ws = {};
        % Coefficients of velocities (for each span)
        Whats = {};
        % type of wavelet (svg,haar,poly)
        WaveType='svg';
        % degree of velocity smoothing for polynomial wavelet
        WaveDeg=2;
        % number of tracks analyzed
        NTrack = {};
        % Which span to use for saving smoothed velocities
        SaveSmVel=0;
        % corrected mean squared displacements (for each span)
        MSD = {};
        % times (as frame numbers) associated with MSDs
        Tbin = {};
        % number of data points in each MSD calculation
        Cnt = {};
        % standard error of each MSD calculation
        Sterr = {};
        % smoothed velocities for each track
        SmVel = []
        % smoothed velocity autocorrelation
        SmVelCorr = {};        
        SmVelCorrCnt = {};
        SmVelCorrErr = {};
        
        % rescaled times, MSD
        Tscaled = {};
        MSDscaled = {};
        % two functions used for rescaling
        % expect MSD = D*A + eps*B
        Afunc = {}; 
        Bfunc = {};      
        
        % function for fitting diffusion coefficient/scaling
        FitFunc=[];
        % Diffusion coefficient, localization error, scaling
        Dfit = [];
        % covariance matrix for fitting;
        CovFit = [];
        
        % fit was done for log-log data
        FitLog = 0;
    end

    methods
        function WL = WaveletAnalysisObj(nvals,varargin)
            %nvals
            WL.Nvals = nvals;
            
            % optional arguments
            WL.WaveType='svg';
            WL.WaveDeg=2;
            WL.SaveSmVel = length(nvals);
            
            for vc = 1:2:length(varargin)
                switch (varargin{vc})
                    case('savesmvel')
                        WL.SaveSmVel=varargin{vc+1};
                    case('wavetype')
                        WL.WaveType=varargin{vc+1};
                    case('wavedeg')
                        WL.WaveDeg=varargin{vc+1};
                    case('Afunc')
                        WL.Afunc = varargin{vc+1};
                    case('Bfunc')
                        WL.Bfunc = varargin{vc+1};
                end
            end                        
        end
        
        function WL = getCoefficients(WL)
            % get wavelet coefficients
            for cc = 1:length(WL.Nvals)
                if (strcmp(WL.WaveType,'svg'))
                    [ws, what] = getoptwavelet_sgolay_pos(WL.Nvals(cc),WL.WaveDeg+1);
                elseif (strcmp(WL.WaveType,'poly'))
                    [ws, what] = getoptwavelet_poly(WL.Nvals(cc),inf,WL.WaveDeg);
                elseif (strcmp(WL.WaveType,'haar'))
                    [ws,what] = haarwavelet(WL.Nvals(cc));
                elseif (strcmp(WL.WaveType,'exp'))
                    [ws, what] = getoptwavelet_tau(WL.Nvals(cc),WL.WaveDeg,1,1);
                elseif (strcmp(WL.WaveType,'const')) 
                    % wavelet that minimizes velocity error, assuming
                    % constant velocity, and locE/D = WaveDeg
                    [ws,what] = getoptwavelet(WL.Nvals(cc),WL.WaveDeg);
                else
                    error('Wavetype %s is invalid. Must be svg, poly, or haar or const')
                end
                WL.Ws{cc} = ws;
                WL.Whats{cc} = what;
            end

        end
        
        function WL = analyzeTracks(WL,tracklist,varargin)
            % run wavelet analysis on tracklist, using preset coefficients to get
            % corrected MSD for each span            

            [WL.MSD, WL.Cnt, WL.Sterr,~,~,~,WL.SmVel] = waveletAnalyzeTracks(tracklist,WL.Nvals,WL.Ws,varargin{:});
            %[WL.MSD, WL.Cnt, WL.Sterr] = waveletAnalyzeTracks(tracklist,WL.Nvals,WL.Ws,varargin{:});
            for nc = 1:length(WL.Nvals)
                WL.Tbin{nc} = 1:length(WL.MSD{nc});
            end
        end
        
        function plotMSD(WL,varargin)
            % plot corrected MSD results post-wavelet analysis
            
            cmap = lines(length(WL.Nvals));
            dottype = '.-';
            plotargs = {};
            
            for vc = 1:2:length(varargin)
                switch (varargin{vc})
                    case('cmap')
                       cmap=varargin{vc+1};                  
                    case('dottype')
                        dottype = varargin{vc+1};
                    case('plotargs')
                        plotargs = varargin{vc+1};
                end
            end                        
            
            for scc = 1:length(WL.Nvals)
                tbinall = WL.Tbin{scc};
                MSDtot = WL.MSD{scc};    
                loglog(tbinall,MSDtot,dottype,'Color',cmap(scc,:),plotargs{:})
                hold all   
            end
            hold off
        end

        function WL = loadWaveletCoeff(WL,coeffile,varargin)
            % load in wavelet coefficients and rescaling functions from
            % file
            % data files are made using savewaveletcoeff.m
            % if loadeddata keyword is supplied, then don't bother loading
            % from file
            % instead use the given values of Afunc, Bfunc,             
            
            load(coeffile,'wavetype','wavedeg','Afunc','Bfunc','ws','what')
            
            % check that wavelet type and degree match up
            if (~strcmp(WL.WaveType,wavetype))
                error('Wavelet coefficient file does not have the right wavelet type: %s %s', WL.WaveType,wavetype)
            end
            
            if (strcmp(wavetype,'svg') || strcmp(wavetype,'poly'))
                if (WL.WaveDeg ~= wavedeg)
                    error('Wavelet coefficient file does not have the right degree: %d %d', WL.WaveDeg, wavedeg)
                end
            end   
            
            for scc = 1:length(WL.Nvals)
                nn = WL.Nvals(scc);
                WL.Afunc{scc} = Afunc(nn,:);
                WL.Bfunc{scc} = Bfunc(nn,:);
                
                WL.Ws{scc} = ws{nn};
                WL.Whats{scc} = what{nn}
            end
            
        end
        
        function WL = rescaleData(WL,varargin)
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
                for scc = 1:length(WL.Nvals)                
                    nn = WL.Nvals(scc);    
                    [WL.Tscaled{scc},WL.MSDscaled{scc},B,A] = rescaleMSD(WL.Tbin{scc},WL.MSD{scc},WL.Whats{scc},nn);        
                    WL.Afunc{scc}(WL.Tbin{scc}) = A;
                    WL.Bfunc{scc}(WL.Tbin{scc}) = B;
                end
            else
                % use precalculated rescaling coefficients                
                for scc = 1:length(WL.Nvals)
                   WL.Tscaled{scc} = WL.Afunc{scc}(WL.Tbin{scc})./WL.Bfunc{scc}(WL.Tbin{scc});
                   WL.MSDscaled{scc} = WL.MSD{scc}./WL.Bfunc{scc}(WL.Tbin{scc});
                end
            end                
        end
        
        function [allxvals, allyvals] = plotRescaled(WL,varargin)
            % plot rescaled data
            
            % default coloring is lines spectrum
            % cmat supercedes color specification
            % plotting color (constant for all n)
            color='';
            % plotting color matrix (different for each n)
            cmat = [];
            
            % additional plotting arguments
            plotargs={};
            
            % scaling for x axis
            dt = 0;
            % plot fitted line if defined
            plotfit = 1;
            % types of dots/lines to use for plotting
            dottype = '.-';
            % which spans to save into allyvals and allxvals
            savevals = 1:length(WL.Nvals);
            % how far out in k to plot (relative to 2*n+1)
            nendscl = 1;
            
            % include plot on loglog axes
            logplot = 1;
            % color of fitted line
            fitlinecolor = 'k';
            
            % scaling factor for MSD
            msdscl = 1;
            % scaling factor for rescaled time
            xscl = 1;
            
            % width for fitted line
            fitlinewidth=2;
           
            for vc = 1:2:length(varargin)
                switch(varargin{vc})
                    case('color')
                        color = varargin{vc+1};                        
                    case('cmat')
                        cmat = varargin{vc+1};
                    case('dt')
                        dt = varargin{vc+1};
                    case('plotfit')
                        plotfit = varargin{vc+1};
                    case('dottype')
                        dottype = varargin{vc+1};
                    case('plotargs')
                        plotargs = varargin{vc+1};
                    case('savevals')
                        savevals = varargin{vc+1};
                    case('nendscl')
                        nendscl = varargin{vc+1};
                    case('logplot')
                        logplot = varargin{vc+1};
                    case('fitlinecolor')
                        fitlinecolor=varargin{vc+1};
                    case('msdscl')
                        msdscl = varargin{vc+1};
                    case('xscl')
                        xscl = varargin{vc+1};
                    case('fitlinewidth')
                        fitlinewidth=varargin{vc+1};
                end
            end
            
            % set up color matrices
            if (isempty(cmat) & isempty(color))                
                cmat = lines(length(WL.Nvals));
            elseif (isempty(cmat))
                for scc = 1:length(WL.Nvals)
                    cmat(scc,:) = color;
                end                
            end
            
            if (logplot)
                nplot=2;
            else
                nplot=1;
            end
            
            allxvals = [];
            allyvals = [];
            
            mint = inf; maxt = -inf;

            for scc = savevals%1:length(WL.Nvals)

                nn = WL.Nvals(scc);
                
                if isempty(WL.Tscaled{scc})
                    % no data, no plotting
                    continue
                end
                
                if (dt>0)
                    xval = WL.Tscaled{scc}*dt;
                else
                    xval = WL.Tscaled{scc};
                end
                MSDtot = WL.MSDscaled{scc};
                sterrtot = WL.Sterr{scc};
                cnttot= WL.Cnt{scc};
                             
                
                if (nendscl>0)
                    nend = round((2*nn+1)*nendscl);
                else
                    [~,nend] = max(xval);
                end
                mint = min([mint,xval(1:nend)]);
                maxt = max([maxt,xval(1:nend)]);                                
                
                if (logplot)
                subplot(1,nplot,1)                  
                loglog(xval(1:nend)*xscl,MSDtot(1:nend)*msdscl,dottype,'Color',cmat(scc,:),plotargs{:})

                drawnow
                hold all
                end
                
                subplot(1,nplot,nplot)     
                plot(xval(1:nend)*xscl,MSDtot(1:nend)*msdscl,dottype,'Color',cmat(scc,:),plotargs{:})
                drawnow
                hold all   
                
                if (ismember(scc,savevals))
                    allxvals = [allxvals, xval(1:nend)];
                    allyvals = [allyvals, MSDtot(1:nend)];
                end
            end            
            for p = 1:nplot
                subplot(1,nplot,p)
                hold off
            end                        
            
            
            if (~isempty(WL.Dfit) & plotfit)
                tlist = linspace(mint,maxt);
                
                if (WL.FitLog)
                    fitvals = exp(WL.FitFunc(WL.Dfit,log(tlist)));
                else
                    fitvals = WL.FitFunc(WL.Dfit,tlist);
                end
                
                for p = 1:nplot
                    subplot(1,nplot,p)
                    hold all
                    plot(tlist*xscl,fitvals*msdscl,'LineWidth',fitlinewidth,'Color',fitlinecolor)
                    hold off
                end
                
            end
        end
        
        function [WL,Dfit,CovFit,binvals,R2] = fitDcoeff(WL,options)
            % fit Dfit coefficients to precalculated rescaled results
            % R2 is sum of residuals squared
            
            % default options
            opt = struct();
            
            % indices of which spans to use
            opt.nvalind = 1:length(WL.Nvals);
            % fit log-log function
            opt.fitlog = 0;
            % how far in delta t to go for fitting (relative to full span 2*n+1)
            opt.nendscl = 1;
            
            % fit to rebinned data
            opt.rebin = 0;
            binvals = {};
            % number of bins to use
            opt.nbin = 20;
            
            % do a linear fit (force alpha=1)
            opt.linfit = 0;
            opt.weight = 1; % use data counts for weighted fits
            opt.weightlist=[];
            
            % supply function to use for fitting
            opt.fitfunc=[];
            
            % starting points for nonlinear fit
            opt.fitd0 = [0.1 0.1 1];
            
            % rescale both t and msd
            opt.bscl = 1;
            
            % input options
            if (nargin>1)
                opt = copyStruct(options,opt);            
            end
            
            
            % concatenate all data
            alldatatbin = [];
            alldataMSD = [];
            alldataweights = [];
            for scc = opt.nvalind
                if isempty(WL.Tscaled{scc})
                    continue
                end
                xval = WL.Tscaled{scc};
                MSDtot = WL.MSDscaled{scc};
                sterrtot = WL.Sterr{scc};
                cnttot= WL.Cnt{scc};
                
                if (opt.nendscl<0)
                    [~,nend] = max(xval);
                else
                    nend = floor((2*WL.Nvals(scc))*opt.nendscl);
                end
                %nend = WL.Nvals(scc);
                alldatatbin = [alldatatbin xval(1:nend)];
                alldataMSD = [alldataMSD MSDtot(1:nend)];
                if (opt.weight==1)
                    alldataweights = [alldataweights sqrt(cnttot(1:nend))];
                elseif (opt.weight==2)
                    alldataweights = [alldataweights 1./sterrtot(1:end)];
                end
            end
            
            WL.FitLog = opt.fitlog;
            
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
                WL.FitFunc = opt.fitfunc;
            end
                        
            alldatatbin = alldatatbin/opt.bscl;
            alldataMSD = alldataMSD/opt.bscl;
            
            if (opt.linfit)
                % do a linear fit (force alpha=1)
                if (isempty(opt.fitfunc)); WL.FitFunc = @(D,t) 4*D(1)*t+2*D(2)^2; end
                [Dfit,R,J,CovFit,MSE,ErrorModelInfo] = nlinfit(alldatatbin(ind),alldataMSD(ind),WL.FitFunc,opt.fitd0(1:2),weights{:});
            elseif (opt.rebin)
                [xnew,ymean,ncnt,sterr] = rebindatamean(alldatatbin,alldataMSD,opt.nbin);
                ind = find(ncnt>0);
                xnew = xnew(ind); ymean = ymean(ind); ncnt = ncnt(ind); sterr = sterr(ind);
                binvals = [xnew', ymean', ncnt', sterr'];
                %loglog(xnew,ymean,'.-')
                %hold all
                %errorbar(xnew,ymean,sterr,'.-')
                %hold off
                if (isempty(opt.fitfunc)); WL.FitFunc = @(D,t) 4*D(1)*t.^D(3)+2*D(2)^2; end
                [Dfit,R,J,CovFit,MSE,ErrorModelInfo] = nlinfit(xnew,ymean,WL.FitFunc,opt.fitd0,'Weights',ncnt);
            else          
                if (opt.fitlog)
                    if (isempty(opt.fitfunc)); WL.FitFunc = @(D,lt) log(4*D(1)*exp(lt).^D(3)+2*D(2)^2); end
                    [Dfit,R,J,CovFit,MSE,ErrorModelInfo] = nlinfit(log(alldatatbin(ind)),log(alldataMSD(ind)),WL.FitFunc,opt.fitd0,weights{:})
                    
                else
                    if (isempty(opt.fitfunc)); WL.FitFunc = @(D,t) 4*D(1)*t.^D(3)+2*D(2)^2; end                 
                    [Dfit,R,J,CovFit,MSE,ErrorModelInfo] = nlinfit(alldatatbin(ind),alldataMSD(ind),WL.FitFunc,opt.fitd0,weights{:});
                end
            end
            
            Dfit(1) = Dfit(1)*opt.bscl^(Dfit(3)-1);
            Dfit(2) = Dfit(2)*sqrt(opt.bscl);
            
            WL.Dfit = Dfit;
            WL.CovFit = CovFit;
            
            R2 = R'*R;
        end
                
    end
end