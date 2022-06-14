% convert old data for multiple cells into Cell objects saved in individual
% files


%% list of file names to convert
dirname = '~/data/Caleb/lysophase20140821/'
matfiles={};
cellct = 0;
for cnum=[5,6]
    cellct = cellct+1;
    matfiles{cellct} = [dirname sprintf('lysoph%d/lysoph%d_segcont.mat',cnum,cnum)]    
end

%%
dirname = '~/data/Caleb/lysophase20141113/'
matfiles={};
cellct = 0;
for cnum=[1,2,4,5]
    cellct = cellct+1;
    if (cnum==1)
         matfiles{cellct} = [dirname sprintf('lysoctrl%d_results.mat',cnum)];
    else
        matfiles{cellct} = [dirname sprintf('lysoctrl%d_segcont.mat',cnum)];
    end
end

%% convert cells
outputdir = '~/data/Cells_lyso/';
for cc = 1
    matfile = matfiles{cc};
    disp(sprintf('Converting data for Cell %d: %s', cc,matfile))
    [Cells,cellmv,tracklist] = cellObjFromData(matfile,'outputdir',outputdir);
    disp(sprintf('%d Cells Found', length(Cells)))
    %delete(Cells)
end

%% load in all cells
outputdir = '~/data/Cells_lyso/';
fileglobs;
for filec = 1:cellct
    matfile = matfiles{filec};
    load(matfile,'prefix')
    fileglobs{filec} = [prefix '*.mat'];
end

Cells = loadCellObj(outputdir,fileglobs);

%% do wavelet analysis
maxntime = 2; % maximum wavelet span in sec

for ccell = 1:2%length(Cells)
    CL = Cells(ccell);
    
    disp(sprintf('Wavelet analysis on Cell %d', CL.Name))
    dt = mean(diff([CL.Frames.time]));
    
    dn = round(0.1/dt);
    nvals = 4:dn:round(maxntime/dt);
    
    if (isempty(CL.Wavelet))
        CL.runWaveletAnalysis(nvals);
        save([outputdir,CL.Name,'.mat'],'CL')
    end
end