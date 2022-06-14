function Cells = loadCellObj(dirname,filenameglobs,excludestr)
% load in an array of cell objects from multiple files
% dirname is the containing directory
% filenameglobs is a cell array of glob strings for the file names
% optionally, exclude all files containing pattern excludestr

Cells = CellObj.empty;
ncell=0;

exclude = {};
if (nargin>2)
    exclude = excludestr;
end

for filec = 1:length(filenameglobs)
    globname = filenameglobs{filec};
   
    
    % list all cells starting with this prefix
    cellfiles = dir([dirname globname]);    
    for cc = 1:length(cellfiles);
        excludefile=0;
        for ec = 1:length(exclude)
            if (strfind(cellfiles(cc).name,exclude{ec}))
                excludefile = 1; break
            end
        end
        if (excludefile); continue; end
        try
            load([dirname cellfiles(cc).name],'CL')
        catch
            warning('Failed to load %s',cellfiles(cc).name)
            continue
        end
        ncell = ncell+1;
        Cells(ncell) = CL;        
    end
end