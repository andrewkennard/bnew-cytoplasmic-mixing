function ind = getCellByName(Cells,str)
    tmp = find(arrayfun(@(x) strcmp(x.Name,str), Cells));
    if (isempty(tmp))
        ind = NaN;
    else
        ind = tmp(1);
    end
end