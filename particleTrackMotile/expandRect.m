function newrect = expandRect(rect,expf,cutx,cuty)
% given a rectangular region [xstart ystart width height]
% expand (or contract) around center by the given factor exp
% expanded if expf>1
% each linear dimension changes by factor expf
% cutx are cutoffs beyond which the rectangle may not extend in the x
% dimension

newrect(3:4) = rect(3:4)*expf;
newrect(1) = rect(1) - rect(3)*(expf-1)/2;
newrect(2) = rect(2) - rect(4)*(expf-1)/2;

% do not allow rectangle boundaries to exceed cutoff
if (nargin>2)
    if (newrect(1)<cutx(1))
        newrect(3) = newrect(3)-(cutx(1)-newrect(1)-1);
        newrect(1) = cutx(1);
    end
    if (newrect(1)+newrect(3)>cutx(2))
        newrect(3) = cutx(2)-newrect(1)+1;
    end
end
if (nargin>3)
    if (newrect(2)<cuty(1))
        newrect(4) = newrect(4)-(cuty(1)-newrect(2)-1);
        newrect(2) = cuty(1);
    end
    if (newrect(2)+newrect(4)>cuty(2))
        newrect(4) = cuty(2)-newrect(2)+1;
    end
end

end