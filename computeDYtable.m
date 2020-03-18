function [DYtable,spillover,fromvar,tovar,net] = computeDYtable(estMdl, steps, useGIRF)
% Return Diebold and Yilmaz Index table and other components.
% Syntax [DYtable,spillover,fromvar,tovar,net] = computeDYtable(estMdl, steps, useGIRF)
%
% Compute Diebold and Yilmaz (2009, 2012) Spillover Index
%
% DYtable has a structure
%
%    (i,j)    var1   var2   ...  ... varK   FROM
%     var1     x%                           Sum(DY(i,j),j <> 1)
%     var2                                  Sum(DY(i,j),j <> 2)
%     ... 
%     ...
%     varK                                  Sum(DY(i,j),j <> K)
%     TO     Sum(    ...    ...  ... i <> K SPILLOVER INDEX
%            DY(i,j), 
%            i <> 1)
%
% DYtable is of size (K+1) x (K+1)
% Column (K+1)th is FROM OTHERS; Row(K+1)th is TO OTHERS. Cell
% DYtable(K+1,K+1) is the SPILLOVER INDEX. While, Net=(TOVAR - FROMVAR).
%
% Reference: Diebold and Yilmaz (2009,2012,2014).
% Author: Binh Thai Pham, PhD. (binhpham79@gmail.com). 2020.

if nargin < 3
    useGIRF = 0;
end

if nargin < 2
    steps = 10;
end

if isobject(estMdl) && isa(estMdl,'varm') && (~isempty(estMdl.AR))
    
    K = estMdl.NumSeries;

    if useGIRF == 0
        
        oMSPE = 100*computeMSPEkj(estMdl,steps);

        fromvar = sum(oMSPE,2) - diag(oMSPE);
        tovar   = sum(oMSPE,1) - diag(oMSPE)';

        spillover = sum(tovar) / K;

        DYtable = [oMSPE fromvar; [tovar spillover]];
        net     = tovar - fromvar;
        
    else % if useGIRF ~= 0
        
        gMSPE = 100*computeGMSPEkj(estMdl,steps,true);
        
        fromvar = sum(gMSPE,2) - diag(gMSPE);
        tovar   = sum(gMSPE,1) - diag(gMSPE)';

        spillover = sum(tovar) / K;

        DYtable = [gMSPE fromvar; [tovar spillover]];
        tovar   = tovar';
        net     = tovar - fromvar;
    end
    
else
    DYtable = [];
end

end