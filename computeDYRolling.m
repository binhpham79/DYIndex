function [DYroll,spillroll,fromroll,toroll,netroll] = ...
    computeDYRolling(dataset, lags, steps, window, useGIRF,checkVAR)
% Return Diebold and Yilmaz Rolling Index.
% Syntax  [DYroll,spillroll,fromroll,toroll,netroll] = 
%                   computeDYRolling(dataset, lags, steps, window, useGIRF)
%
% Input: dataset = data matrix, 
%        lags    = VAR(default p = 1), 
%        steps   = h-FEVD (default = 10),
%        window  = rolling window (default = 200), 
%        useGIRF = 1 (default, Generalized) or 0 (Cholesky)
%        checkVAR= 1 (default, check VAR stability), 0 otherwise
%
% Output: Diebold and Yilmaz (2009, 2012) Rolling Spillover Index
%        
%        DYroll    = (nobs,nvars+1,nvars+1);
%        spillroll = (nobs,1);
%        fromroll  = (nobs,nvars);
%        toroll    = (nobs,nvars);
%        netroll   = (nobs,nvars);
%
% Note: DY(2009,2012) has a structure
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
% Reference: Diebold and Yilmaz (2009,2012,2014).
% Author: Binh Thai Pham, PhD. (binhpham79@gmail.com). 2020.

if nargin < 6
    checkVAR = 1;
end

if nargin < 5 
    useGIRF = 1;
end

if nargin < 4
    window = 200;
end

if nargin < 3
    steps = 10;
end

if nargin < 2
    lags = 1;
end

nvars  = size(dataset,2);
nobs   = size(dataset,1);
rend   = size(dataset,1);
rstart = 1;
jump   = 1;

error_window = round(window/10);

Mdl = varm(nvars,lags);

DYroll    = nan(nobs,nvars+1,nvars+1);
spillroll = nan(nobs,1);
fromroll  = nan(nobs,nvars);
toroll    = nan(nobs,nvars);
netroll   = nan(nobs,nvars);

for step=rstart:jump:(rend-window+1)
    
    estMdl = estimate(Mdl, dataset(step:(step+window-1),:));
    
    AR = [eye(estMdl.NumSeries) estMdl.AR];
    
    if checkVAR
        isStability = isStable(reflect(LagOp(AR)));
    else
        isStability = 1;
    end
    
    fprintf('Start = %d ; End = %d ; VAR stability = %d \n', step,step+window-1,isStability);
    
    if ~isStability
        estMdl = estimate(Mdl, dataset(step-error_window:(step+window-1),:));
        AR = [eye(estMdl.NumSeries) estMdl.AR];
        isStability = isStable(reflect(LagOp(AR)));
        fprintf('Start = %d ; End = %d ; Recheck VAR stability = %d \n', step-error_window,step+window-1,isStability);
    end
    
    if isStability
        pos = step+window-1;
        [   DYroll(pos,:,:), ...
            spillroll(pos), ...
            fromroll(pos,:), ...
            toroll(pos,:), ...
            netroll(pos,:)      ] = computeDYtable(estMdl,steps,useGIRF);
    else
        DYroll(pos,:,:) = nan;
        spillroll(pos) = nan;
        fromroll(pos,:) = nan;
        toroll(pos,:) = nan;
        netroll(pos,:) = nan;
    end
end    
end