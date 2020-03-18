function [IRF, GIRF] = computeIRF(estMdl, steps, fMatrix, getGIRF, asMATLAB)
% Return IRF given estimated 'varm' object and a factor matrix
% 
% Synatx: [IRF, GIRF] = computeIRF(estMdl, steps, fMatrix, getGIRF, asMATLAB)
%
% Input : estMdl (estimated varm object), 
%         steps (horizons), 
%         fMatrix if a factor matrix P (if needed), fMatrix=[] is Cholesky
%         getGIRF is a flag (true = also get Generalized IRF)
%         asMATLAB return matrices in the same order as MATLAB irf or fevd
%
% Output: IRF matrix of size (steps x j-th Response Var x k-th Shock)
%         This output is consistent with Lutkepohl textbook as the author
%         defines Phi(jk,h) or Theta(jk, h) corresponds to response of j-th var
%         to k-th shock at the h-th horizon. Here we have the IRF(h, j, k). 
%         Return empty matrix [] if not an estimated 'varm' object
%
% Note  : Variable order is the same as SeriesName property of estMdl object
%
% Example: gIRF(:,3,2) to get response of var 3 given a impulse on var 2
% Note   : the order above is different from Matlab (steps, imp, res),
% however if asMATLAB = true then the return is as the same as MATLAB.
% 
% Reference: Lutkepohl (2005), and Kilian and Lutkepohl (2017).
% Author: Binh Thai Pham, PhD. (binhpham79@gmail.com). 2020.

if nargin < 5
    asMATLAB = 0;
end

if nargin < 4
    getGIRF = 0;
end

if nargin < 3
    fMatrix = [];
end

if nargin < 2
    steps = 10;
end

if isobject(estMdl) && isa(estMdl,'varm') && (~isempty(estMdl.AR))

    K = estMdl.NumSeries;
    IRF = nan(steps,K,K);
    
    if all(steps <= 0, 'all') || isempty(steps)
        steps = 1;
    end
    
    % if fMatrix =[] then set it to Cholesky
    if isempty(fMatrix)
        fMatrix = chol(estMdl.Covariance,'lower');
    end
    
    % if also get Generalized IRF
    if getGIRF 
        gMatrix = estMdl.Covariance*diag(1./sqrt(diag(estMdl.Covariance)));
        GIRF = zeros(steps,K,K);
    else
        GIRF = [];
    end
    
    for i = 1:steps
        h = i-1;
        % IRF(i,:,:) = computeThetaMatrix(estMdl, i-1, fMatrix);
        Phi_h = computePhiMatrix(estMdl,h);
        IRF(i,:,:) =  Phi_h * fMatrix;
        GIRF(i,:,:) = Phi_h * gMatrix;
    end
    
    % To get same as Matlab output
    if asMATLAB ~= 0
        IRF = permute(IRF,[1,3,2]);
        GIRF= permute(GIRF,[1,3,2]);
    end
    
else
    IRF = [];
    GIRF= [];
end
end
