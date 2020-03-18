function gIRF = computeGIRF(estMdl, steps, Sigma, asMATLAB)
% Return IRF given estimated 'varm' object and (if needed) a covariance matrix
% Syntax: gIRF = computeIRF(estMdl, steps, fMatrix, asMATLAB)
%
% Input : estMdl (estimated varm object), 
%         steps (horizons), 
%         Covariance Matrix (if needed)
%         as is MATLAB order flag

% Note  : Variable order is the same SeriesName property of estMdl object
%         Return empty matrix [] if fails to compute.
%
% Example: gIRF(:,3,2) to get response of var 3 given a impulse on var 2
%
% Note   : the order above is different from Matlab (steps, imp, res),
% however if asMATLAB = true then the return is as the same as MATLAB.
%
% Reference: Pesaran and Shin (1998).
% Author: Binh Thai Pham, PhD. (binhpham79@gmail.com). 2020.

if nargin < 4
    asMATLAB = 0;
end

if nargin < 3
    Sigma = [];
end

if nargin < 2
    steps = 10;
end

if isobject(estMdl) && isa(estMdl,'varm') && (~isempty(estMdl.AR))
    
    K  = estMdl.NumSeries;
    
    if isempty(Sigma)
        Sigma =estMdl.Covariance;
    end
    
    if all(steps <= 0, 'all') || isempty(steps)
        steps = 1;
    end
    
    gIRF = nan(steps,K,K);
    
    for k=1:K % impulse on kth var
        
        % selection vector
        ek = zeros(K,1);
        ek(k) = 1;
        scaleFactor = Sigma(k,k)^-0.5;
        
        for i=1:steps
            % gIRF of size (steps, response, impulse)
            gIRF(i,:,k) = scaleFactor*computePhiMatrix(estMdl, i-1)*Sigma*ek;
        end
        
    end
    
    % To get same as Matlab output
    if asMATLAB ~= 0
        gIRF = permute(gIRF,[1,3,2]);
    end
else
    gIRF = [];
end
    
end