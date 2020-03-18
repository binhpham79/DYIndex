function MSPE = computeMSPE(estMdl, steps, fMatrix, asMATLAB)
% Return MSPE = (steps, var kth, shock jth) at each horizon h=0..steps-1
% Syntax MSPE = computeMSPE(estMdl, steps, fMatrix, asMATLAB)
%
% Each row of squeeze(MSPE(h,:,:)) is the contribution of jth-var shock 
% (j=1..K) to each kth-variable, so that sum over column is of one. That is,
%
% 1 = MSPE(h,k,j=1)/MSPE(h,k,*) + MSPE(h,k,j)/MSPE(h,k,*) +...+ MSPE(h,k,j=K)/MSPE(h,k,*)
% for k = 1..K.
%
% Syntax: MSPE = computeMSPE(estMdl, steps, fMatrix, asMATLAB)
% Input : estMdl (estimated varm object), horizon steps, factor matrix
% (default or [] = Cholesky), asMATLAB=0 (default)
%
% Output: MSPE(k,j) matrix at horizon (steps) h
%
% Reference: Lutkepohl (2005), and Kilian and Lutkepohl (2017).
% Author: Binh Thai Pham, PhD. (binhpham79@gmail.com). 2020.

if nargin < 4
    asMATLAB = 0;
end

if nargin < 3
    fMatrix = [];
end

if nargin < 2
    steps = 2;
end

if isobject(estMdl) && isa(estMdl,'varm') && (~isempty(estMdl.AR))
    
    K = estMdl.NumSeries;
    MSPE = zeros(steps, K, K);
    MSPEkj = zeros(K);
    
    for i=1:steps
        
        % Note: function computeMSPEkj returns MSPEkj for h=0..steps-1 so that
        % we input the function at h = i in this function.
        % MSPEkj = computeMSPEkj(estMdl, i, fMatrix);
%         MSPEk  = repmat(sum(MSPEkj,2), [1 K]);
%         
%         MSPE(i,:,:) = MSPEkj ./ MSPEk;
        MSPE(i,:,:) = computeMSPEkj(estMdl, i, fMatrix);
        
    end
    
    if asMATLAB ~= 0
        MSPE = permute(MSPE, [1,3,2]);
    end
    
else
    MSPE = [];
end

end