function GMSPE = computeGMSPE(estMdl, steps, asMATLAB)
% Return GMSPE(horz h,var kth, shock jth) at each horizon h=0..steps-1
% 
% Each row is the contribution of jth-var shock (j=1..K) to each
% kth-variable.
%
% Input : estMdl (estimated varm object), horizon steps, asMATLAB flag
%
% Output: Generalized FEVD matrix GMSPE(h,k,j)
%
% Reference: Pesaran (2015, eq 24.26), and Pesaran and Shin (1998).
% Author: Binh Thai Pham, PhD. (binhpham79@gmail.com). 2020.

if nargin < 3
    asMATLAB = 0;
end

if nargin < 2
    steps = 2;
end

if isobject(estMdl) && isa(estMdl,'varm') && (~isempty(estMdl.AR))
    
    K = estMdl.NumSeries;
    
    GMSPE = zeros(steps,K,K);
    
    for i=1:steps
        GMSPE(i,:,:) = computeGMSPEkj(estMdl, i);
    end

    if asMATLAB ~= 0
        GMSPE = permute(GMSPE, [1,3,2]);
    end
else
    GMSPE = [];    
end

end