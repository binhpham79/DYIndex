function MSPEkj = computeMSPEkj(estMdl, steps, fMatrix)
% Syntax MSPEkj = computeMSPEkj(estMdl, steps, fMatrix)
% Return MSPE(var kth, shock jth) at the horizon h (0..steps-1)
% 
% Each row is the contribution of jth-var shock (j=1..K) to each kth-variable, 
% so that sum over column is of one. That is,

% 1 = MSPE(k,j=1,h)/MSPE(k,*,h) + MSPE(k,j,h)/MSPE(k,*,h) +...+ MSPE(k,j=K,h)/MSPE(k,*,h)
% for k = 1..K.
%
% Input : estMdl (estimated varm object), horizon steps, factor matrix
% (default or [] = Cholesky)
%
% Output: MSPE(k,j) matrix at horizon (steps) h
%
% Reference: Lutkepohl (2005), and Kilian and Lutkepohl (2017)
% Author: Binh Thai Pham, PhD. (binhpham79@gmail.com). 2020.

if nargin < 3
    fMatrix = [];
end

K = estMdl.NumSeries;

MSPEkj = zeros(estMdl.NumSeries);

for i=1:steps
    
    h = i-1;
    Theta_h = computeThetaMatrix(estMdl, h, fMatrix);
    Theta_h2 = Theta_h .* Theta_h;

    MSPEkj = MSPEkj + Theta_h2;
   
end

MSPEk  = repmat(sum(MSPEkj,2), [1 K]);
        
MSPEkj = MSPEkj ./ MSPEk;

end