function GMSPEkj = computeGMSPEkj(estMdl, steps, scale)
% Syntax GMSPEkj = computeGMSPEkj(estMdl, steps, scale)
%
% Return MSPE(var kth, shock jth) at the horizon h (0..steps-1)
% 
% Each row is the contribution of jth-var shock (j=1..K) to each kth-variable, 
% so that sum over column is of one. That is,
%
% Input : estMdl (estimated varm object), horizon steps, scale flag
%
% Output: GMSPE(k,j) matrix at horizon (steps) h
%
% Reference: Pesaran (2015, eq 24.26), and Pesaran and Shin (1998).
% Author: Binh Thai Pham, PhD. (binhpham79@gmail.com). 2020.

if nargin < 3
    scale = 0;
end

if nargin < 2
    steps = 1;
end

K = estMdl.NumSeries;
Sigma = estMdl.Covariance;

comp1 = zeros(K);
comp2 = zeros(K);
GMSPEkj = zeros(K);

for i=1:steps
    
    h = i-1;
    
    Phi_h = computePhiMatrix(estMdl,h);
    
    comp1 = comp1 + (Phi_h*Sigma) .* (Phi_h*Sigma); %
    comp2 = comp2 + (Phi_h * Sigma * Phi_h'); 
end


% e = eye(K);

for k=1:K
    for j=1:K
        % Formula 24.26 in Pesaran (2015, p.594)
        % GMSPEkj(k,j) = (Sigma(j,j)^-1) * (e(:,k)' * comp1 * e(:,j)) / comp2(k,k);
        GMSPEkj(k,j) = (Sigma(j,j)^-1) * comp1(k,j) / comp2(k,k);
    end
end

if scale == true
    scalemat= repmat(sum(GMSPEkj, 2), [1 K]);
    GMSPEkj = GMSPEkj ./ scalemat;
end

end