function ThetaMatrix = computeThetaMatrix(estMdl, steps, fMatrix)
% Return Theta = PhiMatrix * factorMatrix as in Kilian (2017, p.111)
% 
% Input : estMdl (estimated varm object), steps ith, factor Matrix
% Output: Theta Matrix
%
% Reference: Lutkepohl (2005), and Kilian and Lutkepohl (2017).
% Author: Binh Thai Pham, PhD. (binhpham79@gmail.com). 2020.

if nargin < 3
    fMatrix = [];
end

if nargin < 2
    steps = 0;
end

if isobject(estMdl) && isa(estMdl,'varm') && (~isempty(estMdl.AR))

    [AComp, K, p] = getCompanion(estMdl);

    % if fMatrix is empty then assigns cholesky decomposition of Sigma
    % as default.

    if isempty(fMatrix) && (~isempty(estMdl.Covariance)) && (all(~isnan(estMdl.Covariance),'all'))
        fMatrix = chol(estMdl.Covariance, 'lower');
    else
        fMatrix = eye(K);
    end

    if steps == 0
        ThetaMatrix = eye(K) * fMatrix;
    else
        J = [eye(K) zeros(K, K*(p-1))];
        % Theta = Phi * P if P is a Cholesky decomposition.
        ThetaMatrix = J*(AComp^steps)*J' * fMatrix;
    end
    
    if all(isnan(estMdl.Covariance),'all')
        warning('Covariance matrix is not a number (NaN).');
    end
else
    ThetaMatrix = [];
end

end