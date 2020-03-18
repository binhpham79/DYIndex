function PhiMatrix = computePhiMatrix(estMdl, steps)
% Return Phi = J*A^(i)*J' as in Lutkepohl (2005, eq 2.1.17)
% 
% Input : estMdl (estimated varm object), steps ith
% Output: Phi Matrix
%
% Reference: Lutkepohl (2005), and Kilian and Lutkepohl (2017).
% Author: Binh Thai Pham, PhD. (binhpham79@gmail.com). 2020.

if nargin < 2
    steps = 0;
end

if isobject(estMdl) && isa(estMdl,'varm') && (~isempty(estMdl.AR))
    
    [AComp, K, p] = getCompanion(estMdl);

    % Selection matrix J
    J = [eye(K) zeros(K, K*(p-1))];
    PhiMatrix = J*(AComp^steps)*J';
    
else
    PhiMatrix = [];
end

end