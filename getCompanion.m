function [AComp, K, p] = getCompanion(estMdl)
% Return the companion matrix of an estimated varm object (estMdl).
% Return [] (empty) if estMdl is not a valid estimated varm object.
%
% Input : estMdl (estimated varm object)
% Output: [Acomp, K, p] = [companion matrix, NumSeries, Lag p]
%
% Reference: Lutkepohl (2005), and Kilian and Lutkepohl (2017).
% Author: Binh Thai Pham, PhD. (binhpham79@gmail.com). 2020.

if isobject(estMdl) && isa(estMdl,'varm')
    
    K = 0;
    p = 0;
    
    if ~isempty(estMdl.AR)
        
        K = estMdl.NumSeries;
        p = estMdl.P;
        
        
        AComp = [];
        for i = 1:p
            AComp = [AComp estMdl.AR{i}];
        end
        AComp =[AComp; [kron(eye(p-1), eye(K)) zeros(K*(p-1), K)]];
        
    else
        AComp = [];
    end
    
else
    AComp = [];
end
        
end