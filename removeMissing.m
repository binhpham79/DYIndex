function [cleandata, index] = removeMissing(data)
% Synex [cleandata, index] = removeMissing(data)
% Remcve missing values in any row. Return the clean dataset and the
% selection index which then can be used to select non-missing row from the
% original data.
%
% Author: Pham Thai Binh, PhD. 2020.

    index = all(~ismissing(data),2);
    cleandata = data(index,:);
    
end