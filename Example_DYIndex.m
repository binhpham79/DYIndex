% Example Diebold and Yilmaz (2009,2012) Index
%
% Author: Pham Thai Binh. 2020.

load 'DYdata.mat'

iso2 = dyorder(:,(end-1):end);

% remove missing data

[dy, idx] = removeMissing(dy);

date_dy = datestr(idx,:);

%% DY estimation

startyear = 1991;
endyear   = 2019;
numobs =  (endyear - startyear + 1) * 12;

nvars = size(dy,2);
nlags = 12;
nsteps= 12;

useGIRF = 1;

Mdl = varm(nvars, nlags);

Mdl.SeriesNames = string(iso2);

% estimate VAR model
dy_sub = dy((end-numobs+1):end,:);
dyMdl = estimate(Mdl, dy_sub);

% Compute DY Table
[DYtable, Spillover, From, To, Net] = computeDYtable(dyMdl, nsteps, useGIRF);

DYtable = array2table([DYtable; [Net;NaN]']);

varnames = [iso2 repmat('         ',size(iso2,1),1);'FROM OTHERS'];
rownames = [iso2 repmat('         ',size(iso2,1),1);'TO OTHERS  ';'NET        '];

DYtable.Properties.VariableNames = cellstr(varnames);
DYtable.Properties.RowNames      = cellstr(rownames);

disp(DYtable);
       
%% Rolling DY

nlags = 3;
nsteps= 12;

window = 96;
[~,Spillroll_dy,~,~,~] = computeDYRolling(dy_sub,nlags,nsteps,window);

timeaxis = linspace(1991,2020,size(dy_sub,1));

figure;
plot(timeaxis,Spillroll_dy);
title('DY Rolling Index');
grid on

