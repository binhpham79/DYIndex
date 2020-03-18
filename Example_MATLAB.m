% Example 1

load data_JDanish;

% create new varm object
K = 4;
p = 2;
estMdl = varm(K, p);
estMdl.SeriesNames = DataTable.Properties.VariableNames;

estMdl = estimate(estMdl, DataTable.Series);

%% IRF

% compute othogonalized IRF using MATLAB function
estIRF = irf(estMdl, 'Method', 'orthogonalized', 'NumObs', 100);

% Or use computeIRF function
IRF = computeIRF(estMdl, 100, [], true);

% compute generalized IRF using MATLAB function
estGIRF = irf(estMdl, 'Method', 'generalized', 'NumObs', 100);

% Or use computeGIRF function
gIRF = computeGIRF(estMdl, 100, [], true);


% In MATLAB IRF and FEVD return a matrix of (horz, impulse, response)
% so that IRF(:,2,3) refers to responses of var 3 w.r.t impulse on var 2.
% plot Orthogonalized IRF

armairf(estMdl.AR,[], 'InnovCov', estMdl.Covariance, 'NumObs', 100);

% plot Generalized IRF
armairf(estMdl.AR,[], 'InnovCov', estMdl.Covariance, 'Method', 'generalized', 'NumObs', 100);


%% FEVD
%

% compute orthogonalized FEVD
oFEVD = fevd(estMdl, 'Method', 'orthogonalized', 'NumObs', 30);
% plot FEVD from MATLAB
armafevd(estMdl.AR,[], "InnovCov", estMdl.Covariance, 'Method', 'orthogonalized', 'NumObs', 30);

% compute generalized FEVD
gFEVD = fevd(estMdl, 'Method', 'generalized', 'NumObs', 30);
% Plot g-FEVD from MATLAB
armafevd(estMdl.AR,[], "InnovCov", estMdl.Covariance, 'Method', 'generalized', 'NumObs', 30);

% compute FEVD use function 
oMSPE = computeMSPE(estMdl, 30);

gMSPE = computeGMSPE(estMdl, 30);