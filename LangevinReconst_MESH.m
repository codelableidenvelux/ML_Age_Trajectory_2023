% LangevinReconst_MESH

% Function to compute nonparametric estimators for time series x
% This code was adapted from the R code of Brock & Carpenter 2012 PLoS one.
% Inputs:
%  x0 is the regressor
%  dx is the first difference of x0
%  nx is number of first differences
%  DT is time step
%  bw is the bandwidth for the kernel
%  na is number of a values for computing the kernel
%  avec is the mesh for the kernel

% December 2023. Adapted from R code by
% Enea Ceolini (Leiden University & QuantActions AG)
% Arko Ghosh (Leiden University)

% Notes:
% DT we set it to 1 in our case this corresponds to 1 day.
% bw = 0.3 * std(x0) -> according to Brock & Carpenter 2012 PLoS one paper.
% avec = linspace(min(x0), max(x0), length(x0) / 5) -> /5 needs to be
% double checked depending on the length x0
% na = length(avec)

function out = LangevinReconst_MESH(x0, dx, nx, DT, bw, na, avec)

  % Set up constants and useful preliminaries
  SF = 1/(bw*sqrt(2*pi));  % scale factor for kernel calculation
  x02 = x0.*x0; % second power of x
  dx2 = dx.*dx; % second power of dx
  % Compute matrix of kernel values
  Kmat = zeros(na, nx);
  for i = 1:nx   % loop over columns (x0 values)
    Kmat(:,i) = SF*exp(-0.5*(x0(i)-avec).*(x0(i)-avec)/(bw.*bw));
  end
  % Compute M1, M2, and sum of squares of x0 for each value of a
  M1.a = zeros(na,1);
  M2.a = zeros(na,1);
  mean.a = zeros(na,1);
  SS.a = zeros(na,1);
  for i = 1:na   % loop over rows (a values)
    Ksum = sum(Kmat(i, :));  % sum of weights
    M1.a(i) = (1/DT)*sum(Kmat(i,:)*dx')/Ksum;
    M2.a(i) = (1/DT)*sum(Kmat(i,:)*dx2')/Ksum;
    mean.a(i) = sum(Kmat(i,:)*x0(2:(nx+1))')/Ksum; % Buz removes 1/DT on 17 Nov 2011
    SS.a(i) = sum(Kmat(i,:)*x02(2:(nx+1))')/Ksum; % Buz removes 1/DT on 17 Nov 2011
  end
  % Compute conditional variance, diffusion and drift functions
  S2.x = SS.a - (mean.a.*mean.a); % sum of squares minus squared mean
  diff2.x = M2.a;
  mu.x = M1.a;
  % Return the following functions to the main program:
  %  mu.x is the drift function
  %  diff2.x is the diffusion function
  %  mean.a is the mean function used to compute the conditional variance
  %  S2.x is the conditional variance function

       out = struct('D1', mu.x', 'D2', diff2.x', 'D4',zeros(size(diff2.x')), 'ErrorD1', zeros(size(diff2.x')), 'ErrorD2', zeros(size(diff2.x')), 'N', ones(size(avec)), 'C', double(avec), ...
            'options', struct('datasize', size(x0), 'domain', [min(avec), max(avec)], 'bins', na, 'Tau', 1:5, 'dt', DT, 'method', 'MESH', 'title', 'MESH'));
end
