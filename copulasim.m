function u = copulasim(family,theta,n,varargin)
% Simulates a sample u from a bivariate copula. The function uses the 
% conditioning and inversion technique for simulation.
%
% call: u = copulasim(family,theta,n[,parallel])
%
% input     family              - the copula family: 'gumbel', 'clayton', 
%                                 'frank', 't', 'gauss', 'ind', 'amhaq', 
%                                 'tawn', 'fgm', 'joe', 'plackett', 
%                                 'surclayton', 'surgumbel', 'surjoe'
%           theta               - copula parameters: for t-copula [rho, nu]
%           n                   - number of simulated points
%           parallel (optional) - switch parallelization on (=1) 
%                                 or off (=0); default: 0
%
% output    u 	- the simulated sample points
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

% some parsing
p = inputParser;
p.addRequired('family',@isstr);
p.addRequired('theta',@isvector);
p.addRequired('n',@isscalar);
p.addOptional('para',0,@isscalar);
p.parse(family,theta,n,varargin{:});

% sanity checks
if ~cpcheck(family,theta)
    error('copulasim:InvalidParameter',['invalid parameter for ',family,' copula']);
end

% determine parallelization mode
para = 0;
if nargin > 3
    if varargin{1} ~= 1 && varargin{1} ~= 0
        error('optional input argument parallel has to be 0 or 1')
    end
    para = varargin{1};
end

% initialize variables
w = rand(n,2);
u = zeros(n,2);

if para
    
    ppool = gcp; 
    poolsize = ppool.NumWorkers;

    w1 = w(:,1);
    w2 = w(:,2);

    % simulate
    parfor (ii = 1:n,poolsize)

        u(ii,:) = [w1(ii) hinv(w2(ii),w1(ii),family,theta)];

    end
    
else
    
    % simulate
    for ii = 1:1:n

        u(ii,:) = [w(ii,1) hinv(w(ii,2),w(ii,1),family,theta)];

    end
    
end

% shut down parallel pool
if para
    delete(ppool);
end


end