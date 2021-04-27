function [tstat, pval] = goftest_a2(u,A,family,theta,varargin)
% Conducts a goodness-of-fit test of type A2 from Berg(2009) for  a 
% specified simplified vine copula structure and given data. The test is 
% based on a Cramer-von Mises statistic of the difference between empirical 
% copula and estimated copula.
% 
% call: [tstat, pval] = goftest_a2(u,A,family,theta[,parallel])
%
% input     u                   - nxd data matrix of pseudo-observations
%           A                   - a vine array; note that a feasible  
%                                 structure has to be used, since the 
%                                 function does not check this 
%           family              - a (d-1)x(d-1) cell variable determining 
%                                 the copula families; possible families: 
%                                 'gumbel', 'clayton', 'frank', 't', 
%                                 'gauss', 'ind', 'amhaq', 'tawn', 'fgm',  
%                                 'plackett', 'joe', 'surclayton', 
%                                 'surgumbel', 'surjoe'
%           theta               - a (d-1)x(d-1) cell variable of copula 
%                                 parameters; for t-copula insert [rho nu]  
%                                 in cell element
%           parallel (optional) - switch parallelization on (=1) 
%                                 or off (=0); default: 0
%
% output    tstat           - the test statistic value
%           pval            - the corresponding p-value
%
%
% How does it work?
% The function conducts a goodness-of-fit test of an estimated vine
% structure on a given data set (estimation of the vine has to be done in a
% previous step, see function ssp). Under the null hypothesis the estimated 
% copula model has appropriate fit.
%
% Note that for the function to  work, the vine array provided by the user 
% has to be a feasible vine array in the first place.  The function will 
% not check feasibilty on its own! For c- and d-vines the function 
% cdvinearray can be used to generate a feasible vine array.
% 
% Structure of the input is demonstrated for a 5-dimensional r-vine copula:
% 
% Let the sample r-vine structure be
%
%           4
%          /
% 1 - 2 - 3 
%          \
%           5
%
% 12 - 23 - 34 - 35
%       .
%       .
%       .
%
% , where the numbers correspond to the columns of the matrix u. In this 
% case 
%
%     1 1 2 3 3
%       2 1 2 4
% A =     3 1 2
%           4 1
%             5
%
% is the corresponding vine array.
%
% In order for the function to work, the user has to input information on
% the following bivariate copulas: 12, 23, 34, 35, 13|2, 24|3, 45|3, 14|23,
% 25|34, 15|234, where '|' represents conditioning. Note that this system 
% corresponds to the appearance of the copula in the vine array from left 
% to right. Input family and theta cell variables for the copulas like 
% this:
%
%          family12     family23    family34   family35
% family = family13|2   family24|3  family45|3    0
%          family14|23  family25|34    0          0
%          family15|234    0           0          0
%
% Matlab syntax: 
%    family = {'family12','family23','family34','family35'; 'family13|2','family24|3','family45|3',0; 'family14|23','family25|34',0,0;'familiy15|234',0,0,0}
%
%          theta12     theta23    theta34   theta35
% theta =  theta13|2   theta24|3  theta45|3    0
%          theta14|23  theta25|34    0         0
%          theta15|234    0          0         0
%
% Matlab syntax: 
%    theta = {theta12,theta23,theta34,theta35; theta13|2,theta24|3,theta45|3,0; theta14|23,theta25|34,0,0;theta15|234,0,0,0}
% 
% Additionally, input the vine array A as stated above.
%
% 
% References:
% Berg (2009), Copula Goodness-of-fit Testing: An Overview and Power
% Comparison, The European Journal of Finance, Vol. 15(7-8), 675-701.
% Genest et al (2009), Goodness-of-fit Tests for Copulas: A Review and
% Power Study, Insurance: Mathematics and Economics, Vol. 44, 199-213.
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

% some parsing
p = inputParser;
p.addRequired('u',@ismatrix);
p.addRequired('A',@ismatrix);
p.addRequired('family',@iscell);
p.addRequired('theta',@iscell);
p.addOptional('para',0,@isscalar);
p.parse(u,A,family,theta,varargin{:});

% sanity checks
for ii = 1:1:size(family,1)
    
    for jj = 1:1:size(family,2)-ii+1
        
        if ~cpcheck(family{ii,jj},theta{ii,jj})
            error('goftest_a2:InvalidParameter',['invalid parameter for ',family{ii,jj},' copula at (',num2str(ii),',',num2str(jj),')']);    
        end
        
    end % jj
    
end % ii

% some sanity checks for vine array A
if (size(A,1) ~= size(A,2))
    error('goftest_a2:InvalidVineArray','vine array A has to be a quadratic matrix');
end

for jj = 1:1:size(A,1)
    
    if (length(unique(A(1:jj,jj))) ~= jj)
        error('goftest_a2:InvalidVineArray','input A is not a vine array');
    end
    
end % jj

% determine parallelization mode
para = 0;
if nargin > 4
    if varargin{1} ~= 1 && varargin{1} ~= 0
        error('goftest_a2:InvalidParallelizationMode','optional input argument parallel has to be 0 or 1')
    end
    para = varargin{1};
end

% initialize variables
n = size(u,1);

k = 200; % number of bootstraps

if 2*n < 5000 % number of resampled points from vine structure
    n_b = 5000; 
elseif 2*n >= 1000000 
    n_b = n;
else
    n_b = 2*n;
end

tstatbs = zeros(k,1);

% start testing procedure
% 1 - empirical copula
c_hat = ecopula(u); 

% 2 - analytical expression unknown, therefore approximate copula CDF
x_star = simrvine(n_b,A,family,theta);
c_thetahat = ecopula(u,x_star);

% 3 - Cramer-von Mises test statistic
tstat = sum((c_hat - c_thetahat).^2);

% 4 - approximate corresponding p value
if para
    
    ppool = gcp; 
    poolsize = ppool.NumWorkers;
    
    parfor (ii = 1:k,poolsize)

        % draw bootstrap sample   
        x_bs = simrvine(n,A,family,theta);
        [thetabs,~,~,~] = ssp(x_bs,A,family);
        c_hatbs = ecopula(x_bs);

        % approximate CDF
        x_starbs = simrvine(n_b,A,family,thetabs);
        c_thetahatbs = ecopula(x_bs,x_starbs);

        % save bootstrapped test statistic
        tstatbs(ii) = sum((c_hatbs - c_thetahatbs).^2);

    end % ii
    
else
    
    for ii = 1:1:k

        % draw bootstrap sample   
        x_bs = simrvine(n,A,family,theta);
        [thetabs,~,~,~] = ssp(x_bs,A,family);
        c_hatbs = ecopula(x_bs);

        % approximate CDF
        x_starbs = simrvine(n_b,A,family,thetabs);
        c_thetahatbs = ecopula(x_bs,x_starbs);

        % save bootstrapped test statistic
        tstatbs(ii) = sum((c_hatbs - c_thetahatbs).^2);

    end % ii
    
end

pval = sum(tstatbs > tstat)/(k+1);

% shut down parallel pool
if para
    delete(ppool);
end


end