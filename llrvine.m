function [loglik,sll] = llrvine(u,A,family,theta)
% This function evaluates the loglikelihood function of a simplified
% r-vine copula.
%
% call: [loglik,sll] = llrvine(u,A,family,theta)
%
% input     u       - nxd data matrix of pseudo-observations for which the 
%                     loglikelihood function is evaluated
%           A       - a vine array; note that a feasible structure has to
%                     be used, since the function will not check this
%           family  - a (d-1)x(d-1) cell variable determining the copula 
%                     families used in the r-vine-structure
%                     possible families: 'gumbel', 'clayton', 'frank', 't',
%                                        'gauss', 'ind', 'amhaq', 'tawn',
%                                        'fgm', 'plackett', 'joe', 
%                                        'surclayton', 'surgumbel', 
%                                        'surjoe'
%           theta   - a (d-1)x(d-1) cell variable of copula parameters;
%                     for t-copula insert [rho nu] in cell element
%
% output    loglik  - nx1 vector of loglikelihoods of the r-vine copula 
%                     applied to each data point in u
%           sll     - sum of loglikelihoods
%
%
% How does it work?
% The function evaluates the loglikelihood function of a d-dimensional 
% simplified r-vine of nxd data matrix u. Note that for the function to 
% work, the vine array provided by the user has to be a feasible vine array 
% in the first place. The function will not check feasibilty on its own! 
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
% , where the numbers correspond to the columns of input matrix u. In this 
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
%
% References: 
% Joe (2015), Dependence Modeling with Copulas, CRC Press.
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

% some parsing
p = inputParser;
p.addRequired('x',@ismatrix);
p.addRequired('vineArray',@ismatrix);
p.addRequired('family',@iscell);
p.addRequired('theta',@iscell);
p.parse(u,A,family,theta);

% some sanity checks
for ii = 1:1:size(family,1)
    
    for jj = 1:1:size(family,2)-ii+1
        
        if ~cpcheck(family{ii,jj},theta{ii,jj})
            error('llrvine:InvalidParameter',['invalid parameter for ',family{ii,jj},' copula at (',num2str(ii),',',num2str(jj),')']);    
        end
        
    end % jj
    
end % ii

if (size(A,1) ~= size(A,2))
    error('llrvine:InvalidVineArray','vine array A has to be a quadratic matrix');
end

for jj = 1:1:size(A,1)

    if (length(unique(A(1:jj,jj))) ~= jj)
        error('llrvine:InvalidVineArray','input A is not a vine array');
    end
    
end % jj

%initialize variables
n = size(u,1);
d = size(u,2);
M = zeros(d);
loglik = zeros(n,1);
v = zeros(n,d);
v_prime = zeros(n,d);
w = zeros(n,d);
w_prime = zeros(n,d);
s = zeros(n,d);

% permute A, such that a_jj = jj
[A,perm,~] = transforma(A);
u = u(:,perm(:,1));

% compute matrix M
for jj = 2:1:d
    
    for kk = 1:1:jj-1
        
        M(kk,jj) = max(A(1:kk,jj));
        
    end % kk
    
end % jj

% compute indicator matrix I
I = iarray_rvine(A);

% first tree (T1)
for ii = 2:1:d
    
    loglik = loglik + log(copulapdfadv(family{1,ii-1},[u(:,A(1,ii)) u(:,ii)],theta{1,ii-1}));
    
end % ii

% trees 2-d (T2, T3,...,Td)
for jj = 2:1:d
    
    if I(1,jj) == 1
        v_prime(:,jj) = hfunc(u(:,A(1,jj)),u(:,jj),family{1,jj-1},theta{1,jj-1});
    end
    
    v(:,jj) = hfunc(u(:,jj),u(:,A(1,jj)),family{1,jj-1},theta{1,jj-1});
    
end % jj

for jj=3:1:d
    
    if A(2,jj) == M(2,jj)
        s(:,jj) = v(:,M(2,jj));
    else
        s(:,jj) = v_prime(:,M(2,jj));
    end
    
    loglik = loglik + log(copulapdfadv(family{2,jj-2},[s(:,jj) v(:,jj)],theta{2,jj-2}));
    
end % jj

w(:,3:d) = v(:,3:d);
w_prime(:,3:d) = v_prime(:,3:d);

for ll = 3:1:d-1
    
    for jj = ll:1:d
        
        if I(ll-1,jj) == 1
            v_prime(:,jj) = hfunc(s(:,jj),w(:,jj),family{ll-1,jj-ll+1},theta{ll-1,jj-ll+1});
        end
        
        v(:,jj) = hfunc(w(:,jj),s(:,jj),family{ll-1,jj-ll+1},theta{ll-1,jj-ll+1});
        
    end % jj
    
    for jj = ll+1:1:d
        
        if A(ll,jj) == M(ll,jj)
            s(:,jj) = v(:,M(ll,jj));
        else
            s(:,jj) = v_prime(:,M(ll,jj));
        end
        
        loglik = loglik + log(copulapdfadv(family{ll,jj-ll},[s(:,jj) v(:,jj)],theta{ll,jj-ll}));
        
    end % jj
    
    w(:,ll+1:d) = v(:,ll+1:d);
    w_prime(:,ll+1:d) = v_prime(:,ll+1:d);
    
end % ll

sll = sum(loglik);


end