function value = vinecdf(u,A,family,theta)
% Computes the cdf of a simplified vine copula via Monte-Carlo integration.
%
% call: value = vinecdf(u,A,family,theta)
%
% input     u               - nxd data matrix of pseudo-observations
%           A               - a vine array; note that a feasible structure 
%                             has to be used, since the function will not 
%                             check this
%           family          - a (d-1)x(d-1) cell variable determining the
%                             copula families used in the calculation; 
%                             possible families:  'gumbel', 'clayton', 
%                             'frank', 't', 'gauss', 'ind', 'amhaq', 
%                             'tawn', 'fgm', 'plackett', 'joe', 
%                             'surclayton', 'surgumbel', 'surjoe'
%           theta           - a (d-1)x(d-1) cell variable of copula
%                             parameters; for t-copula insert [rho nu] in
%                             cell element
%           
%
% output    value           - column vector of cdf values for each point in 
%                             u
%
%
% How does it work?
% This function computes the cdf values of points from a given simplified
% vine copula.
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
% , where the numbers correspond to the columns of the output. In this case
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
% to right. Input family cell variable for the copulas like this:
%
%          family12     family23    family34   family35
% family = family13|2   family24|3  family45|3    0
%          family14|23  family25|34    0          0
%          family15|234    0           0          0
%
% Matlab syntax:
%    family = {'family12','family23','family34','family35'; 'family13|2','family24|3','family45|3',0; 'family14|23','family25|34',0,0;'familiy15|234',0,0,0}
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
p.parse(u,A,family,theta);

% sanity checks
for ii = 1:1:size(family,1)
    
    for jj = 1:1:size(family,2)-ii+1
        
        if ~cpcheck(family{ii,jj},theta{ii,jj})
            error(['invalid parameter for ',family{ii,jj},' copula at (',num2str(ii),',',num2str(jj),')']);    
        end
        
    end % jj
    
end % ii

% some sanity checks for vine array A    
if (size(A,1) ~= size(A,2))
    error('vine array A has to be a quadratic matrix');
end

for jj = 1:1:size(A,1)
    
    if (length(unique(A(1:jj,jj))) ~= jj)
        error('input A is not a vine array');
    end
    
end % jj


% Initialize variables
n = size(u,1);
d = size(u,2);
n_mc = 1000000;
value = ones(n,1);

% start calculation
for ii = 1:1:n
    
    u = rand(n_mc,d).*u(ii,:);
    value(ii) = prod(u(ii,:))/n_mc*sum(vinepdf(u,A,family,theta));
   
end % ii


end