function [thetahat,loglik,sll,fam] = ssp(u,A,family)
% Implements the stepwise semiparametric estimator for simplified vine 
% copulas.
%
% call: [thetahat,loglik,sll,fam] = ssp(u,A,family)
%
% input     u               - nxd data matrix of pseudo-observations
%           A               - a vine array; note that a feasible structure 
%                             has to be used, since the function will not 
%                             check this 
%           family          - a (d-1)x(d-1) cell variable determining the 
%                             copula families used in the estimation 
%                             process; possible families: 
%                             'gumbel', 'clayton', 'frank', 't', 'gauss', 
%                             'ind', 'amhaq', 'tawn', 'fgm', 'plackett', 
%                             'joe', 'surclayton', 'surgumbel', 'surjoe';
%                             if instead of a valid family 'aic', 'bic', or
%                             'sll' is chosen, the program will select the
%                             best copula family on its own based on the
%                             chosen criterion;
%                             alternatively, a 1x1 cell variable with 
%                             'aic', 'bic', or 'sll' can be provided if all 
%                             copulas should be estimated according to the 
%                             chosen criterion                          
%
% output    thetahat        - (d-1)x(d-1) cell variable of estimated  
%                             parameters in the given vine copula structure
%           loglik          - column vector of loglikelihoods for each data
%                             point
%           sll             - the sum of logliklihoods 
%           fam             - a (d-1)x(d-1) cell variable indicating the
%                             copula families in the vine structure; this 
%                             will be different to the input 'family' only
%                             if some copula families are not prespecified
%                             by the user but selected by the program
%
%
% How does it work?
% This function performs a semiparametric stepwise estimation for
% simplified vine copulas. Depending on the input family, the function 
% chooses either one copula among different families of copulas (the one 
% yielding the best selection criterion), estimates the parameter of one 
% copula prespecified by the user, or selects the best fitting copula for a 
% given set of copulas. If more than one family is prespecified the 
% function selects the best fitting copula among the prespecified set.
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
% The outputs thetahat and fam will have the same structure as the input 
% family.
%
%
% References: 
% Diﬂmann et al (2013), Selecting and Estimating Regular Vine Copulae and 
% Applications to Financial Returns, Computational Statistics and Data 
% Analysis, Vol. 59, 52-69.
% Hobaek Haff (2012), Parameter Estimation for Pair-Copula Constructions,
% Bernoulli, Vol. 19(2), 462-491.
% Joe (2015), Dependence Modeling with Copulas, CRC Press.
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

% some parsing
p = inputParser;
p.addRequired('x',@ismatrix);
p.addRequired('A',@ismatrix);
p.addRequired('family',@iscell);
p.parse(u,A,family);

% sanity checks
if (size(A,1) ~= size(A,2))
    error('ssp:InvalidVineArray','vine array A has to be a quadratic matrix');
end

for jj = 1:1:size(A,1)
    
    if (length(unique(A(1:jj,jj))) ~= jj)
        error('ssp:InvalidVineArray','input A is not a vine array');
    end
    
end % jj

% Initialize variables
n = size(u,1);
d = size(u,2);

thetahat = cell(d-1,d-1);
thetahat(:) = {0};
loglik = zeros(n,1);
fam = cell(d-1,d-1);
fam(:) = {0};

% check, whether input family is {'aic'}, {'bic'}, or {'sll'}
if sum(size(family) == [1 1]) == 2
    
    straux = family{1,1};
    family = cell(d-1,d-1);
    family(:) = {straux};
    
end

% start estimation

% permute A, such that a_jj = jj
[A,perm,~] = transforma(A);
u = u(:,perm(:,1));

M = zeros(d);
% compute matrix M
for jj = 2:1:d
    
    for kk = 1:1:jj-1
        
        M(kk,jj) = max(A(1:kk,jj));
        
    end % kk
    
end % jj

v = cell(d,d);
v_prime = cell(d,d);
for jj = 1:1:d
    v{1,jj} = u(:,jj);
    v_prime{1,jj} = u(:,jj);
end % jj

% levels 1 to d
for kk = 2:1:d
    
    for ii = 1:1:kk-1
        
        % select correct variables
        z1 = v{ii,kk};
        if M(ii,kk) == A(ii,kk)
            z2 = v{ii,M(ii,kk)};
        else
            z2 = v_prime{ii,M(ii,kk)};
        end
        
        % estimate copula
        if ~iscell(family{ii,kk-ii})
            
            [fam{ii,kk-ii}, thetahat{ii,kk-ii}, loglik_h, ~] = copulaselect([z2 z1],family{ii,kk-ii});
            
        else
            
            famaux = ['copulaselect([z2 z1],''',family{ii,kk-ii}{1}];
            for jj = 2:1:size(family{ii,kk-ii},2)
                
                famaux = strcat(famaux,''',''',family{ii,kk-ii}{jj});
                
            end % kk
            famaux = strcat(famaux,''')');
            [fam{ii,kk-ii}, thetahat{ii,kk-ii}, loglik_h, ~] = eval(famaux);
            
        end
        
        loglik = loglik + loglik_h;
        
        % compute pseudo-observations
        v{ii+1,kk} = hfunc(z1,z2,fam{ii,kk-ii},thetahat{ii,kk-ii});
        v_prime{ii+1,kk}  = hfunc(z2,z1,fam{ii,kk-ii},thetahat{ii,kk-ii});
        
    end % ii
    
end % kk

% compute sum of loglikelihoods
sll = sum(loglik);


end