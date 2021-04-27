function [fam, thetahat, loglik, mcrit] = copulaselect(u, varargin)
% Performs a bivariate copula selection on data u based on Maximum 
% Likelihood Estimation. 
%
% call: [fam, thetahat, loglik, mcrit] = copulaselect(u[,crit,family1,family2,...])
%
% input     u                               - nx2 matrix of
%                                             pseudo-observations
%           crit (optional)                 - the selection criterion: 
%                                             'aic' (Akaike's information 
%                                             criterion), 'bic' (Bayesian 
%                                             information criterion), 'sll' 
%                                             (sum of loglikelihoods); 
%                                             default is 'aic'
%           family1...familyk (optional)    - the copula family: 'gumbel', 
%                                             'clayton', 'frank', 't', 
%                                             'gauss', 'amhaq', 'ind',  
%                                             'tawn', 'fgm', 'plackett', 
%                                             'joe', 'surclayton', 
%                                             'surgumbel', 'surjoe'                                        
%
% output    fam                             - the estimated copula family: 
%                                             'gumbel', 'clayton', 'frank', 
%                                             't', 'gauss', 'amhaq', 'ind',  
%                                             'tawn', 'fgm', 'plackett', 
%                                             'joe', 'surclayton', 
%                                             'surgumbel', 'surjoe' 
%           thetahat                        - estimated copula parameters; 
%                                             for t-copula [rho, nu]
%           loglik                          - column vector of 
%                                             loglikelihoods for each data
%                                             point
%           mcrit                           - the selection criterion value 
%                                             for the chosen model; in case 
%                                             the optional input family is
%                                             specified by the user this is 
%                                             the sum of loglikelihoods
%
%
% How does it work:
% The function uses Maximum Likelihood Estimation for copula selection.
% Depending on the number of inputs, the function chooses either one copula
% among different families of copulas (the one yielding the best selection 
% criterion), estimates the parameter of one copula prespecified by the 
% user, or selects the best fitting copula for a given set of copulas. If 
% more than one family is prespecified the function selects the best
% fitting copula among the prespecified set.
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

% some parsing
p = inputParser;
p.addRequired('u',@ismatrix);
p.parse(u);

% sanity checks
if (nargin < 2) 
    error('copulaselect:InvalidNumberOfInputs','Number of input arguments has to be at least 2. Please specify a selection criterion or at least one copula family.');
end

tau = corr(u,'type','Kendall');
tau = tau(1,2);

if nargin < 3
    
    if (tau < 0 && sum(strcmpi(varargin{1},{'gumbel','clayton','tawn','joe','surclayton','surgumbel','surjoe'})))
        warning([varargin{1},' copula can only be estimated on positively dependent data']);
    end
    
else
    
    for ii = 1:1:nargin-1
        if(tau < 0 && sum(strcmpi(varargin{ii},{'gumbel','clayton','tawn','joe','surclayton','surgumbel','surjoe'})))
            warning([varargin{ii},' copula can only be estimated on positively dependent data']);
        end
    end % ii
    
end

% begin estimation
if nargin < 3
    
    if sum(strcmpi(varargin{1},{'aic','bic','sll'})) % input argument crit
                                                     % is specified,
                                                     % hence find best
                                                     % copula
        
        crit = varargin{1};
        
        % independence copula
        fam = 'ind';
        thetahat = 0;
        loglik = log(ones(size(u,1),1));
        mcrit = cscrit(crit,loglik,0);
        
        % all other copulas
        for ii = {'gauss','t','frank','gumbel','clayton','amhaq','tawn','fgm','plackett','joe','surclayton','surgumbel','surjoe'}
           
            [thetahat_aux, loglik_aux] = copmle(u,ii{1});
            
            if cscrit(crit,loglik_aux,nopcount(ii(1))) < mcrit
                fam = ii{1};
                thetahat = thetahat_aux;
                loglik = loglik_aux;
                mcrit = cscrit(crit,loglik,nopcount(ii(1)));
            end
            
        end % ii
            
        if strcmp(crit,'sll')
            mcrit = -1*mcrit;
        end
        
    else % input argument family1 is specified, hence find best fit for chosen copula
        
        fam = varargin{1};
        [thetahat, loglik] = copmle(u,fam);
        mcrit = sum(loglik);
        
    end
    
else % a set of copulas to select from is given, AIC is assumed
    
    fam = varargin{1};
    [thetahat, loglik] = copmle(u,fam);
    mcrit = cscrit('aic',loglik,nopcount(varargin(1)));
    
    for ii = 2:1:nargin-1
        
        [thetahat_aux, loglik_aux] = copmle(u,varargin{ii});
        
        if cscrit('aic',loglik_aux,nopcount(varargin(ii))) < mcrit
            fam = varargin{ii};
            thetahat = thetahat_aux;
            loglik = loglik_aux;
            mcrit = cscrit('aic',loglik,nopcount(varargin(ii)));
        end
        
    end % ii
    
end


end