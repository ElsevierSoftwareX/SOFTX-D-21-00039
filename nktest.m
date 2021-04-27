function [tstat, pval] = nktest(u,parallel,family1,theta1,A1,family2,theta2,A2,varargin)
% Conducts a model comparison test of k models according to 
% Nikoloulopoulos & Karlis (2008) adapted for vine copulas.
% 
% call: [tstat, pval] = nktest(u,parallel,family1,theta1,A1,family2,theta2,A2[,family3,theta3,A3,...])
%
% input     u                       - nxd data matrix of
%                                     pseudo-observations
%           parallel                - switch parallelization on (=1) or 
%                                     off (=0)
%           family1                 - a (d-1)x(d-1) cell variable  
%                                     determining the copula families for 
%                                     model1; possible families: 'gumbel',
%                                     'clayton', 'frank', 't', 'gauss', 
%                                     'ind', 'amhaq', 'tawn', 'fgm', 'joe',
%                                     'plackett', 'surclayton', 'surjoe',
%                                     'surgumbel' 
%           theta1                  - a (d-1)x(d-1) cell variable of copula 
%                                     parameters for model1; for t-copula  
%                                     insert [rho nu] in cell element
%           A1                      - vine array for mode11; note that a 
%                                     feasible structure has to be used, 
%                                     since the function will not check 
%                                     this
%           family2                 - a (d-1)x(d-1) cell variable  
%                                     determining the copula families for 
%                                     model2; possible families: 'gumbel',
%                                     'clayton', 'frank', 't', 'gauss', 
%                                     'ind', 'amhaq', 'tawn', 'fgm', 'joe',
%                                     'plackett', 'surclayton', 'surjoe',
%                                     'surgumbel' 
%           theta2                  - a (d-1)x(d-1) cell variable of copula 
%                                     parameters for model2; for t-copula  
%                                     insert [rho nu] in cell element
%           A2                      - vine array for mode12; note that a 
%                                     feasible structure has to be used, 
%                                     since the function will not check 
%                                     this
%           (optional) if more than 2 models are compared:
%           family3,...,familiyk    - a (d-1)x(d-1) cell variable  
%                                     determining the copula families;  
%                                     possible families: 'gumbel',  
%                                     'clayton', 'frank', 't', 'gauss', 
%                                     'ind', 'amhaq', 'tawn', 'fgm', 'joe',
%                                     'plackett', 'surclayton', 'surjoe',
%                                     'surgumbel' 
%           theta3,...,thetak       - a (d-1)x(d-1) cell variable of copula 
%                                     parameters; for t-copula insert  
%                                     [rho nu] in cell element
%           A3,...,Ak               - a vine array; note that a feasible 
%                                     structure has to be used, since the 
%                                     function will not check this
%
% output    tstat                   - kx1 vector of test statistic values
%           pval                    - kx1 vector of corresponding p values
%
% 
% How does it work?
% The function conducts a model comparison test according to 
% Nikoloulopoulos & Karlis (2008) adapted for vine copulas. Note that the 
% models have to be fitted in advance. The minimum number of models 
% compared is two. The function can compare k models at once. The null
% hypothesis is that model k is correct. For each model beyond the first 
% two, specify the family, theta and A variables analoguesly to the 
% first two. 
% 
% Note that for the function to  work, the vine arrays provided by the user 
% have to be a feasible vine arrays in the first place.  The function will 
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
% Additionally, input the vine array A as a matrix as stated above.
%
%
% References:
% Allcroft & Glasbey (2003), A Simulation-Based Method for Model
% Evaluation, Statistical Modelling, Vol. 3, 1-13.
% Nikoloulopoulos & Karlis (2008), Copula Model Evaluation Based on
% Parametric Bootstrap, Computational Statistics and Data Analysis, Vol.
% 52, 3342-3353.
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

% some parsing
p = inputParser;
p.addRequired('x',@ismatrix);
p.addRequired('para',@isscalar);
p.addRequired('family1',@iscell);
p.addRequired('theta1',@iscell);
p.addRequired('A1',@ismatrix);
p.addRequired('family2',@iscell);
p.addRequired('theta2',@iscell);
p.addRequired('A2',@ismatrix);
p.parse(u,parallel,family1,theta1,A1,family2,theta2,A2);

% sanity checks
if mod(nargin-2,3) ~= 0
    error('nktest:InvalidNumberOfArguments','please specify the correct number of optional input arguments (3 per model)');
end

for ii = 1:1:size(family1,1)
    
    for jj = 1:1:size(family1,2)-ii+1
        
        if ~cpcheck(family1{ii,jj},theta1{ii,jj})
            error('nktest:InvalidParameter',['invalid parameter for ',family1{ii,jj},' copula at (',num2str(ii),',',num2str(jj),') in family1']);    
        end
        
    end % jj
    
end % ii

for ii = 1:1:size(family2,1)
    
    for jj = 1:1:size(family2,2)-ii+1
        
        if ~cpcheck(family2{ii,jj},theta2{ii,jj})
            error('nktest:InvalidParameter',['invalid parameter for ',family2{ii,jj},' copula at (',num2str(ii),',',num2str(jj),') in family2']);    
        end
        
    end % jj
    
end % ii

% some sanity checks for vine array A1
if (size(A1,1) ~= size(A1,2))
    error('nktest:InvalidVineArray','vine array A1 has to be a quadratic matrix');
end

for jj = 1:1:size(A1,1)
    
    if (length(unique(A1(1:jj,jj))) ~= jj)
        error('nktest:InvalidVineArray','input A1 is not a vine array');
    end
    
end % jj

% some sanity checks for vine array A2
if (size(A2,1) ~= size(A2,2))
    error('nktest:InvalidVineArray','vine array A2 has to be a quadratic matrix');
end

for jj = 1:1:size(A2,1)
    
    if (length(unique(A2(1:jj,jj))) ~= jj)
        error('nktest:InvalidVineArray','input A2 is not a vine array');
    end
    
end % jj

if nargin > 8
    
    m = 3;
    for ii = 1:3:nargin-8
        
        for kk = 1:1:size(varargin{ii},1)
            
            for jj = 1:1:size(varargin{ii},2)-kk+1
                
                if ~cpcheck(varargin{ii}{kk,jj},varargin{ii+1}{kk,jj})
                    error('nktest:InvalidParameter',['invalid parameter for ',varargin{ii}{kk,jj},' copula at (',num2str(kk),',',num2str(jj),') for family',num2str(m)]);
                end
                
            end % jj
            
        end % kk
        
        % some sanity checks for vine arrays A3...Ak
        aux = varargin{ii+2};
        
        if (size(aux,1) ~= size(aux,2))
            error('nktest:InvalidVineArray',['vine array A', num2str(m), ' has to be a quadratic matrix']);
        end
        
        for jj = 1:1:size(aux,1)
            
            if (length(unique(aux(1:jj,jj))) ~= jj)
                error('nktest:InvalidVineArray',['input A', num2str(m), ' is not a vine array']);
            end
            
        end % jj
        
        m = m+1;
        
    end % ii
    
end

% parallelization mode
if parallel ~= 1 && parallel ~= 0
        error('nktest:InvalidParallelizationMode','input argument parallel has to be 0 or 1')
end
if parallel
    ppool = gcp; 
    poolsize = ppool.NumWorkers;
end

% initialize some variables
n = size(u,1);
k = (nargin-2)/3;
loglik = zeros(k,1);
tstat = zeros(k,1);

B = 100; % number of bootstraps
bsample = cell(k,B);

lambda = zeros(k,B,k);

% 1 - compute log-likelihood of data x for each model 
% model 1
[~,loglik(1)] = llrvine(u,A1,family1,theta1);

% model 2
[~,loglik(2)] = llrvine(u,A2,family2,theta2);

if k > 2 % models 3 to k
    
    m = 1;
    
    for ii = 3:1:k

        [~,loglik(ii)] = llrvine(u,varargin{2+m},varargin{m},varargin{1+m});

        m = m + 3;
        
    end % ii
    
end

% 2 - Simulate B samples from each fitted model
% At the end of this step there is a series of B samples for each of the k
% models:
% model1: sample1 sample2 ... sampleB
% model2: sample1 sample2 ... sampleB
% ...
% modelk: sample1 sample2 ... sampleB
if parallel
    
    parfor (jj = 1:B,poolsize) % model 1

        bsample{1,jj} = simrvine(n,A1,family1,theta1);

    end % jj
    
    parfor (jj = 1:B,poolsize) % model 2

        bsample{2,jj} = simrvine(n,A2,family2,theta2);

    end % jj
    
else
    
    for jj = 1:1:B % models 1 & 2

        bsample{1,jj} = simrvine(n,A1,family1,theta1);
        bsample{2,jj} = simrvine(n,A2,family2,theta2);

    end % jj
    
end

if k > 2 % models 3 to k
    
    m = 1;
    for ii = 3:1:k
        
        if parallel
            
            parfor (jj = 1:B,poolsize)

                bsample{ii,jj} = simrvine(n,varargin{2+m},varargin{m},varargin{1+m});

            end % jj
            
        else
            
            for jj = 1:1:B

                bsample{ii,jj} = simrvine(n,varargin{2+m},varargin{m},varargin{1+m});

            end % jj
            
        end
        
        m = m + 3;
        
    end % ii
    
end

% 3 - For each of the series of B samples 1) either refit all the k 
% competing models, 2) or evaluate the log-likelihood function of the k 
% models at the sample. Here the models are refitted and the loglikelihoods
% evaluated.
% At the end of this step there is a kxB matrix of
% log-likelihoods for each of the k series.
if parallel
    
    parfor (jj = 1:B,poolsize) % model 1 
        
        lambda_aux = zeros(k,1);
        for ll = 1:1:k

            [~,~,lambda_aux(ll),~] = ssp(bsample{ll,jj},A1,family1);

        end % ll
        
        lambda(:,jj,1) = lambda_aux;

    end % jj
    
else
    
    for jj = 1:1:B % model 1 

        for ll = 1:1:k

            [~,~,lambda(ll,jj,1),~] = ssp(bsample{ll,jj},A1,family1);

        end % ll

    end % jj
    
end

if parallel
    
    parfor (jj = 1:B,poolsize) % model 2 
        
        lambda_aux = zeros(k,1);
        for ll = 1:1:k

            [~,~,lambda_aux(ll),~] = ssp(bsample{ll,jj},A2,family2);

        end % ll
        
        lambda(:,jj,2) = lambda_aux;

    end % jj
    
else
    
    for jj = 1:1:B % model 2 

        for ll = 1:1:k

            [~,~,lambda(ll,jj,2),~] = ssp(bsample{ll,jj},A2,family2);

        end % ll

    end % jj
    
end

if k > 2 % models 3 to k
    
    m = 1;
    for ii = 3:1:k
        
        if parallel
            
            parfor (jj = 1:B,poolsize)
                
                lambda_aux = zeros(k,1);
                for ll = 1:1:k

                    [~,~,lambda_aux(ll),~] = ssp(bsample{ll,jj},varargin{m+2},varargin{m});

                end % ll
                
                lambda(:,jj,ii) = lambda_aux;

            end % jj
            
        else
            
            for jj = 1:1:B

                for ll = 1:1:k

                    [~,~,lambda(ll,jj,ii),~] = ssp(bsample{ll,jj},varargin{m+2},varargin{m});

                end % ll

            end % jj
            
        end
        
        m = m + 3;
        
    end % ii
    
end

% 4 - test statistics and p values
if parallel
    
    parfor (ii = 1:k,poolsize)

        Sinv = inv(cov(lambda(:,:,ii)'));

        tstat(ii) = (loglik-mean(lambda(:,:,ii),2))'*Sinv*(loglik-mean(lambda(:,:,ii),2))/k;

    end % ii
    
else
    
    for ii = 1:1:k

        Sinv = inv(cov(lambda(:,:,ii)'));

        tstat(ii) = (loglik-mean(lambda(:,:,ii),2))'*Sinv*(loglik-mean(lambda(:,:,ii),2))/k;

    end % ii
    
end

pval = fcdf(tstat,k,B-1);


end