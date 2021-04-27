function [pval, tstat] = vuongtest(loglik1,loglik2,nop1,nop2)
% Conducts a Vuong's test (1989) for non-nested models.
%
% call: [pval, tstat] = vuongtest(loglik1,loglik2,nop1,nop2)
%
% input     loglik1 - column vector of loglikelihoods of model 1
%           loglik2 - column vector of loglikelihoods of model 2
%           nop1    - number of parameters of model 1
%           nop2    - number of parameters of model 2
%
% output    pval    - p value according to a standard normal distribution
%           tstat   - the test statistic
%
% 
% How does it work:
% The Vuong test is a form of likelihood ratio test for model
% comparison/selection. This function implements the strictly non-nested
% version of the test. The output can be interpreted as follows: If pval is
% below your assumed error rate alpha, the test is significant. A positive
% tstat indicates that model 1 is superior compared to model 2. Vice versa,
% a negative tstat indicates that model 2 is superior compared to model 1.
%
%
% References:
% Vuong (1989), Likelihood Ratio Tests for Model Selection and Non-Nested
% Hypotheses, Econometrica, Vol. 57(2), 307-333.
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

% some parsing
p = inputParser;
p.addRequired('loglik1',@isvector);
p.addRequired('loglik2',@isvector);
p.addRequired('nop1',@isscalar);
p.addRequired('nop2',@isscalar);
p.parse(loglik1,loglik2,nop1,nop2);

% sanity check
if size(loglik1,1) ~= size(loglik2,1)
    error('number of loglikelihood values for model 1 and for model 2 have to agree');
end

% initialize parameters
n = size(loglik1,1);

% compute test statistic
omega = std(loglik1-loglik2);
lr = sum(loglik1)-sum(loglik2)-(nop1-nop2)/2*log(n);
tstat = lr/(sqrt(n)*omega);

% compute p-value
if tstat < 0
    pval = normcdf(tstat);
else
    pval = 1-normcdf(tstat);
end


end