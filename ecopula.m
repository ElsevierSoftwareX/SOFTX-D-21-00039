function ret = ecopula(u,varargin)
% Calculates the empirical copula values for multivariate data u based on 
% data v.
%
% call: ret = ecopula(u[,v])
%
% input     u                   - nxp data matrix with n observations over 
%                                 p dimensions, for which the empirical  
%                                 copula value is calculated
%           v (optional)        - mxp data matrix, on which the empirical 
%                                 copula is based
%
% output    ret                 - vector of the empirical copula values for 
%                                 data u
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

% some parsing
p = inputParser;
p.addRequired('X',@ismatrix);
p.addOptional('Y',0,@ismatrix);
p.parse(u,varargin{:});

% initialization of variables
[n,p] = size(u);

if nargin == 2 % in case u is evaluated on empirical copula based on v
    v = varargin{1};
    m = size(v,1);
else % else u is evaluated on the empirical copula of itself
    v = u;
    m = n;
end

Fn = zeros(n,1);

% compute empirical copula    
for ii = 1:1:n
    
    Fn(ii) = sum(all(reshape(v<=repmat(u(ii,:),m,1),m,p),2))/(m+1);
    
end % ii
    
ret = Fn;


end