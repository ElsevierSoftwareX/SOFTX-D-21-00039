function u = pobs(X)
% Rank-transforms observation matrix X to pseudo-observation matrix u. The 
% pseudo-observations can be used to work with copulas.
%
% call: u = pobs(X)
%
% input     X   - nxd data matrix, where n is the sample size and d is the 
%                 number of dimensions                             
%
% output    u   - the rank-transformed nxd matrix of pseudo-observations
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

% some parsing
p = inputParser;
p.addRequired('X',@ismatrix);
p.parse(X);

% initialization of variables
n = size(X,1); % sample size
d = size(X,2); % dimension
u = zeros(size(X));

% rank-transform data
for ii = 1:1:d
    
    u(:,ii) = tiedrank(X(:,ii))/(n+1);
    
end % ii


end