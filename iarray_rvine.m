function I = iarray_rvine(A)
% Helper function to speed up algorithms for d-dimensional simplified 
% r-vine copulas.
%
% call: I = iarray_rvine(A)
%
% input     A       - a vine array; note that a feasible structure has to
%                     be used, since the function will not check this
%
% output    I       - indicator matrix needed for r-vine copula algorithms 
%                     (e.g. simulation) 
%
%
% References: 
% Joe (2015), Dependence Modeling with Copulas, CRC Press.
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

% some sanity checks for vine array A
if (size(A,1) ~= size(A,2))
    error('iarray_rvine:InvalidVineArray','vine array A has to be a quadratic matrix');
end

for jj = 1:1:size(A,1)
    
    if (length(unique(A(1:jj,jj))) ~= jj)
        error('iarray_rvine:InvalidVineArray','input A is not a vine array');
    end
    
end % jj

% initialze variables
d = size(A,2);
M = zeros(d);
I = zeros(d);

% compute matrix M
for jj = 2:1:d
    
    for kk = 1:1:jj-1
        
        M(kk,jj) = max(A(1:kk,jj));
        
    end
    
end % jj

% now compute indicator matrix I
for kk = 2:1:d-1
    
    for jj = kk+1:1:d
        
        if A(kk,jj) < M(kk,jj)
            
            I(kk-1,M(kk,jj)) = 1;
            
        end
                    
    end % jj
    
end % kk


end