function A = cdvinearray(vine,d)
% Generates the vine array A for a standard c- or d-vine.
%
% call: A = cdvinearray(vine,d)
%
% input     vine    - the type of vine; options are 'c' or 'd'
%           d       - dimension
%
% output    A       - dxd vine array; note that the diagonal elements a_jj
%                     fulfill a_jj = j
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

% some parsing
p = inputParser;
p.addRequired('vine',@isstr);
p.addRequired('d',@isscalar);
p.parse(vine,d);

% sanity checks
if (sum(strcmpi(vine,{'c','d'})) ~= 1)
    error('cdvinearray:WrongVineType','only vine arrays for c- or d-vine can be generated');
end

% initialize variables
A = diag(1:d);

% generate vine array
switch lower(vine)
    
    case 'c'
        for jj = 1:1:d-1
            
            A = A + diag(1:d-jj,jj);
            
        end % jj
        
    case 'd'
        for jj = 1:1:d-1
            
            A = A + diag(jj*ones(1,d-jj),jj);
            
        end % jj
        
end % switch


end