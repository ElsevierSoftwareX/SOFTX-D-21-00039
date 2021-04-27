function [At,p,pinv] = transforma(A)
% Transforms an arbitrary feasible vine array A to an ordered vine array, 
% i.e., such that the diagonal entries a_jj of A satisfy a_jj = j.
%
% call: [At,p,pinv] = transforma(A)
%
% input     A       - feasible vine array
%
% output    At      - ordered vine array, diagonal entries satisfy
%                     at_jj = j
%           p       - the permutation of variables that yields the
%                     transformed vine array
%           pinv    - the inverse of permutation p
%
% How does it work?
% The function transforms an arbitrary vine array A to an ordered vine 
% array which satisfies a_jj = j on its diagonal entries a_jj. For this 
% the entries in A are permuted. This is equivalent to a 
% renumbering/relabeling of the vine structure. The permutation for the 
% renumbering is given in the output p, which contains the unsorted 
% diagonal elements of A in the first column and the new labels in the 
% second column.
%
% In order to transform the columns of data u to reflect the renumbering,
% one can use 
% ut = u(:,p(:,1));
%
% Compliance between u and the diagonal elements of A is, e.g., needed for
% llrvine when comparing models.
%
% To reverse the permutation, i.e., to recover the original u, one can use
% u = ut(:,pinv(:,2));
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

% some parsing
p = inputParser;
p.addRequired('vineArray',@ismatrix);
p.parse(A);

% some sanity checks for vine array A
if (size(A,1) ~= size(A,2))
    error('transforma:InvalidVineArray','vine array A has to be a quadratic matrix');
end

for jj = 1:1:size(A,1)
    
    if (length(unique(A(1:jj,jj))) ~= jj)
        error('transforma:InvalidVineArray','input A is not a vine array');
    end
    
end % jj

% initialize variables
d = size(A,1);
At = A;

% generate permutation
p = [diag(A) (1:d)'];
pinv = sortrows(p);

% transform A to At
for ii = 1:1:d
    
    At(1:ii,ii) = pinv(A(1:ii,ii),2);
    
end % ii


end