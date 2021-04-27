function list = primmaxst(A)
% Implements Prim's algorithm (1957) for adjacency matrix A to determine a 
% maximum(!) spanning tree.
%
% call: list = primmaxst(A)
%
% input     A       - nxn adjacency matrix of graph
%
% output    list    - a list of the n-1 edges in the maximum spanning tree
%
%
% References:
% Prim, R. C. (1957), Shortest Connection Networks And Some 
% Generalizations, Bell System Technical Journal, Vol. 36(6), 1389–1401.
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

% some parsing
p = inputParser;
p.addRequired('A',@ismatrix);
p.parse(A);

% sanity checks
if norm(A-A','fro') ~= 0
    error('primmaxst:MatrixAsymmetric','adjacency matrix has to be symmetric');
end
if sum(A(:)<0) ~= 0
    error('primmaxst:InvalidMatrix','adjacency matrix may contain only positive weights');
end

% initialize variables
n = size(A,1);
list = zeros(n-1,2);

n_intree = 2;
n_notintree = n-n_intree;
n_edges = 1;

% start with a pair of highest weight
B = triu(A,1);
[~, max_idx] = max(B(:));
[list(1,1),list(1,2)]=ind2sub(size(B),max_idx);

% initialize lists to keep track which vertices are already in the tree
% and which vertices are not 
intree = [list(1,1);list(1,2);zeros(1,n-2)'];
notintree = [1:list(1,1)-1 list(1,1)+1:list(1,2)-1 list(1,2)+1:n]';

% now add all other nodes to tree, consecutively
while n_intree < n
    
    % determine largest weight from a vertex in tree to the other vertices
    % not yet in the tree
    maxcost = -1;
    for kk = 1:n_intree
        
        for mm = 1:n_notintree
            
            i = intree(kk);  
            j = notintree(mm);
            
            if A(i,j) > maxcost 
                
                maxcost = A(i,j); 
                iaux = i; 
                jaux = j;
                
            end 
            
        end % mm
        
    end % kk
    
    % add new edge to the tree
    n_edges = n_edges + 1;
    list(n_edges,1) = iaux;    
    list(n_edges,2) = jaux;
    
    % update intree and notintree variables
    intree(n_edges+1) = jaux;
    notintree = notintree(notintree~=jaux);
    
    % update loop variables
    n_notintree = n_notintree - 1; 
    n_intree = n_intree + 1;
    
end % while


end