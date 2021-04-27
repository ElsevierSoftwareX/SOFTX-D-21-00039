function [A,fam,theta] = rvineselect(u,varargin)
% Estimates an arbitrary simplified r-vine copula from data u according to 
% the algorithm in Diﬂmann et al (2013).
%
% call: [A,fam,thetahat] = rvineselect(u[,crit])
%
% input     u                   - nxp data matrix of pseudo-observations
%           crit (optional)     - the selection criterion: 'aic' (Akaike's 
%                                 information criterion), 'bic' (Bayesian 
%                                 information criterion), 'sll' (sum of 
%                                 loglikelihoods); default is 'aic'
%
% output    A               - the estimated vine array
%           fam             - a (d-1)x(d-1) cell variable of copula 
%                             families of the r-vine structure
%           thetahat        - a (d-1)x(d-1) cell variable of estimated 
%                             copula parameters corresponding to the 
%                             copulas in fam
%
%
% How does it work?
% This function estimates an arbitrary simplified r-vine. The following 
% gives an example of output interpretation:
% 
% Let's consider a 5-dimensional sample u and let the output A be 
% 
%     1 1 2 3 3
%     0 2 1 2 4
% A = 0 0 3 1 2
%     0 0 0 4 1
%     0 0 0 0 5
%
% This represents the following simplified r-vine copula:
%
%           4
%          /
% 1 - 2 - 3 
%          \
%           5
%
% 12 - 23 - 34 - 35
%
% 13|2 - 24|3 - 45|3
%
% 14|23 - 25|34
%
% 15|234
%
% , where the numbers correspond to the columns of input u. 
%
% The output fam stores the copula families like this:
%
%          family12     family23    family34   family35
% family = family13|2   family24|3  family45|3    0
%          family14|23  family25|34    0          0
%          family15|234    0           0          0
%
% The output thetahat is structured in the same fashion.
%
%
% References:
% Diﬂmann et al (2013), Selecting and Estimating Regular Vine Copulae and
% Application to Financial Returns, Computational Statistics and Data
% Analysis, Vol. 59, 52-69.
% Czado (2019), Analyzing Dependent Data with Vine Copulas, Springer.
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

% some parsing
p = inputParser;
p.addRequired('u',@ismatrix);
p.addOptional('crit',0,@isstr);
p.parse(u,varargin{:});

% sanity checks
if size(u,2) < 3
    error('rvineselect:InvalidNumberOfDimensions','r-vine selection is only possible for 3 or more dimensions');
end

% initialize variables
if nargin > 2
    crit = varargin{1};
else
    crit = 'aic';
end

d = size(u,2);

% the vine array
A = zeros(d);

% the estimated copula families in the vine
fam = cell(d-1,d-1);
fam(:) = {0};

% the estimated copula parameters in the vine
theta = cell(d-1,d-1);
theta(:) = {0};

% used to save all the edges in the vine
edge_array = cell(d-1,1);
for ii = 1:1:d-1
    
    edge_array{ii,1} = cell(d-ii,1);
    
end % ii


% 1. Estimate first tree T_1

% 1.1 estimate empirical Kendall's tau and prepare it as the adjacency  
%     matrix for Prim's algorithm
taumat = abs(corr(u,'type','Kendall'));
tauadj = taumat-eye(d);

% 1.2 now use Prim's algorithm to determine the maximal spanning tree
mst = primmaxst(tauadj);

% 1.3 estimate copula families and parameters in first tree 
%     and transform pseudo-observations for next r-vine level
for ii = 1:1:d-1

    mst_sort = sort(mst(ii,:));

    % estimate copula family and parameter(s)
    [fam_aux,theta_aux,~,~] = copulaselect(u(:,mst_sort),crit);

    % transform pseudo-observations
    v1 = hfunc(u(:,mst_sort(1)),u(:,mst_sort(2)),fam_aux,theta_aux);
    v2 = hfunc(u(:,mst_sort(2)),u(:,mst_sort(1)),fam_aux,theta_aux);

    % save everything for this edge in edge_array 
    % each edge has the following structure:
    % {conditioned set, 
    %  conditioning set, 
    %  family, 
    %  theta, 
    %  [pseudo-obs for conditioned_set(1)|(conditioned_set(2) & conditioning_set),
    %   pseudo-obs for conditioned_set(2)|(conditioned_set(1) & conditioning_set)]}
    edge_array{1,1}{ii,1} = {mst_sort,[],fam_aux,theta_aux,[v1 v2]};

end % ii

% 2. Estimate trees 2 to d-1
for ii = 2:1:d-1
    
    % 2.1 prepare adjacency matrix for Prim's algorithm
    tauadj = zeros(d+1-ii);
    for jj = 1:1:d-ii+1
        
        edge_jj = edge_array{ii-1,1}{jj,1};
        for kk = jj+1:1:d-ii+1
            
            edge_kk = edge_array{ii-1,1}{kk,1};
            conditioning_vars = intersect(union(edge_jj{1},edge_jj{2}), union(edge_kk{1},edge_kk{2}));
            
            if length(conditioning_vars) == ii-1 % this is a reformulated proximity condition
                
                idx_jj = edge_jj{1}==setdiff(conditioning_vars,edge_jj{2});
                idx_kk = edge_kk{1}==setdiff(conditioning_vars,edge_kk{2});
                tobs = [edge_jj{5}(:,idx_jj) edge_kk{5}(:,idx_kk)];
                ktau = abs(corr(tobs,'type','Kendall'));
                tauadj(jj,kk) = ktau(1,2);
                tauadj(kk,jj) = ktau(1,2);
                
            end 
            
        end % kk
        
    end % jj
    
    % 2.2 determine the maximal spanning tree
    mst = primmaxst(tauadj);
    
    % 2.3 estimate copula families and parameters
    %     and transform pseudo-observations for next r-vine level
    for jj = 1:1:d-ii
        
        edge1 = edge_array{ii-1,1}{mst(jj,1),1};
        edge2 = edge_array{ii-1,1}{mst(jj,2),1};
        
        % compute conditioning set and conditioned set
        conditioning_vars = intersect(union(edge1{1},edge1{2}), union(edge2{1},edge2{2}));
        conditioned_vars = setxor(union(edge1{1},edge1{2}), union(edge2{1},edge2{2}));
        
        % check which data has to be used for the further estimation
        if ismember(conditioned_vars(1),edge1{1})
            idx1 = edge1{1}==conditioned_vars(1);
            idx2 = edge2{1}==conditioned_vars(2);
            tobs = [edge1{5}(:,idx1) edge2{5}(:,idx2)];
        else
            idx1 = edge1{1}==conditioned_vars(2);
            idx2 = edge2{1}==conditioned_vars(1);
            tobs = [edge2{5}(:,idx2) edge1{5}(:,idx1)];
        end
    
        % estimate copula family and parameter(s)
        [fam_aux,theta_aux,~,~] = copulaselect(tobs,crit);

        % transform pseudo-observations
        v1 = hfunc(tobs(:,1),tobs(:,2),fam_aux,theta_aux);
        v2 = hfunc(tobs(:,2),tobs(:,1),fam_aux,theta_aux);

        % save everything for this edge in edge_array 
        edge_array{ii,1}{jj,1} = {conditioned_vars,conditioning_vars,fam_aux,theta_aux,[v1 v2]};
    
    end % jj
    
end % ii

% 3. Retrieve vine array A and the estimated copulas from edge_array
for ii = d-1:-1:1 % vine trees from bottom to top

    edge = edge_array{ii,1}{1};
    lead = edge{1}(1);
    
    % update A, fam, and theta
    A(ii+1,ii+1) = lead;
    A(ii,ii+1) = edge{1}(2);
    fam{ii,1} = edge{3};
    theta{ii,1} = edge{4};
    
    for jj = ii-1:-1:1 % vine trees above current vine tree
        
        for kk = 1:1:length(edge_array{jj,1}) % edges in current tree
            
            % check whether lead element is in the conditioned set of the
            % current edge
            if ismember(lead,edge_array{jj,1}{kk,1}{1}) 
                
                % update A, fam, and theta
                A(jj,ii+1) = setdiff(edge_array{jj,1}{kk,1}{1},lead);
                fam{jj,ii-jj+1} = edge_array{jj,1}{kk,1}{3};
                theta{jj,ii-jj+1} = edge_array{jj,1}{kk,1}{4};
                
                % delete edge from vine
                edge_array{jj,1}(kk,:) = [];
                
                break
                
            end
            
        end % kk
        
    end % jj
    
    if ii == 1
        
        A(ii,ii) = setdiff(1:d, diag(A));
        
    end
    
end % ii


end