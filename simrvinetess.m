function u = simrvinetess(n,A,family,theta)
% Simulates a sample of size n from a d-dimensional non-simplified r-vine 
% copula with tesselated conditioning spaces.
%
% call: u = simrvinetess(n,A,family,theta)
%
% input     n                   - number of sample points
%           A                   - a vine array; note that a feasible 
%                                 structure has to be used, since the 
%                                 function will not check this
%           family              - a (d-1)x(d-1) cell variable determining 
%                                 the copula families used in the r-vine 
%                                 copula; possible families: 'gumbel', 
%                                 'clayton', 'frank', 't', 'gauss', 'ind', 
%                                 'amhaq', 'tawn', 'fgm', 'plackett', 
%                                 'joe', 'surclayton', 'surgumbel', 
%                                 'surjoe'
%           theta               - a (d-1)x(d-1) cell variable of copula 
%                                 parameters; for t-copula insert [rho nu] 
%                                 in cell element
%
% output    u                   - simulated sample
%
%
% How does it work?
% The function simulates points from a d-dimensional non-simplified r-vine 
% copula of arbitrary structure. The algorithm employed here lends its core 
% from the one given in Joe (2015) for simplified r-vine copulas. Note that 
% for the function to work, the vine array provided by the user has to be a 
% feasible vine array in the first place. The function will not check 
% feasibilty on its own!
%
% Structure of the input is demonstrated for a 5-dimensional non-simplified
% r-vine copula:
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
% , where the numbers correspond to the columns of the output. In this case 
%
%     1 1 2 3 3
%       2 1 2 4
% A =     3 1 2
%           4 1
%             5
%
% is the corresponding vine array.
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
% The inputs for copula parameters theta13|2, theta24|3, theta45|3, 
% theta14|23, theta25|34 and theta15|234 in the theta variable set the 
% break points that define the conditional copula parameters of the 
% corresponding copulas.
% If there is only one conditioning variable (e.g. theta13|2, theta24|3, 
% theta45|3) it has to be a 2xk cell variable, where k >= 2. As the last 
% break point the user always has to include a 1. The first row defines the 
% k break points xi, i=1:1:k, used in ascending order 
% 0 < x1 < x2 <...< xk-1 < xk=1. The second row contains the copula 
% parameters used. For example,
%
%             0.2 0.7 1
% theta13|2 =  2   3  4
%
% Matlab syntax:
%    theta13|2 = {0.2 0.7 1; 2 3 4}
%
% specifies two break points at 0.2 and 0.7. However, in order to
% complete the line segment [0,1], a parameter has to be specified for 
% 0.7 < u <= 1, thus yielding 3 columns in theta13|2.
% If there are two or more conditioning variables (e.g. theta14|23, 
% theta25|34, theta15|234) the break points represent p-dimensional 
% (hyper-)rectangles, where p is the number of conditioning variables. The 
% inputs have to be 2xk cell variables, where k represents the number of 
% (hyper-)rectangles and k >= 2. The first "row" of the cell variable 
% contains the lower left and upper right corners of the (hyper-)rectangles 
% (orientation along the diagonal) as a px2 matrix each. Thus, each corner 
% is represented as a p-dimensional column vector. The second "row" 
% contains the corresponding copula parameter. For example
%
%               [0 0.5   [0.5 1
% theta14|23 =   0  1 ]    0  1]
%                  1         2
%
% Matlab syntax:
%    theta14|23 = {[0 0.5;0 1],[0.5 1; 0 1]; 1, 2}
%
% specifies two "break points". One as the rectangle with lower left corner
% (0,0) and upper right corner (0.5,1) and one as the rectangle with lower
% left corner (0.5,0) and upper right corner (1,1). Thus, the 2-dimensional
% unit cube is partitioned as follows
%
%        (0.5,1)   (1,1)
%            |     |
%      ------+-----*
%      |     |     |
%      |     |     |
%      |     |     |
%      |     |     |
%      +-----*------
%      |     |
%  (0,0) (0.5,0)
%
% ,where + indicates the corners of the first rectangle and * indicates the 
% corners of the second rectangle, respectively.
%
% Note that the user has to input rectangles such that
% 1. the whole p-dimensional unit (hyper-)cube is partitioned by the k
%    rectangles;
% 2. the p-dimensional (hyper-)rectangles do not overlap.
% The function will check these conditions not exhaustively!
%
% Also note that the conditioning variables are assumed to be in ascending 
% order.
%
% In addition to that, the user can specify different copula families at
% the break points. This is done in the family input. For example, 
%            
% family14|23 =  {'gumbel' 'clayton'}                
%
% Matlab syntax:
%    family = {'family12','family23','family34','family35'; 'family13|2','family24|3','family45|3',0; {'gumbel' 'clayton'},'family25|34',0,0;'familiy15|234',0,0,0}
%
% specifies a Gumbel copula at the first break point and a Clayton copula
% at the second break point. If only one family is specified (e.g., 
% family14|23 = 'gumbel'), the procedure will use this family for all break 
% points. Note that in all other cases the number of break points and the 
% number of copulas specified in the family element have to align. 
%
%
% References:
% Joe (2015), Dependence Modeling with Copulas, CRC Press.
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

% some parsing
p = inputParser;
p.addRequired('n',@isscalar);
p.addRequired('A',@ismatrix);
p.addRequired('family',@iscell);
p.addRequired('theta',@iscell);
p.parse(n,A,family,theta);

% some sanity 
if (size(A,1) ~= size(A,2))
    error('simrvinetess:InvalidVineArray','vine array A has to be a quadratic matrix');
end

for jj = 1:1:size(A,1)
    
    if (length(unique(A(1:jj,jj))) ~= jj)
        error('simrvinetess:InvalidVineArray','input A is not a vine array');
    end
    
end % jj

% sanity checks
for ii = 1:1:size(family,1)
    
    for jj = 1:1:size(family,2)-ii+1
        
        if isvector(theta{ii,jj})
            
            if ~cpcheck(family{ii,jj},theta{ii,jj})
                error('simrvinetess:InvalidParameter',['invalid parameter for ',family{ii,jj},' copula at (',num2str(ii),',',num2str(jj),')']);
            end
            
        else % tesselation is given
            
            if ~check_cube(theta{ii,jj})
                
                error('simrvinetess:InvalidTesselation',['invalid tesselation for theta at (',num2str(ii),',',num2str(jj),')']);
                
            end
            
            if ischar(family{ii,jj})
                
                for kk = 1:1:size(theta{ii,jj},2)
                    if ~cpcheck(family{ii,jj},theta{ii,jj}{2,kk})
                        error('simrvinetess:InvalidParameter',['invalid parameter for ',family{ii,jj},' copula at (',num2str(ii),',',num2str(jj),')']);
                    end
                end % kk
                
            else
                
                for kk = 1:1:size(theta{ii,jj},2)
                    if ~cpcheck(family{ii,jj}{kk},theta{ii,jj}{2,kk})
                        error('simrvinetess:InvalidParameter',['invalid parameter for ',family{ii,jj}{kk},' copula at (',num2str(ii),',',num2str(jj),')']);
                    end
                end % kk
                
            end
            
        end
        
    end % jj
    
end % ii

% initialze variables
d = size(A,2);
u = zeros(n,d);
M = zeros(d);

theta_aux = cell(d-1,d-1);
theta_aux(:) = {0};
theta_aux(1,:) = theta(1,:);

fam_aux = cell(d-1,d-1);
fam_aux(:) = {0};
fam_aux(1,:) = family(1,:);

% permute A, such that a_jj = jj
[A,~,pinv] = transforma(A);

% compute matrix M
for jj = 2:1:d
    
    for kk = 1:1:jj-1
        
        M(kk,jj) = max(A(1:kk,jj));
        
    end % kk
    
end % jj

% compute indicator matrix I
I = iarray_rvine(A);

% draw nxd uniform random variables used in the simulation
p = rand(n,d);

for kk = 1:1:n
    
    % initialize arrays q,v,z
    q = zeros(d);
    v = zeros(d);
    z = zeros(d);
    
    % the first sample component
    u(kk,1) = p(kk,1);
    
    % the second sample component
    u(kk,2) = hinv(p(kk,2),p(kk,1),fam_aux{1,1},theta_aux{1,1});
    q(2,2) = p(kk,2);
    if I(1,2) == 1
        v(1,2) = hfunc(u(kk,1),u(kk,2),fam_aux{1,1},theta_aux{1,1});
    end
    
    % now select correct copula parameters
    if isvector(theta{2,1})
        
        theta_aux{2,1} = theta{2,1};
        fam_aux{2,1} = family{2,1};
        
    else
        
        interval = theta{2,1};
        ll = find(u(kk,A(1,3))<=[interval{1,:}],1);
        theta_aux{2,1} = interval{2,ll};
        
        if ischar(family{2,1}) % check, whether multiple families are specified
            fam_aux{2,1} = family{2,1};
        else
            fam_aux{2,1} = family{2,1}{1,ll};
        end
        
    end
    
    % sample components 3 to d
    for jj = 3:1:d
        
        q(jj,jj) = p(kk,jj);
        
        for ll = jj-1:-1:2
            
            if A(ll,jj)==M(ll,jj)
                s = q(ll,A(ll,jj));
            else
                s = v(ll-1,M(ll,jj));
            end
            
            z(ll,jj) = s;
            q(ll,jj) = hinv(q(ll+1,jj),s,fam_aux{ll,jj-ll},theta_aux{ll,jj-ll});
            
        end % ll
        
        q(1,jj) = hinv(q(2,jj),u(kk,A(1,jj)),fam_aux{1,jj-1},theta_aux{1,jj-1});
        u(kk,jj) = q(1,jj);
        v(1,jj) = hfunc(u(kk,A(1,jj)),u(kk,jj),fam_aux{1,jj-1},theta_aux{1,jj-1});
        
        % now select correct copula parameters
        if jj ~= d % stop, if last u component is sampled
            for ii = 1:1:jj-1
                
                if isvector(theta{ii+1,jj-ii})
                    theta_aux(ii+1,jj-ii) = theta(ii+1,jj-ii);
                    fam_aux(ii+1,jj-ii) = family(ii+1,jj-ii);
                else
                    
                    if ii == 1 % only one conditioning variable
                        
                        interval = theta{ii+1,jj-ii};
                        ll = find(u(kk,A(ii,jj))<=[interval{1,:}],1);
                        theta_aux{ii+1,jj-ii} = interval{2,ll};
                        
                        if ischar(family{ii+1,jj-ii}) % check, whether multiple families are specified
                            fam_aux{ii+1,jj-ii} = family{ii+1,jj-ii};
                        else
                            fam_aux{ii+1,jj-ii} = family{ii+1,jj-ii}{1,ll};
                        end
                        
                    else % two or more conditioning variables
                        
                        u_aux2 = u(kk,sort(A(1:ii,jj+1)))';
                        for ll = 1:1:size(theta{ii+1,jj-ii},2) % loop rectangles
                            
                            if all(all([u_aux2 >= theta{ii+1,jj-ii}{1,ll}(:,1) u_aux2 <= theta{ii+1,jj-ii}{1,ll}(:,2)]))
                                
                                theta_aux{ii+1,jj-ii} = theta{ii+1,jj-ii}{2,ll};
                                
                                if ischar(family{ii+1,jj-ii}) % check, whether multiple families are specified
                                    fam_aux{ii+1,jj-ii} = family{ii+1,jj-ii};
                                else
                                    fam_aux{ii+1,jj-ii} = family{ii+1,jj-ii}{1,ll};
                                end
                                
                                break
                                
                            end
                            
                        end % ll
                        
                    end
                    
                end
                
            end % ii
            
        end
        
        for ll = 2:1:jj-1
            
            if I(ll,jj) == 1
                v(ll,jj) = hfunc(z(ll,jj),q(ll,jj),fam_aux{ll,jj-ll},theta_aux{ll,jj-ll});
            end
            
        end % ll
        
    end % jj
    
end % kk

% prevent numerical trouble
if sum(sum(isnan(u))) > 0
    
    u_nan = simrvinebp(sum(sum(isnan(u))),family,theta); 
    u(any(isnan(u),2),:) = u_nan;
    
end 

u(u<=0) = 0+0.00001;
u(u>=1) = 1-0.00001;

% recover original order of variables
u = u(:,pinv(:,2));


end