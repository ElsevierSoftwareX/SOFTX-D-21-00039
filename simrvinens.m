function u = simrvinens(n,A,family,theta,varargin)
% Simulates a sample of size n from a d-dimensional non-simplified r-vine 
% copula using function handles for conditional copula parameters.
%
% call: u = simrvinens(n,A,family,theta[,parallel])
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
%           parallel (optional) - switch parallelization on (=1) 
%                                 or off (=0); default: 0
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
% theta14|23, theta25|34 and theta15|234 in the theta variable define the 
% functions used to compute the conditional copula parameters of the 
% corresponding copulas. They have to be either a String that defines a 
% function handle - e.g., for a quadratic function '@(x)3*x(1)^2+3' - or a 
% scalar value, in which case it is simplified. In case a function handle 
% is used, always specify the variables being used - even when it is only 
% the first one.
%
% Note that 
% 1. the function handle always gets handed the full vector of variables. 
% So the variables which appear in the condition function always have to  
% be specified - even in the case when only the first one is used! Make 
% sure you choose the correct conditioning variables, i.e. indices, in the 
% function handle corresponding to your vine structure, as the function 
% will not check correctness on its own!
% 2. the possible values returned by the function handle have to be
% in line with the boundaries of the copula parameter of the employed  
% copula. The function will not check this! If you are not sure, whether
% the parameter boundaries are fulfilled, wrap your function by thetalink,
% which will transform function values into the correct parameter space,  
% e.g. use '@(x)(thetalink(3*x(1)^2+3,family))' instead of 
% '@(x)(3*x(1)^2+3,family)'. 
% 3. if the condition includes more than one variable this has to be 
% reflected in the function handle. For example use '@(x)3*x(1)^2+x(2)+3'
% for two conditioning variables.
% 4. even if the condition includes more than one variable, you can choose
% to make it depend on less than the full amount. The variables, which
% are not used can be left out in the function itself. For example use 
% '@(x)3*x(1)^2+x(3)+3', if you want the condition depend on variable one
% and three only.
% 
% Also note that the conditioning variables are assumed to be in ascending 
% order.
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
p.addOptional('para',0,@isscalar);
p.parse(n,A,family,theta,varargin{:});

% some sanity checks for vine array A
if (size(A,1) ~= size(A,2))
    error('simrvinens:InvalidVineArray','vine array A has to be a quadratic matrix');
end

for jj = 1:1:size(A,1)
    
    if (length(unique(A(1:jj,jj))) ~= jj)
        error('simrvinens:InvalidVineArray','input A is not a vine array');
    end
    
end % jj

% sanity checks
for ii = 1:1:size(family,1)
    
    for jj = 1:1:size(family,2)-ii+1
        
        if isscalar(theta{ii,jj})
            if ~cpcheck(family{ii,jj},theta{ii,jj})
                error('simrvinens:InvalidParameter',['invalid parameter for ',family{ii,jj},' copula at (',num2str(ii),',',num2str(jj),')']);
            end
        end
        
    end % jj
    
end % ii

% determine parallelization mode
para = 0;
if nargin > 4
    if varargin{1} ~= 1 && varargin{1} ~= 0
        error('simrvinens:InvalidParallelizationMode','optional input argument parallel has to be 0 or 1')
    end
    para = varargin{1};
end

% initialze variables
d = size(A,2);
u = zeros(n,d);
M = zeros(d);

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

if para
    
    ppool = gcp; 
    poolsize = ppool.NumWorkers;
    
    parfor(kk = 1:n,poolsize)
        
        u_aux = zeros(1,d);
        
        % initialize arrays a,b,c anew
        q = zeros(d);
        v = zeros(d);
        z = zeros(d);

        % the first sample component
        u_aux(1) = p(kk,1);

        % the second sample component
        u_aux(2) = hinv(p(kk,2),p(kk,1),family{1,1},theta{1,1});
        q(2,2) = p(kk,2);
        if I(1,2) == 1
            v(1,2) = hfunc(u_aux(1),u_aux(2),family{1,1},theta{1,1});
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

                if ischar(theta{ll,jj-ll}) % function for theta of conditional copula is given

                    func = str2func(theta{ll,jj-ll});
                    q(ll,jj) = hinv(q(ll+1,jj),s,family{ll,jj-ll},func(u_aux(:)));

                else % theta is a scalar

                    q(ll,jj) = hinv(q(ll+1,jj),s,family{ll,jj-ll},theta{ll,jj-ll});

                end

            end % ll

            q(1,jj) = hinv(q(2,jj),u_aux(A(1,jj)),family{1,jj-1},theta{1,jj-1});
            u_aux(jj) = q(1,jj);
            v(1,jj) = hfunc(u_aux(A(1,jj)),u_aux(jj),family{1,jj-1},theta{1,jj-1});

            for ll = 2:1:jj-1

                if I(ll,jj) == 1

                    if ischar(theta{ll,jj-ll}) % function for theta of conditional copula is given

                        func = str2func(theta{ll,jj-ll});
                        v(ll,jj) = hfunc(z(ll,jj),q(ll,jj),family{ll,jj-ll},func(u_aux(:)));

                    else % theta is a scalar

                        v(ll,jj) = hfunc(z(ll,jj),q(ll,jj),family{ll,jj-ll},theta{ll,jj-ll});

                    end

                end

            end % ll

        end % jj
        
        u(kk,:) = u_aux;
        
    end
    
else

    for kk = 1:1:n

        % initialize arrays a,b,c anew
        q = zeros(d);
        v = zeros(d);
        z = zeros(d);

        % the first sample component
        u(kk,1) = p(kk,1);

        % the second sample component
        u(kk,2) = hinv(p(kk,2),p(kk,1),family{1,1},theta{1,1});
        q(2,2) = p(kk,2);
        if I(1,2) == 1
            v(1,2) = hfunc(u(kk,1),u(kk,2),family{1,1},theta{1,1});
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

                if ischar(theta{ll,jj-ll}) % function for theta of conditional copula is given

                    func = str2func(theta{ll,jj-ll});
                    q(ll,jj) = hinv(q(ll+1,jj),s,family{ll,jj-ll},func(u(kk,:)));

                else % theta is a scalar

                    q(ll,jj) = hinv(q(ll+1,jj),s,family{ll,jj-ll},theta{ll,jj-ll});

                end

            end % ll

            q(1,jj) = hinv(q(2,jj),u(kk,A(1,jj)),family{1,jj-1},theta{1,jj-1});
            u(kk,jj) = q(1,jj);
            v(1,jj) = hfunc(u(kk,A(1,jj)),u(kk,jj),family{1,jj-1},theta{1,jj-1});

            for ll = 2:1:jj-1

                if I(ll,jj) == 1

                    if ischar(theta{ll,jj-ll}) % function for theta of conditional copula is given

                        func = str2func(theta{ll,jj-ll});
                        v(ll,jj) = hfunc(z(ll,jj),q(ll,jj),family{ll,jj-ll},func(u(kk,:)));

                    else % theta is a scalar

                        v(ll,jj) = hfunc(z(ll,jj),q(ll,jj),family{ll,jj-ll},theta{ll,jj-ll});

                    end

                end

            end % ll

        end % jj

    end % kk

end

% prevent numerical trouble
if sum(sum(isnan(u))) > 0
    
    u_nan = simrvinens(sum(sum(isnan(u))),family,theta); 
    u(any(isnan(u),2),:) = u_nan;
    
end 

u(u<=0) = 0+0.00001;
u(u>=1) = 1-0.00001;

% recover original order of variables
u = u(:,pinv(:,2));

% shut down parallel pool
if para
    delete(ppool);
end


end