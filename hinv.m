function hi = hinv(w,u,family,theta)
% Computes the inverse of the h-function of Aas et al (2009), which is the 
% inverse of the conditional distribution C(v|u), i.e. C^-1(v|u). 
%
% call: hi = hinv(w,u,family,theta)
%
% input     w       - conditioned variable 
%           u       - conditioning variable 
%           family  - the copula family: 'gumbel', 'clayton', 'frank', 't',
%                                        'gauss', 'ind', 'amhaq', 'tawn', 
%                                        'fgm', 'plackett', 'joe', 
%                                        'surclayton', 'surgumbel',
%                                        'surjoe'
%           theta   - vector of copula parameters; for t-copula [rho, nu]
%
% output    hi      - the inverse h-function value
% 
%
% References: 
% Aas et al (2009), Pair-copula constructions of multiple dependence, 
% Insurance: Mathematics and Economics, Vol. 44, 182-198.
% Joe (2015), Dependence Modeling with Copulas, CRC Press.
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

% some parsing
p = inputParser;
p.addRequired('w',@isvector);
p.addRequired('u',@isvector);
p.addRequired('family',@isstr);
p.addRequired('theta',@isvector);
p.parse(w,u,family,theta);

% sanity checks
if ~cpcheck(family,theta)
    error('hinv:InvalidParameter',['invalid parameter for ',family,' copula']);
end

% computation 
switch lower(family)
    
    case('gumbel') % Gumbel needs numerical procedure
        if theta > 120 
            theta = 120;
        end
        fun = @(t)(w-copulacdf('Gumbel',[t u],theta)*1/u*(-log(u))^(theta-1)*...
            ((-log(t))^(theta)+(-log(u))^(theta))^(1/theta-1)) ;
        try
            hi = fzero(fun,[0.00001 0.99999]);
        catch
            if abs(fun(0.00001)) < abs(fun(0.99999))
                hi = 0.00001;
            else
                hi = 0.99999;
            end
        end % try
        if hi < 0
            hi = 0;
        elseif hi > 1
            hi = 1;
        end
        
    case('clayton')
        if theta < 0.00001
            theta = 0.00001;
        elseif theta > 150
            theta = 150;
        end
        hi = ((w*u^(theta+1))^(-theta/(1+theta))+1-u^(-theta))^(-1/theta);
        if hi < 0
            hi = 0;
        elseif hi > 1
            hi = 1;
        end
        
    case('frank')
        if theta == 0 
            theta = 0.00001;
        elseif theta < -700
            theta = -700;
        elseif theta > 700
            theta = 700;
        end
        hi = -1/theta*log((-exp(-theta*u)-w*(exp(-theta)-exp(-theta*u)))/...
            (w*exp(-theta*u)-w-exp(-theta*u)));
        if hi < 0
            hi = 0;
        elseif hi > 1
            hi = 1;
        end
        
    case('t')
        if theta(1) == -1
            theta(1) = -0.99999;
        elseif theta(1) == 1
            theta(1) = 0.99999;
        end
        hi = tcdf(tinv(w,theta(2)+1)*...
            sqrt((theta(2)+tinv(u,theta(2))^2)*(1-theta(1)^2)/(theta(2)+1))+...
            theta(1)*tinv(u,theta(2)),theta(2));
        if hi < 0
            hi = 0;
        elseif hi > 1
            hi = 1;
        end
        
    case('gauss')
        if theta == -1
            theta = -0.99999;
        elseif theta == 1
            theta = 0.99999;
        end
        hi = normcdf(norminv(w)*sqrt(1-theta^2)+theta*norminv(u));
        if hi < 0
            hi = 0;
        elseif hi > 1
            hi = 1;
        end
        
    case('ind')
        hi = w;
        if hi < 0
            hi = 0;
        elseif hi > 1
            hi = 1;
        end
        
    case('amhaq') % Ali-Mikhail-Haq needs numerical procedure
        fun = @(t)(w-(t-t*theta*(1-t)*(1-u)-t*u*theta*(1-t))/(1-theta*(1-t)*...
            (1-u))^2);
        try
            hi = fzero(fun,[0.00001 0.99999]);
        catch
            if abs(fun(0.00001)) < abs(fun(0.99999))
                hi = 0.00001;
            else
                hi = 0.99999;
            end
        end
        if hi < 0
            hi = 0;
        elseif hi > 1
            hi = 1;
        end
        
    case('tawn') % Tawn needs numerical procedure
        fun = @(t)(w-t*exp(-theta*log(u)*log(t)/log(u*t))*(1-theta*log(t)*log(t)/(log(u*t)*log(u*t))));
        try
            hi = fzero(fun,[0.00001 0.99999]);
        catch
            if abs(fun(0.00001)) < abs(fun(0.99999))
                hi = 0.00001;
            else
                hi = 0.99999;
            end
        end
        if hi < 0
            hi = 0;
        elseif hi > 1
            hi = 1;
        end
        
    case('joe') % Joe needs numerical procedure
        if theta > 150
            theta = 150;
        end
        fun = @(t)(w-(1-u)^(theta-1)*(1-(1-t)^theta)*...
            ((1-u)^theta+(1-t)^theta-(1-u)^theta*(1-t)^theta)^(1/theta-1));
        try
            hi = fzero(fun,[0.00001 0.99999]);
        catch
            if abs(fun(0.00001)) < abs(fun(0.99999))
                hi = 0.00001;
            else
                hi = 0.99999;
            end
        end
        if hi < 0
            hi = 0;
        elseif hi > 1
            hi = 1;
        end
        
    case('fgm')
        if theta == 0
            theta = 0.00001;
        end
        if u == 0.5
            hi = w;
        else
            hi = (-1-theta*(1-2*u)+sqrt((1+theta*(1-2*u))^2+4*theta*w*(2*u-1)))/(2*theta*(2*u-1));
            if hi < 0
                hi = 0;
            elseif hi > 1
                hi = 1;
            end
        end
        
    case('plackett') % Plackett needs numerical procedure
        if theta > 1000000
            theta = 1000000;
        end
        fun = @(t)(w-0.5+0.5*(1+(-theta-1)*t+(theta-1)*u)/sqrt((1+(theta-1)*(u+t))^2-4*theta*(theta-1)*u*t));
        try
            hi = fzero(fun,[0.00001 0.99999]);
        catch
            if abs(fun(0.00001)) < abs(fun(0.99999))
                hi = 0.00001;
            else
                hi = 0.99999;
            end
        end
        if hi < 0
            hi = 0;
        elseif hi > 1
            hi = 1;
        end
        
    case{'surclayton','surgumbel','surjoe'}
        hi = 1-hinv(1-w,1-u,erase(lower(family),'sur'),theta);
        
 
    % add new copulas here    
%     case('XXX') % numerical procedure
%         fun = @(t)(w-formula);
%         try
%             hi = fzero(fun,[0.0001 0.9999]);
%         catch
%             if abs(fun(0.0001)) < abs(fun(0.9999))
%                 hi = 0.0001;
%             else
%                 hi = 0.9999;
%             end
%         end
%         if hi < 0
%             hi = 0;
%         elseif hi > 1
%             hi = 1;
%         end
%
%     case('XXX') % formula is given
%         hi = formula;
%         if hi < 0
%             hi = 0;
%         elseif hi > 1
%             hi = 1;
%         end        
        
end % switch

% prevent numerical trouble
hi(hi<0.00001) = 0+0.00001;
hi(hi>0.99999) = 1-0.00001;


end