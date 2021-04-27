function h = hfunc(u,v,family,theta)
% Computes the h-function of Aas et al (2009), which is the conditional 
% copula C(u|v). 
%
% call: h = hfunc(u,v,family,theta)
%
% input     u       - conditioned variable (can be nx1 vector)
%           v       - conditioning variable (can be nx1 vector)
%           family  - the copula family: 'gumbel', 'clayton', 'frank', 't',
%                                        'gauss', 'ind', 'amhaq', 'tawn',
%                                        'fgm', 'plackett', 'joe',
%                                        'surclayton', 'surgumbel',
%                                        'surjoe'
%           theta   - vector of copula parameters; for t-copula [rho, nu]
%
% output    h       - the h-function value
%
%
% References: 
% Aas et al (2009), Pair-copula constructions of multiple dependence, 
% Insurance: Mathematics and Economics, Vol. 44, 182-198.
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

% some parsing
p = inputParser;
p.addRequired('u',@isvector);
p.addRequired('v',@isvector);
p.addRequired('family',@isstr);
p.addRequired('theta',@isvector);
p.parse(u,v,family,theta);

% sanity checks
if ~cpcheck(family,theta)
    error('hfunc:InvalidParameter',['invalid parameter for ',family,' copula']);
end

% computation 
switch lower(family)
    
    case('gumbel')
        h = copulacdf('Gumbel',[u v],theta).*1./v.*(-log(v)).^(theta-1).*...
            ((-log(u)).^(theta)+(-log(v)).^(theta)).^(1/theta-1);
        
    case('clayton')
        h = v.^(-theta-1) .*(u.^(-theta)+v.^(-theta)-1).^(-1-1/theta);
        
    case('frank')
        h = exp(-theta*v).*(exp(-theta*u)-1)./...
            (exp(-theta)-1+(exp(-theta*u)-1).*(exp(-theta*v)-1));
        
    case('t')
        h = tcdf((tinv(u,theta(2))-theta(1).*tinv(v,theta(2)))./...
            sqrt((theta(2)+tinv(v,theta(2)).^2).*(1-theta(1)^2)./(theta(2)+1)),theta(2)+1);
        
    case('gauss')
        h = normcdf((norminv(u)-theta*norminv(v))./sqrt(1-theta^2));
        
    case('ind')
        h = u;
        
    case('amhaq')
        h = (u-u.*theta.*(1-u).*(1-v)-u.*v.*theta.*(1-u))./(1-theta.*(1-u).*(1-v)).^2;
        
    case('tawn')
        h = u.*exp(-theta*log(u).*log(v)./log(u.*v)).*(1-theta*log(u).*log(u)./(log(u.*v).*log(u.*v)));
        
    case('joe')
        h = (1-v).^(theta-1).*(1-(1-u).^theta).*...
            ((1-u).^theta+(1-v).^theta-(1-u).^theta.*(1-v).^theta).^(1/theta-1);
        
    case('fgm')
        h = u.*(1+theta*(1-u).*(1-2*v));
        
    case('plackett')
        h = 0.5 - 0.5*(1+(-theta-1)*u+(theta-1)*v)./...
            sqrt((1+(theta-1)*(u+v)).^2-4*theta*(theta-1)*u.*v);
        
    case{'surclayton','surgumbel','surjoe'}
        h = 1-hfunc(1-u,1-v,erase(lower(family),'sur'),theta);
        
        
    % add new copulas here
%     case('XXX')
%         h = formula;
        
end % switch

% prevent numerical trouble
h(h<0.00001) = 0+0.00001;
h(h>0.99999) = 1-0.00001;


end