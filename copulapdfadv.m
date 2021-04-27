function pdf = copulapdfadv(family,u,theta)
% Computes the pdf of a copula at u. 
%
% call: pdf = copulapdfadv(family,u,theta)
%
% input     family  - the copula family: 'gumbel', 'clayton', 'frank', 't',
%                                        'gauss', 'ind', 'amhaq', 'tawn', 
%                                        'fgm', 'plackett', 'joe', 
%                                        'surclayton', 'surgumbel',
%                                        'surjoe'
%           u       - nx2 matrix of points to be evaluated
%           theta   - copula parameter: for t-copula [rho, nu]
%
% output    pdf     - the pdf at u
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

% some parsing
p = inputParser;
p.addRequired('family',@isstr);
p.addRequired('u',@ismatrix);
p.addRequired('theta',@isvector);
p.parse(family,u,theta);

% sanity checks
if ~cpcheck(family,theta)
    error('copulapdfadv:InvalidParameter',['invalid parameter for ',family,' copula']);
end

% compute pdf
switch lower(family)
    
    case{'gauss','gumbel','clayton','frank'} 
        pdf = copulapdf(family,u,theta);
        
    case('t')
        pdf = copulapdf(family,u,theta(1),theta(2));
    
    case('ind')
        pdf = ones(size(u,1),1);
        
    case('amhaq')
        pdf = (1-theta+2*theta*u(:,1).*u(:,2)-theta*(1-theta)*(1-u(:,1)).*(1-u(:,2)))./...
              (1-theta*(1-u(:,1)).*(1-u(:,2))).^3;
          
    case('tawn')
        pdf = exp(-theta*log(u(:,1)).*log(u(:,2))./log(u(:,1).*u(:,2))).*...
              ((1-theta*(log(u(:,2))).^2./(log(u(:,1).*u(:,2))).^2).*(1-theta*(log(u(:,1))).^2./(log(u(:,1).*u(:,2))).^2)-...
              theta*(2*log(u(:,1)).*log(u(:,2)))./(log(u(:,1).*u(:,2))).^3);
          
    case('joe')
        pdf = ((1-u(:,1)).^theta+(1-u(:,2)).^theta-(1-u(:,1)).^theta.*(1-u(:,2)).^theta).^(1/theta-2).*...
              (1-u(:,1)).^(theta-1).*(1-u(:,2)).^(theta-1).*(theta-1+(1-u(:,1)).^theta+(1-u(:,2)).^theta-...
              (1-u(:,1)).^theta.*(1-u(:,2)).^theta);
          
    case('fgm')
        pdf = 1+theta*(1-2*u(:,1)).*(1-2*u(:,2));
        
    case('plackett')
        pdf = theta*(1+(theta-1)*(u(:,1)+u(:,2)-2*u(:,1).*u(:,2)))./...
              ((1+(theta-1)*(u(:,1)+u(:,2))).^2-4*theta*(theta-1)*u(:,1).*u(:,2)).^(3/2);
          
    case{'surclayton','surgumbel','surjoe'}
        pdf = copulapdfadv(erase(lower(family),'sur'),1-u,theta);
      
          
   % add new copulas here      
%     case('XXX')
%         pdf = formula;
    
end % switch


end