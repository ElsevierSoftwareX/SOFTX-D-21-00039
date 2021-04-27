function negll = mle_copula(u,family,theta)
% Helper function for fminsearch in order to perform Maximum Likelihood 
% Estimation for the parameter theta of the copula given in family.
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

switch(lower(family))
    
    case('amhaq')
        ll = (1-theta+2*theta*u(:,1).*u(:,2)-theta*(1-theta)*(1-u(:,1)).*(1-u(:,2)))./...
            (1-theta*(1-u(:,1)).*(1-u(:,2))).^3;
        
        if (theta > 1 || theta < -1)
            negll = -1*sum(log(ll)) + 1000000000;
        else
            negll = -1*sum(log(ll));
        end
        
    case('fgm')
        ll = 1+theta*(1-2*u(:,1)).*(1-2*u(:,2));
        
        if (theta > 1 || theta < -1)
            negll = -1*sum(log(ll)) + 1000000000;
        else
            negll = -1*sum(log(ll));
        end
        
    case('joe')
        ll = ((1-u(:,1)).^theta+(1-u(:,2)).^theta-(1-u(:,1)).^theta.*(1-u(:,2)).^theta).^(1/theta-2).*...
            (1-u(:,1)).^(theta-1).*(1-u(:,2)).^(theta-1).*(theta-1+(1-u(:,1)).^theta+(1-u(:,2)).^theta-(1-u(:,1)).^theta.*(1-u(:,2)).^theta);
        
        if (theta < 1)
            negll = -1*sum(log(ll)) + 1000000000;
        else
            negll = -1*sum(log(ll));
        end
        
    case('plackett')
        ll = theta*(1+(theta-1)*(u(:,1)+u(:,2)-2*u(:,1).*u(:,2)))./...
            ((1+(theta-1)*(u(:,1)+u(:,2))).^2-4*theta*(theta-1)*u(:,1).*u(:,2)).^(3/2);
        
        if (theta < 0)
            negll = -1*sum(log(ll)) + 1000000000;
        else
            negll = -1*sum(log(ll));
        end
        
    case('tawn')
        ll = exp(-theta*log(u(:,1)).*log(u(:,2))./log(u(:,1).*u(:,2))).*...
            ((1-theta*(log(u(:,2))).^2./(log(u(:,1).*u(:,2))).^2).*(1-theta*(log(u(:,1))).^2./(log(u(:,1).*u(:,2))).^2)-...
            theta*(2*log(u(:,1)).*log(u(:,2)))./(log(u(:,1).*u(:,2))).^3);
        
        if (theta < 0 || theta > 1)
            negll = -1*sum(log(ll)) + 1000000000;
        else
            negll = -1*sum(log(ll));
        end
        
        % add new copula here
        % case('XXX')
        %         ll = XXX;
        %
        %         if (theta < XXX || theta > XXX)
        %             negll = -1*sum(log(ll)) + 1000000000;
        %         else
        %             negll = -1*sum(log(ll));
        %         end
         
end


end