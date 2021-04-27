function [thetahat, loglik] = copmle(u,family)
% Helper function used in copulaselect for maximum likelihood estimation.
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

tau = corr(u,'type','Kendall');
tau = tau(1,2);

if tau >= 0 % positively dependent data
    
    switch(lower(family))
        
        % copulas implemented in Matlab
        case{'ind'}
            thetahat = 0;
            loglik = log(u(:,1).*u(:,2));
        
        case{'gumbel','clayton','frank'} 
            thetahat = copulafit(family,u);
            loglik = log(copulapdfadv(family,u,thetahat));
            
        case('gauss')
            thetahat = copulafit(family,u);
            thetahat = thetahat(1,2);
            loglik = log(copulapdfadv(family,u,thetahat));
            
        case('t')
            [rho, nu] = copulafit(family,u);
            thetahat = [rho(1,2) nu];
            loglik = log(copulapdfadv(family,u,[rho(1,2) nu]));
        
         % copulas not implemented in Matlab    
        case{'plackett','joe'}
            thetahat = fminsearch(@(theta)mle_copula(u,family,theta),2);
            loglik = log(copulapdfadv(family,u,thetahat));
            
        case{'amhaq','tawn','fgm',}
            theta0 = 0.2;
            thetahat = fminsearch(@(theta)mle_copula(u,family,theta),theta0);
            loglik = log(copulapdfadv(family,u,thetahat));
            
        case{'surclayton','surgumbel'}
            thetahat = copulafit(erase(lower(family),'sur'),1-u);
            loglik = log(copulapdfadv(family,u,thetahat));
            
        case{'surjoe'}
            thetahat = fminsearch(@(theta)mle_copula(1-u,erase(lower(family),'sur'),theta),2);
            loglik = log(copulapdfadv(family,u,thetahat));
            
        % add new copulas for positive dependence here
            %         thetahat_h = fminsearch(@(theta)mle_copula(u,family,theta),XXX);
            %         loglik_h = log(copulapdfadv('XXX',u,thetahat_h));
            %

                  
    end % switch

else % negatively dependent data
    
    switch(lower(family))
        
        case{'ind'}
            thetahat = 0;
            loglik = log(u(:,1).*u(:,2));
        
        case{'gumbel','clayton','tawn','joe','surclayton','surgumbel','surjoe'}
            thetahat = 1;
            loglik = -1000000000;
        
        case{'frank'} 
            thetahat = copulafit(family,u);
            loglik = log(copulapdfadv(family,u,thetahat));
            
        case('gauss')
            thetahat = copulafit(family,u);
            thetahat = thetahat(1,2);
            loglik = log(copulapdfadv(family,u,thetahat));            
            
        case('t')
            [rho, nu] = copulafit(family,u);
            thetahat = [rho(1,2) nu];
            loglik = log(copulapdfadv(family,u,[rho(1,2) nu]));
            
        case{'plackett'} 
            thetahat = fminsearch(@(theta)mle_copula(u,family,theta),0.5);
            loglik = log(copulapdfadv(family,u,thetahat));
            
        case{'amhaq','fgm',}
            theta0 = -0.2;
            thetahat = fminsearch(@(theta)mle_copula(u,family,theta),theta0);
            loglik = log(copulapdfadv(family,u,thetahat));
            
        % add new copulas for negative dependence here
            %         thetahat_h = fminsearch(@(theta)mle_copula(u,family,theta),XXX);
            %         loglik_h = log(copulapdfadv('XXX',u,thetahat_h));
            %
            
    end % switch
    
end


end
