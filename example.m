%% settings
clear; clc;
rng(123456789);


%% basic functionalities

% build a 5-dimensional R-vine
A = [1 1 2 3 3;
     0 2 1 2 2;
     0 0 3 1 4;
     0 0 0 4 1;
     0 0 0 0 5
    ];
family1 = {'gauss', 'surjoe', 'clayton', 'gumbel';'clayton', 'clayton', 'frank',0;'amhaq','clayton',0,0;'t',0,0,0};
theta1 = {-0.7,2,3,4;2,3,5,0;0.8,3,0,0;[-0.4 5],0,0,0};

% calculate pdf at points (0.1,0.2,0.3,0.4,0.5) and (0.3,0.4,0.5,0.6,0.7)
points = [0.1,0.2,0.3,0.4,0.5;0.3,0.4,0.5,0.6,0.7];

valuepdf = vinepdf(points,A,family1,theta1);

% calculate cdf at these points
valuecdf = vinecdf(points,A,family1,theta1);


%% simulation

% sample from the simplified R-vine above and plot
u = simrvine(1000,A,family1,theta1,0);

figure
plotmatrix(u);

% sample from a non-simplified vine  
family2 = {'gauss','clayton','clayton','gumbel'; 'gauss','gumbel','gauss',0;'gauss','clayton',0,0;'gauss',0,0,0};
theta2 = {-0.6,3,2,3;0.6,'@(u)(thetalink((4*u(3) - 2).^2,''gumbel''))','@(u)(thetalink((4*u(3) - 2).^2,''gauss''))',0; '@(u)(thetalink((4*u(2)^2+u(3) - 2).^2,''gauss''))',4,0,0;'@(u)(thetalink((4*u(2)+u(3)+u(4) - 2).^2,''gauss''))',0,0,0};

u2 = simrvinens(1000,A,family2,theta2,0);

figure
plotmatrix(u2);

% sample from a non-simplified vine with tessellated conditioning spaces
family3 = {'gauss','clayton','clayton','gumbel'; 'gauss',{'gumbel' 'clayton' 'gumbel'},'gauss',0;{'gauss' 'frank'},'clayton',0,0;'gauss',0,0,0};
theta3 = {-0.6,3,2,3;0.6,{0.2 0.7 1; 2 3 4},{0.2 0.5 1; 0.5 -0.7 -0.1},0; {[0 0.5;0 1],[0.5 1; 0 1]; -0.5, 0.5},4,0,0;{[0 0.5;0 1;0 1],[0.5 1;0 1;0 1]; -0.5, 0.5},0,0,0};

u3 = simrvinetess(1000,A,family3,theta3);

figure
plotmatrix(u3);

figure
subplot(3,1,1)
plotmatrix(u);
title('simplified vine')
subplot(3,1,2)
plotmatrix(u2);
title('non-simplified vine')
subplot(3,1,3)
plotmatrix(u3);
title('non-simplified vine with tessellated conditioning spaces')


%% estimation

% estimate with vine structure given
familyssp = {'AIC'};
[thetahat,loglik,~,famhat] = ssp(u,A,familyssp);
uhat = simrvine(1000,A,famhat,thetahat);

figure
subplot(1,2,1)
plotmatrix(u)
title('original data')
subplot(1,2,2)
plotmatrix(uhat)
title('estimated data')

% also estimate the vine structure
[Ahat2, famhat2, thetahat2] = rvineselect(u);
uhat2 = simrvine(1000,Ahat2,famhat2,thetahat2);

figure
subplot(1,2,1)
plotmatrix(u)
title('original data')
subplot(1,2,2)
plotmatrix(uhat2)
title('estimated data')


%% testing

% conduct a goodness-of-fit test of type A2 from Berg (2009) for the
% estimated model on the original data
[tstat, pval] = goftest_a2(u,A,famhat,thetahat,1);

% conduct a goodness-of-fit test of type A4 from Berg (2009) for the
% estimated model on the original data
[tstat2, pval2] = goftest_a4(u,A,famhat,thetahat,1);

% compare the two estimated models with the test by Nikoloulopoulos & Karlis (2008)
[tstat3, pval3] = nktest(u,1,famhat,thetahat,A,famhat2,thetahat2,Ahat2);

% now use a Vuong model comparison test for non-nested models 
ll1 = llrvine(u,A,famhat,thetahat); % loglikelihoods of first model
ll2 = llrvine(u,Ahat2,famhat2,thetahat2); % loglikelihoods of second model
[pval4, tstat4] = vuongtest(ll1,ll2,nopcount(famhat),nopcount(famhat2));


%% further functionalities

% simulate any bivariate copula from the toolbox
u = copulasim('surgumbel',4,1000);

figure
plot(u(:,1),u(:,2),'k.')

% rank-transform the data
u2 = pobs(u);

figure
plot(u2(:,1),u2(:,2),'k.')

figure
subplot(1,2,1)
plotmatrix(u)
title('original sample')
subplot(1,2,2)
plotmatrix(u2)
title('rank-transformed sample')

% and calculate the empirical copula
cophat = ecopula(u2);

% create a generic (ordered) D- or C-vine array
AC = cdvinearray('c',10);
AD = cdvinearray('d',15);

