function tc = thetalink(beta,family)
% Helper function that calculates the corresponding copula parameter value 
% in simrvinens for restricting parameter values when using function 
% handles.
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

switch lower(family)
    
    case{'gauss','amhaq','fgm'}
        tc = (exp(beta)-1)/(exp(beta)+1);
        
    case('t')
        tc = [(exp(beta)-1)/(exp(beta)+1) exp(beta)];
        
    case{'clayton','plackett','surclayton'}
        tc = exp(beta);
        
    case{'gumbel','joe','surgumbel','surjoe'}
        tc = 1 + exp(beta);
        
    case('frank')
        tc = beta;
        
    case('tawn')
        tc = abs((exp(beta)-1)/(exp(beta)+1));
        
end % switch


end