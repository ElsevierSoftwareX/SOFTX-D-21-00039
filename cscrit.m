function val = cscrit(crit,ll,nop)
% Helper function that computes the selection criterion for function 
% copulaselect.
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

switch lower(crit)
    
    case('aic')
        val = -2*sum(ll) + 2*nop;
        
    case('bic')
        val = -2*sum(ll) + log(size(ll,1))*nop;
        
    case('sll')
        val = -1*sum(ll);
        
end % switch


end