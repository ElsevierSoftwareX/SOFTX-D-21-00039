function bool = cpcheck(family,theta)
% Checks whether the copula parameter theta is valid for a given copula 
% family.
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

switch(lower(family))
    
    case{'gauss'}
        if (theta <= -1 || theta >= 1)
            bool = 0;
        else
            bool = 1;
        end
        
    case{'t'}
        if (theta(1) <= -1 || theta(1) >= 1)
            bool = 0;
        elseif (theta(2) < 1)
            bool = 0;
        else
            bool = 1;
        end
        
    case{'clayton','plackett','surclayton'}
        if (theta < 0)
            bool = 0;
        else
            bool = 1;
        end
        
    case{'gumbel','joe','surgumbel','surjoe'}
        if (theta < 1)
            bool = 0;
        else
            bool = 1;
        end
        
    case{'ind','frank'}
        bool = 1;
        
    case{'amhaq','fgm'}
        if (theta < -1 || theta > 1)
            bool = 0;
        else
            bool = 1;
        end
        
    case{'tawn'}
        if (theta < 0 || theta > 1)
            bool = 0;
        else
            bool = 1;
        end
    
end % switch

end