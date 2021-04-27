function bool = check_cube(rectangles)
% Checks whether the hyperrectangles provided in cell array rectangles
% 1. cover the unit hypercube and
% 2. do not overlap.
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

bool = 1;

if isscalar(rectangles{1,1}) % edge case line
    
    % zero at the beginning not allowed, line has to end with a 1
    if rectangles{1,1} == 0 || rectangles{1,size(rectangles,2)} ~= 1
        
        bool = 0;
        return
        
    end
    
    for ii = 1:1:size(rectangles,2)-1
        
        % check that line is ordered
        if rectangles{1,ii} > rectangles{1,ii+1}
            
            bool = 0;
            return
            
        end
        
    end % ii
    
else % all other cases (square, cube, hypercube)
    
    % check volume is exactly 1 and that lower left corner and upper right
    % corner are present
    check_0 = 0;
    check_1 = 0;
    vol = 0;
    for ii = 1:1:size(rectangles,2)
        
        cube_part = rectangles{1,ii};
        
        % cube_part is not ordered (lower left corner has to be in the left
        % column, upper right corner in the right column)
        if any(cube_part(:,1) >= cube_part(:,2))
            
            bool = 0;
            return
            
        end
        
        if isequal(cube_part(:,1),zeros(size(cube_part(:,1),1),1))
            
            check_0 = 1;
            
        end
        
        if isequal(cube_part(:,2),ones(size(cube_part(:,2),1),1))
            
            check_1 = 1;
            
        end
        
        vol = vol + prod(cube_part(:,2) - cube_part(:,1));
        
    end % ii
    
    % lowermost corner or uppermost corner is missing
    if check_0 == 0 || check_1 == 0
        
        bool = 0;
        return
        
    end
    
    % volume does not add up to 1
    if vol > 1 + eps || vol < 1 - eps
        
        bool = 0;
        return
        
    end
    
    % check no cube parts are equal
    for ii = 1:1:size(rectangles,2)-1

            for jj = ii+1:1:size(rectangles,2)
                
                if isequal(rectangles{1,ii}, rectangles{1,jj})
                    
                    bool = 0;
                    return
                    
                end
                
            end % jj
            
    end % ii
     
end


end