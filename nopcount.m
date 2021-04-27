function nop = nopcount(family)
% Counts the number of parameters of a given d-dimensional arbitrary vine 
% copula defined in input family.
%
% call: nop = nopcount(family)
%
% input     family  - a (d-1)x(d-1) cell variable determining the copula 
%                     families used in the vine structure
%                     possible families: 'gumbel', 'clayton', 'frank', 't',
%                                        'gauss', 'ind', 'amhaq', 'tawn', 
%                                        'fgm', 'plackett', 'joe', 
%                                        'surclayton', 'surgumbel', 
%                                        'surjoe'
%
% output    nop     - number of parameters in the vine copula
%
%
% Copyright 2020, Maximilian Coblenz
% This code is released under the 3-clause BSD license.
%

% some parsing
p = inputParser;
p.addRequired('family',@iscell);
p.parse(family);

% one parameter for each copula
nop = size(family,1)*(size(family,1)+1)/2;

% now incorporate additional parameters for t-copula
nop = nop + sum(sum(strcmpi(family,'t')));


end