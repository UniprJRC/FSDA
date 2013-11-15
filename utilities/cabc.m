function cabc
%cabc closes all open figures except the one in foreground (the current)
%
%
% Copyright 2008-2011.
% Written by Marco Riani, Domenico Perrotta, Francesca Torti
%            and Vytis Kopustinskas (2009-2010)
%

% Given a vector containing the open figure handles, all_openfigs, the
% figure in foreground is all_openfigs(1)
all_openfigs = findobj(0, 'type', 'figure');
current_fig = all_openfigs(1);

% then, those to delete are:
delete(setdiff(all_openfigs, current_fig));
end