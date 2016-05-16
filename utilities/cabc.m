function cabc()
%cabc closes all open figures except the one in foreground (the current)
%
%
%
%
%<a href="matlab: docsearchFS('cabc')">Link to the help page for this function</a>
%
%
% Required input arguments:
%
%
% Optional input arguments:
%
% Output:
%
% See also:     close
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('cabc')">Link to the help function</a>
% Last modified 06-Feb-2015
%
% Examples:

%{
    %% Plot sin, cos and atan in the interval 0 2pi
    x = 0:pi/1000:2*pi;
    plot(x,sin(x));
    figure
    plot(x,cos(x));
    figure
    plot(x,atan(x))
    cascade
    % Now highlight a figure and then digit
    cabc
    % All the other figures but the one which has been selected will be
    % closed
%}

%{
    % Use of routine cabc inside GUI brushRES
    brushRES
%}

%% Beginning of code
% Given a vector containing the open figure handles, all_openfigs, the
% figure in foreground is all_openfigs(1)
all_openfigs = findobj(0, 'type', 'figure');
current_fig = all_openfigs(1);

% then, those to delete are:
delete(setdiff(all_openfigs, current_fig));
end
%FScategory:UTIGEN
