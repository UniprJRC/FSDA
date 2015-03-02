function [out1,out2,out3] = ginputFS(arg1,varargin)
%ginputFS extends the MATLAB ginput to encounter specific FSDA interaction needs.
% Currently ginputFS is used by just two plotting functions, malfwdplot and
% resfwdplot, to allow user to select interactively a step of the Forward
% Search in correspondence to which the trajectories of units in the subset
% will be highlighted. 
%
% We use ginputFS to get a single point from the current axes, but we added
% as an optional argument a pointerstyle variable which by default is
% 'fullcrosshair', coherently with the style used by ginput, but which can
% be any of the following:
%         crosshair | {arrow} | watch | topl |
%         topr | botl | botr | circle | cross |
%         fleur | left | right | top | bottom |
%         fullcrosshair | ibeam | custom | hand
% Implicit calls made to ginputFS by the current FSDA functions use the
% 'right' style.
%
% See also: ginput
%
% ginput was modified into ginputFS by FSDA team
% The modified code segments have been marked with string FSDAmodif ... FSDAmodifEnd.
%
% Copyright 2008-2015.
% Written by FSDA team
%
% Last modified 06-Feb-2015
%
%Examples
%
%{
% Interactive_example

disp('ginputFS is called inside resfwdplot function')
disp('to get the FS step selected by the user in the interactive session')

load('multiple_regression.txt');
y=multiple_regression(:,4);
X=multiple_regression(:,1:3);
yXplot(y,X);
[out]=LXS(y,X,'nsamp',10000);
[out]=FSReda(y,X,out.bs);
out1=out;
out1.RES=out.RES.^2;
datatooltip = struct;
datatooltip.SubsetLinesColor = FSColors.purplish.RGB;
resfwdplot(out,'datatooltip',datatooltip);

%}

%% Beginning of code

% FSDAmodif
if nargin == 1
    pointerstyle = 'fullcrosshair';
else
    pointerstyle = varargin{1};
end
% FSDAmodifEnd

out1 = []; out2 = []; out3 = []; y = [];
c = computer;
if ~strcmp(c(1:2),'PC')
    tp = get(0,'TerminalProtocol');
else
    tp = 'micro';
end

if ~strcmp(tp,'none') && ~strcmp(tp,'x') && ~strcmp(tp,'micro'),
    if nargout == 1,
        if nargin == 1,
            out1 = trmginput(arg1);
        else
            out1 = trmginput;
        end
    elseif nargout == 2 || nargout == 0,
        if nargin == 1,
            [out1,out2] = trmginput(arg1);
        else
            [out1,out2] = trmginput;
        end
        if  nargout == 0
            out1 = [ out1 out2 ];
        end
    elseif nargout == 3,
        if nargin == 1,
            [out1,out2,out3] = trmginput(arg1);
        else
            [out1,out2,out3] = trmginput;
        end
    end
else
    
    fig = gcf;
    figure(gcf);
    
    if nargin == 0
        how_many = -1;
        b = [];
    else
        how_many = arg1;
        b = [];
        if  ischar(how_many) ...
                || size(how_many,1) ~= 1 || size(how_many,2) ~= 1 ...
                || ~(fix(how_many) == how_many) ...
                || how_many < 0
            error('FSDA:ginputFS:NeedPositiveInt', 'Requires a positive integer.')
        end
        if how_many == 0
        % FSDAmodif
            ptr_fig = 0;
            while(ptr_fig ~= fig)
                ptr_fig = get(0,'PointerWindow');
            end
            scrn_pt = get(0,'PointerLocation');
            loc = get(fig,'Position');
            pt = [scrn_pt(1) - loc(1), scrn_pt(2) - loc(2)];
            out1 = pt(1); y = pt(2);
        elseif how_many < 0
            error('FSDA:ginputFS:InvalidArgument', 'Argument must be a positive integer.')
        end
        % FSDAmodifEnd
    end
    
    % Suspend figure functions
    state = uisuspend(fig);
    
    toolbar = findobj(allchild(fig),'flat','Type','uitoolbar');
    if ~isempty(toolbar)
        ptButtons = [uigettool(toolbar,'Plottools.PlottoolsOff'), ...
            uigettool(toolbar,'Plottools.PlottoolsOn')];
        ptState = get (ptButtons,'Enable');
        set (ptButtons,'Enable','off');
    end
    
    % FSDAmodif
    try
        exception = {};
        set(fig,'pointer',pointerstyle);
    catch exception
    end
    if ~isempty(exception)
        set(fig,'pointer','fullcrosshair');
        warning('FSDA:ginputFS:WrongPointer','Wrong pointer style selected. Set to fullcrosshair');
    end
    
    fig_units = get(fig,'Units');
    char = 0;
    % FSDAmodifEnd
    
    % We need to pump the event queue on unix
    % before calling WAITFORBUTTONPRESS
    drawnow
    
    while how_many ~= 0
        % Use no-side effect WAITFORBUTTONPRESS
        waserr = 0;
        try
            keydown = wfbp ; 
        catch %#ok<CTCH>
            waserr = 1;
        end
        if(waserr == 1)
            if(ishghandle(fig))
                set(fig,'Units',fig_units);
                uirestore(state);
                error('FSDA:ginputFS:Interrupted', 'Interrupted');
            else
                % FSDAmodif
                disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
                disp('The figure you are interacting with has been deleted: PROCESS ENDED');
                disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
                out3 = 99;
                return;
                % FSDAmodifEnd
            end
        end
        % g467403 - ginput failed to discern clicks/keypresses on the figure it was
        % registered to operate on and any other open figures whose handle
        % visibility were set to off
        figchildren = allchild(0);
        if ~isempty(figchildren)
            ptr_fig = figchildren(1);
        else
            error('FSDA:ginputFS:FigureUnavailable','No figure available to process a mouse/key event');
        end
        %         old code -> ptr_fig = get(0,'CurrentFigure'); Fails when the
        %         clicked figure has handlevisibility set to callback
        if(ptr_fig == fig)
            if keydown
                char = get(fig, 'CurrentCharacter');
                button = abs(get(fig, 'CurrentCharacter'));
                scrn_pt = get(0, 'PointerLocation');
                set(fig,'Units','pixels');
                loc = get(fig, 'Position');
                % We need to compensate for an off-by-one error:
                pt = [scrn_pt(1) - loc(1) + 1, scrn_pt(2) - loc(2) + 1];
                set(fig,'CurrentPoint',pt);
            else
                button = get(fig, 'SelectionType');
                if strcmp(button,'open')
                    button = 1;
                elseif strcmp(button,'normal')
                    button = 1;
                elseif strcmp(button,'extend')
                    button = 2;
                elseif strcmp(button,'alt')
                    button = 3;
                else
                    error('FSDA:ginpuFS:InvalidSelection', 'Invalid mouse selection.')
                end
            end
            pt = get(gca, 'CurrentPoint');
            
            how_many = how_many - 1;
            
            if(char == 13) % & how_many ~= 0)
                % if the return key was pressed, char will == 13,
                % and that's our signal to break out of here whether
                % or not we have collected all the requested data
                % points.
                % If this was an early breakout, don't include
                % the <Return> key info in the return arrays.
                % We will no longer count it if it's the last input.
                break;
            end
            
            out1 = [out1;pt(1,1)]; %#ok<AGROW>
            y = [y;pt(1,2)]; %#ok<AGROW>
            b = [b;button]; %#ok<AGROW>
        end
    end
    
    uirestore(state);
    if ~isempty(toolbar) && ~isempty(ptButtons)
        set (ptButtons(1),'Enable',ptState{1});
        set (ptButtons(2),'Enable',ptState{2});
    end
    set(fig,'Units',fig_units);
    
    if nargout > 1
        out2 = y;
        if nargout > 2
            out3 = b;
        end
    else
        out1 = [out1 y];
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = wfbp
%WFBP   Replacement for WAITFORBUTTONPRESS that has no side effects.

fig = gcf;
current_char = []; %#ok<NASGU>

% Now wait for that buttonpress, and check for error conditions
waserr = 0;
try
    h=findall(fig,'Type','uimenu','Accelerator','C');   % Disabling ^C for edit menu so the only ^C is for
    set(h,'Accelerator','');                            % interrupting the function.
    keydown = waitforbuttonpress;
    current_char = double(get(fig,'CurrentCharacter')); % Capturing the character.
    if~isempty(current_char) && (keydown == 1)          % If the character was generated by the
        if(current_char == 3)                           % current keypress AND is ^C, set 'waserr'to 1
            waserr = 1;                                 % so that it errors out.
        end
    end
    
    set(h,'Accelerator','C');                           % Set back the accelerator for edit menu.
catch %#ok<CTCH>
    waserr = 1;
end
drawnow;
if(waserr == 1)
    set(h,'Accelerator','C');                          % Set back the accelerator if it errored out.
    error('MATLAB:ginputFS:Interrupted', 'Interrupted');
end

if nargout>0, key = keydown; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

