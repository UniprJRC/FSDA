classdef FSColors
%FScolors is a MATLAB class encapsulating some common color definitions  used in the FSDA toolbox. 
%
%  <a href="matlab: docsearchFS('fscolors')">Link to the help function</a>
%
%The associated method returns a cell formed by 
%the RGB vector of the color and a short name for it.
% The colors currently addressed by the method include the standard MATLAB
% colors, i.e.:
%   yellow, magenta, cyan, red, green, blue, white, black
% and the following newly defined colors:
%   blueish, reddish, greenish, purplish, yellowish, greysh, lightblue  
%
% REMARK: This represents a first test in view of making a more extensive
% use of MATLAB classes/methods constructs in FSDA.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%$LastChangedDate::                      $: Date of the last commit
%
%  <a href="matlab: docsearchFS('fscolors')">Link to the help function</a>
%
% Examples:
%
%{
% Usage for, e.g., reddish:
%
RGB_vector = FSColors.reddish.RGB
short_name = FSColors.reddish.ShortName
%
%}

%% Beginning of code
   properties
      RGB = [];
      ShortName = '';
   end
   methods(Static)
       % black and white
       function c = white
          c.RGB =  [1, 1, 1];
          c.ShortName = 'w';
       end
       function c = black
          c.RGB =  [0, 0, 0];
          c.ShortName = 'k';
       end
       
       % standard colors
       function c = yellow
          c.RGB =  [1, 1, 0];
          c.ShortName = 'y';
       end
       function c = magenta
          c.RGB =  [1, 0, 1];
          c.ShortName = 'm';
       end
       function c = cyan
          c.RGB =  [0, 1, 1];
          c.ShortName = 'c';
       end
       function c = red
          c.RGB =  [1, 0, 0];
          c.ShortName = 'r';
       end
       function c = green
          c.RGB =  [0, 1, 0];
          c.ShortName = 'g';
       end
       function c = blue
          c.RGB =  [0, 0, 1];
          c.ShortName = 'b';
       end
       function c = brown
          c.RGB =  [172/255, 115/255, 57/255];
          c.ShortName = 'br';
       end

       % faint colors
       function c = blueish
          c.RGB =  [18/255,104/255,179/255];
          c.ShortName = 'bs';
       end
       function c = reddish
          c.RGB =  [237/255,36/255,38/255];
          c.ShortName = 'rs';
       end
       function c = greenish
          c.RGB =  [155/255,190/255,61/255];
          c.ShortName = 'gs';
       end
       function c = purplish
          c.RGB =  [123/255,45/255,116/255];
          c.ShortName = 'ps';
       end
       function c = yellowish
          c.RGB =  [1,199/255,0];
          c.ShortName = 'ys';
       end
       function c = greysh
          c.RGB =  [0.9, 0.9, 0.9];
          c.ShortName = 'gs';
       end
       function c = lightblue
          c.RGB =  [77/255,190/255,238/255];
          c.ShortName = 'ls';
       end
       function c = lightbrown
          c.RGB =  [223/255, 191/255, 15/255];
          c.ShortName = 'lbr';
       end
       
       % darker colors
       function c = darkblue
          c.RGB =  [0, 0, 128/255];
          c.ShortName = 'db';
       end
       function c = darkred
          c.RGB =  [180/255,0,0];
          c.ShortName = 'dr';
       end
       function c = darkgreen
          c.RGB =  [0,77/255,0];
          c.ShortName = 'dg';
       end
       function c = darkpurpl
          c.RGB =  [100/255,0,130/255];
          c.ShortName = 'dp';
       end
       function c = darkyellow
          c.RGB =  [204/255,204/255,0];
          c.ShortName = 'dy';
       end
       function c = darkgrey
          c.RGB =  [0.6, 0.6, 0.6];
          c.ShortName = 'dg';
       end
       function c = darkbrown
          c.RGB =  [115/255, 77/255, 38/255];
          c.ShortName = 'dbr';
       end
       function c = darkcyan
          c.RGB =  [0, 139/255, 139/255];
          c.ShortName = 'dc';
       end

   end
   
%% REMARK: 
%  From MATLAB Release 2010b classes can contain an enumeration block
%  defining enumeration members. If you have this or later releses, you can
%  comment the static methods definitions above and uncomment the methods
%  and enumeration definitions below. This would be more elegant to read
%  and efficient to execute.
%
%    methods
%       function c = FSColors(r, g, b, sn)
%          c.RGB = [r g b];
%          c.ShortName = sn;
%       end
%    end
%    enumeration
%       yellow (1, 1, 0, 'y')  
%       magenta (1, 0, 1, 'm')  
%       cyan (0, 1, 1, 'c')  
%       red (1, 0, 0, 'r')  
%       green (0, 1, 0, 'g')  
%       blue (0, 0, 1, 'b')  
%       white (1, 1, 1, 'w')  
%       black (0, 0, 0, 'k')  
%       blueish   (18/255,104/255,179/255,'bs')
%       reddish   (237/255,36/255,38/255,'rs')
%       greenish  (155/255,190/255,61/255,'gs')
%       purplish  (123/255,45/255,116/255,'ps')
%       yellowish (1,199/255,0,'ys')
%       greysh    (0.9, 0.9, 0.9, 'gs')
%       lightblue (77/255,190/255,238/255,'lb')
%    end

end
