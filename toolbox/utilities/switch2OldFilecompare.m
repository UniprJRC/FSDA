function switch2OldFilecompare(oldfc)
% switch2OldFilecompare allows the user to switch back and forth to  
% the previous version (ver < R2021b) of the menu functionality
% 'compare selected files and folders' that allows to perform reverse lookup
% on the code of the selected line.
%
% Required input arguments:
%
%    oldfc:     Switch to old file compare, boolean. If olfc is true
%               it is possible to use old file comparison tool.


settingsRoot = settings();
textSetting = settingsRoot.comparisons.text.UseNoJava;
if oldfc == true
    % switch to ver < R2021b version
    textSetting.PersonalValue = false;
else
    % switch back to ver >= R2022a version 
    textSetting.PersonalValue = true;
end