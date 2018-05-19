function lte5gPostInstallCmd
%LTE5GPOSTINSTALLCMD displays post installation dialog

% Copyright 2016-2018 The MathWorks, Inc.

dialogMessage = getString(message('lte5g:lte5gPostInstallCmd:DialogMessage'));
dialog = msgbox(dialogMessage, 'modal');
uiwait(dialog);

help lte5g;

end