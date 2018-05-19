%getPathDelays get path delays

%   Copyright 2017 The MathWorks, Inc.

%#codegen

function pathDelays = getPathDelays(theStruct)

    objinfo = lte.internal.CDLChannelInfo(theStruct);
    pathDelays = objinfo.PathDelays;

end
