%hasLOSCluster indicates whether or not a delay profile has an LOS cluster

%   Copyright 2017 The MathWorks, Inc.

%#codegen

function has = hasLOSCluster(DelayProfile,HasLOSCluster)

    if (strcmpi(DelayProfile,'Custom'))
        has = HasLOSCluster;
    else
        has = lte.internal.hasLOSPath(DelayProfile);
    end

end
