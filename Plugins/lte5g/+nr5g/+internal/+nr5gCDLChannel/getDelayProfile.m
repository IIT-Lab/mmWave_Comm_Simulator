%getDelayProfile get delay profile

%   Copyright 2017 The MathWorks, Inc.

%#codegen

function [pathDelays_out,pathGains_out,AoD_out,AoA_out,ZoD_out,ZoA_out] = getDelayProfile(DelayProfile,PathDelays,AveragePathGains,AnglesAoD,AnglesAoA,AnglesZoD,AnglesZoA,HasLOSCluster,KFactorFirstCluster,KFactorScaling,KFactor,DelaySpread)

    coder.extrinsic('lte.internal.getCDLProfile');
    coder.extrinsic('lte.internal.scaleDelaysAndKFactor');

    if(strcmpi(DelayProfile,'Custom'))
        if (nr5g.internal.nr5gCDLChannel.hasLOSCluster(DelayProfile,HasLOSCluster))
            % Split the first path into a LOS part and Rayleigh
            % part according to K_1
            K_1dB = KFactorFirstCluster;
            K_1 = 10^(K_1dB/10);
            P_1dB = AveragePathGains(1);
            P_1 = 10^(P_1dB/10);
            pathDelays_out = [PathDelays(1) PathDelays(1) PathDelays(2:end)];
            pathGains_out = [(10*log10(P_1 * K_1 / (1 + K_1)) + [0 -K_1dB]) AveragePathGains(2:end)];
            AoD_out = [AnglesAoD(1) AnglesAoD];
            AoA_out = [AnglesAoA(1) AnglesAoA];
            ZoD_out = [AnglesZoD(1) AnglesZoD];
            ZoA_out = [AnglesZoA(1) AnglesZoA];
        else
            pathDelays_out = PathDelays;
            pathGains_out = AveragePathGains;
            AoD_out = AnglesAoD;
            AoA_out = AnglesAoA;
            ZoD_out = AnglesZoD;
            ZoA_out = AnglesZoA;
        end
    else
        desiredKFactor = NaN;
        if (nr5g.internal.nr5gCDLChannel.hasLOSCluster(DelayProfile,HasLOSCluster))
            if (KFactorScaling)
                desiredKFactor = KFactor;
            end
        end
        pdp = coder.const(double(lte.internal.getCDLProfile(DelayProfile)));
        pdp = coder.const(double(lte.internal.scaleDelaysAndKFactor(pdp,desiredKFactor,DelaySpread)));
        pathDelays_out = pdp(:,1).'; % 1st column is delay
        pathGains_out = pdp(:,2).';  % 2nd column is power
        AoD_out = pdp(:,3).';
        AoA_out = pdp(:,4).';
        ZoD_out = pdp(:,5).';
        ZoA_out = pdp(:,6).';
    end

    % At this point in the code, if a Rician path is present it is
    % split into a LOS part and a Rayleigh part, whether the delay
    % profile was from the standard tables or is custom

end
