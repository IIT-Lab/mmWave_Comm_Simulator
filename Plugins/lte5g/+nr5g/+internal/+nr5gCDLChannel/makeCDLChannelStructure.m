%makeCDLChannelStructure make CDL channel structure

%   Copyright 2017 The MathWorks, Inc.

%#codegen

function model = makeCDLChannelStructure(NormalizePathGains,NormalizeChannelOutputs,MaximumDopplerShift,UTDirectionOfTravel,CarrierFrequency,Seed,DelaySpread,SampleDensity,SampleRate,DelayProfile,PathDelays,AveragePathGains,AnglesAoD,AnglesAoA,AnglesZoD,AnglesZoA,HasLOSCluster,KFactorFirstCluster,KFactorScaling,KFactor,AngleScaling,AngleSpreads,XPR,ClusterDelaySpread,NumStrongestClusters,MeanAngles,channelFilterDelay,TransmitAntennaArray,ReceiveAntennaArray)
            
    coder.extrinsic('lte.internal.makeAntennaArray');
    coder.extrinsic('lte.internal.getCDLPerClusterParameters');

    model = struct();
    model.NormalizePathGains = NormalizePathGains;
    model.NormalizeChannelOutputs = NormalizeChannelOutputs;
    model.TransmitAntennaArray = lte.internal.emptyAntennaArray(TransmitAntennaArray.Size,lte.internal.makeElementPattern(TransmitAntennaArray));
    model.ReceiveAntennaArray = lte.internal.emptyAntennaArray(ReceiveAntennaArray.Size,lte.internal.makeElementPattern(ReceiveAntennaArray));
    taa = coder.const(lte.internal.makeAntennaArray(TransmitAntennaArray));
    raa = coder.const(lte.internal.makeAntennaArray(ReceiveAntennaArray));
    model.TransmitAntennaArray.Position = taa.Position;
    model.TransmitAntennaArray.Orientation = taa.Orientation;
    model.TransmitAntennaArray.ElementPositions = taa.ElementPositions;
    model.TransmitAntennaArray.ElementOrientations = taa.ElementOrientations;
    model.ReceiveAntennaArray.Position = raa.Position;
    model.ReceiveAntennaArray.Orientation = raa.Orientation;
    model.ReceiveAntennaArray.ElementPositions = raa.ElementPositions;
    model.ReceiveAntennaArray.ElementOrientations = raa.ElementOrientations;
    model.MaximumDopplerShift = MaximumDopplerShift;
    model.UTDirectionOfTravel = UTDirectionOfTravel;
    model.CarrierFrequency = CarrierFrequency;
    model.Seed = Seed;
    model.DelaySpread = DelaySpread;
    model.SampleDensity = SampleDensity;
    model.SampleRate = SampleRate;

    [pathDelays,pathGains,AoD,AoA,ZoD,ZoA] = nr5g.internal.nr5gCDLChannel.getDelayProfile(DelayProfile,PathDelays,AveragePathGains,AnglesAoD,AnglesAoA,AnglesZoD,AnglesZoA,HasLOSCluster,KFactorFirstCluster,KFactorScaling,KFactor,DelaySpread);
    model.AveragePathGains = pathGains;
    model.PathDelays = pathDelays;
    model.AnglesAoA = AoA;
    model.AnglesAoD = AoD;
    model.AnglesZoA = ZoA;
    model.AnglesZoD = ZoD;
    model.HasLOSCluster = nr5g.internal.nr5gCDLChannel.hasLOSCluster(DelayProfile,HasLOSCluster);

    if (strcmpi(DelayProfile,'Custom') || AngleScaling)
        model.DesiredASD = AngleSpreads(1);
        model.DesiredASA = AngleSpreads(2);
        model.DesiredZSD = AngleSpreads(3);
        model.DesiredZSA = AngleSpreads(4);
    end

    if (strcmpi(DelayProfile,'Custom'))
        model.XPR = XPR;
        model.ClusterDelaySpread = ClusterDelaySpread;
        model.AngleScaling = false;
        model.NumStrongestClusters = NumStrongestClusters;
    else
        model.ClusterDelaySpread = NaN;
        per_cluster = coder.const(lte.internal.getCDLPerClusterParameters(DelayProfile));
        model.AngleScaling = AngleScaling;
        if (AngleScaling)
            model.DesiredMeanAoD = MeanAngles(1);
            model.DesiredMeanAoA = MeanAngles(2);
            model.DesiredMeanZoD = MeanAngles(3);
            model.DesiredMeanZoA = MeanAngles(4);
        else
            model.DesiredASD = per_cluster.C_ASD;
            model.DesiredASA = per_cluster.C_ASA;
            model.DesiredZSD = per_cluster.C_ZSD;
            model.DesiredZSA = per_cluster.C_ZSA;
        end
        model.XPR = per_cluster.XPR;
        model.NumStrongestClusters = 0;
    end

    model.InitPhase = '38.901';
    model.RayCoupling = 'Random';

    model.ChannelFilterDelay = channelFilterDelay;
    
end
