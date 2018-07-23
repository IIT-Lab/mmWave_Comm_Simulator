clear; clc; close all; warning ('off','all');

[~, hostmachine] = system('hostname');

% To run this program faster, it is suggested to use computer that has more
% than 48G of memory, and don't open up more than 4 concurrent workers
% (might blow up memory)
% if strcmp(hostmachine, 'Zeta') || strcmp(hostmachine, 'Iota')
%     error('Do not run on other machines except for Zeta and Iota!')
% end

addpath('../mmWave-MU-MIMO');
addpath('../mmWave-MU-MIMO/utilities');
addpath('../mmWave-MU-MIMO/data');
addpath('../mmWave-MU-MIMO/PER_TABLES');
addpath('./Plugins/lte5g/');

problem = o_read_input_problem('../mmWave-MU-MIMO/data/metaproblem_test.dat');
problem.DEBUG = false;
conf = o_read_config('../mmWave-MU-MIMO/data/config_test.dat');
% overwrite configurations...
conf.verbosity = 0;
conf.NumPhaseShifterBits = 0;                                       % Number of bits to control heuristic solution
conf.FunctionTolerance_Data = 1e-10;                                % Heuristics stops when not improving solution by this much
conf.PopSizeList = 30;                                              % Number of genes in heuristics
conf.Maxgenerations_Data = 150;                                     % Maximum generations (iterations)
conf.EliteCount_Data = ceil(conf.PopulationSize_Data/5);            % Number of genes unchanged from one generation to the next
conf.MaxStallgenerations_Data = ceil(conf.Maxgenerations_Data/4);   % Stop if fitness value does not improve

totPkt = 200;
psduLength = 128;
cdlProfile = {'CDL-C'};
nUsers = 4;
mcsIndex = 1;

nAntennaList = [16, 32, 64, 128, 144];

PER_LCMV = zeros(length(cdlProfile), length(nAntennaList));
PER_CBF = zeros(length(cdlProfile), length(nAntennaList));
PER_HEU = zeros(length(cdlProfile), length(nAntennaList));

expID = 1;    % stands for LCMV etc
% expID = 2;      % stands for Heuristics

overwriteLocation = true;

for i = 1 : length(cdlProfile)
    for j = 1 : length(nAntennaList)
        
        % overwrite configurations...
        problem.nUsers = nUsers;
        problem.N_Antennas = nAntennaList(j);
        problem.MinObjF = 1.*ones(1,problem.nUsers);
        
        if nAntennaList(j) == 32
            problem.NxPatch = 8;
            problem.NyPatch = 4;
        elseif nAntennaList(j) == 128
            problem.NxPatch = 16;
            problem.NyPatch = 8;
        else
            problem.NxPatch = sqrt(nAntennaList(j));
            problem.NyPatch = sqrt(nAntennaList(j));
        end
        conf.DelayProfile = cell2mat(cdlProfile(i));
        
        candSet = 1 : nUsers;
        PSDULENGTH = psduLength * ones(1, nUsers);
        
        [problem,~,~] = f_configuration(conf, problem);
        
        if expID == 1
            if overwriteLocation
                problem.phiUsers = randi([-180 180], 1, problem.nUsers);   % Azimuth
                problem.thetaUsers = randi([-90 90], 1, problem.nUsers);   % Elevation
            end
            % Call Beamforming
            [W_LCMV,W_CBF,arrayHandle_old,estObj_LCMV,estObj_CBF] = ...
                f_conventionalBF(problem,conf,candSet);
            MCS_LCMV = mcsIndex * ones(1, nUsers);
%             MCS_CBF = mcsIndex * ones(1, nUsers);
            % Evaluate PER
            PER_LCMV(i, j) = f_PER_stats(candSet, problem, W_LCMV, PSDULENGTH, MCS_LCMV, problem.fullChannels, arrayHandle_old, totPkt);
%             PER_CBF(i, j) = f_PER_stats(candSet, problem, W_CBF, PSDULENGTH, MCS_CBF, problem.fullChannels, arrayHandle_old, totPkt);
            fprintf('Solved for %d antenna, profile = %s, PER_LCMV = %.3f\n', nAntennaList(j), cell2mat(cdlProfile(i)), PER_LCMV(i, j));
        else
            % Call Beamforming
            [~,W_heu,arrayHandle_heu,~] = f_heuristics(problem,conf,candSet);
            MCS_heu = mcsIndex * ones(1, nUsers);
            PER_HEU(i, j) = f_PER_stats(candSet, problem, W_heu, PSDULENGTH, MCS_heu, problem.fullChannels, arrayHandle_heu, totPkt);
            fprintf('Solved for %d antenna, profile = %s, PER_HEU = %.3f\n', nAntennaList(j), cell2mat(cdlProfile(i)), PER_HEU(i, j));
        end
    end
    
    if expID == 1
        save(['ANT_PER_TRD_' num2str(nUsers) '_USERS_' num2str(mcsIndex) '_MCS.mat']);
    else
        save(['ANT_PER_HEU_' num2str(nUsers) '_USERS_' num2str(mcsIndex) '_MCS.mat']);
    end
    
end

if expID == 1
    save(['ANT_PER_TRD_' num2str(nUsers) '_USERS_' num2str(mcsIndex) '_MCS.mat']);
else
    save(['ANT_PER_HEU_' num2str(nUsers) '_USERS_' num2str(mcsIndex) '_MCS.mat']);
end

