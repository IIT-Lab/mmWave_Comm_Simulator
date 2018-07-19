clear; clc; close all; warning ('off','all');

[~, hostmachine] = system('hostname');

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
cdlProfile = {'CDL-A', 'CDL-B', 'CDL-C', 'CDL-D', 'CDL-E'};
nAntenna = 128;
mcsIndex = 1;

nUserList = [2, 4, 8];

PER_LCMV = zeros(length(cdlProfile), length(nUserList));
PER_CBF = zeros(length(cdlProfile), length(nUserList));
PER_HEU = zeros(length(cdlProfile), length(nUserList));

expID = 1;    % stands for LCMV etc
% expID = 2;      % stands for Heuristics

for i = 1 : length(cdlProfile)
    for j = 1 : length(nUserList)
        
        % overwrite configurations...
        problem.nUsers = nUserList(j);
        problem.N_Antennas = nAntenna;
        problem.MinObjF = 1.*ones(1,problem.nUsers);
        
        if nAntenna == 32
            problem.NxPatch = 8;
            problem.NyPatch = 4;
        elseif nAntenna == 128
            problem.NxPatch = 16;
            problem.NyPatch = 8;
        else
            problem.NxPatch = sqrt(nAntenna);
            problem.NyPatch = sqrt(nAntenna);
        end
        conf.DelayProfile = cell2mat(cdlProfile(i));
        
        candSet = 1 : problem.nUsers;
        PSDULENGTH = psduLength * ones(1, problem.nUsers);
        
        [problem,~,~] = f_configuration(conf, problem);
        
        if expID == 1
            % Call Beamforming
            [W_LCMV,W_CBF,arrayHandle_old,estObj_LCMV,estObj_CBF] = ...
                f_conventionalBF(problem,conf,candSet);
            MCS_LCMV = mcsIndex * ones(1, nUserList(j));
            MCS_CBF = mcsIndex * ones(1, nUserList(j));
            % Evaluate PER
            PER_LCMV(i, j) = f_PER_stats(candSet, problem, W_LCMV, PSDULENGTH, MCS_LCMV, problem.fullChannels, arrayHandle_old, totPkt);
            PER_CBF(i, j) = f_PER_stats(candSet, problem, W_CBF, PSDULENGTH, MCS_CBF, problem.fullChannels, arrayHandle_old, totPkt);
            fprintf('Solved for %d user, profile = %s, PER_LCMV = %.3f, PER_CBF = %.3f\n', nUserList(j), cell2mat(cdlProfile(i)), PER_LCMV(i, j), PER_CBF(i, j));
        else
            % Call Beamforming
            [~,W_heu,arrayHandle_heu,~] = f_heuristics(problem,conf,candSet);
            MCS_heu = mcsIndex * ones(1, nUserList(j));
            PER_HEU(i, j) = f_PER_stats(candSet, problem, W_heu, PSDULENGTH, MCS_heu, problem.fullChannels, arrayHandle_heu, totPkt);
            fprintf('Solved for %d user, profile = %s, PER_HEU = %.3f\n', nUserList(j), cell2mat(cdlProfile(i)), PER_HEU(i, j));
        end
    end
    
    if expID == 1
        save(['USER_PER_TRD_' num2str(nAntenna) '_ANT_' num2str(mcsIndex) '_MCS.mat']);
    else
        save(['USER_PER_HEU_' num2str(nAntenna) '_ANT_' num2str(mcsIndex) '_MCS.mat']);
    end
    
end

if expID == 1
    save(['USER_PER_TRD_' num2str(nAntenna) '_ANT_' num2str(mcsIndex) '_MCS.mat']);
else
    save(['USER_PER_HEU_' num2str(nAntenna) '_ANT_' num2str(mcsIndex) '_MCS.mat']);
end

