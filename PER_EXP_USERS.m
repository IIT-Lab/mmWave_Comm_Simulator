% COMMENT
% Make sure the following requirement is met:
% ceil(sqrt(problem.N_Antennas)) > problem.nUsers

clear; clc; close all; warning ('off','all');
addpath('../mmWave-MU-MIMO');
addpath('../mmWave-MU-MIMO/utilities');
addpath('../mmWave-MU-MIMO/data');
addpath('../mmWave-MU-MIMO/PER_TABLES');
addpath('./Plugins/');

problem = o_read_input_problem('../mmWave-MU-MIMO/data/metaproblem_test.dat');
problem.DEBUG = false;
conf = o_read_config('../mmWave-MU-MIMO/data/config_test.dat');
conf.verbosity = 0;

totPkt = 200;
psduLength = 128;
cdlProfile = {'CDL-A', 'CDL-B', 'CDL-C', 'CDL-D', 'CDL-E'};
nAntenna = 64;

nUserList = [2, 4, 8];

PER_LCMV = zeros(length(cdlProfile), length(nUserList));
PER_CBF = zeros(length(cdlProfile), length(nUserList));
PER_HEU = zeros(length(cdlProfile), length(nUserList));

expID = 1; % stands for LCMV etc
% expID = 2; % stands for Heuristics


for i = 1 : length(cdlProfile)
    for j = 1 : length(nUserList)
        
        problem.nUsers = nUserList(j);
        problem.N_Antennas = nAntenna;
        problem.MinObjF = 1.*ones(1,problem.nUsers);
        problem.NxPatch = sqrt(nAntenna);
        problem.NyPatch = sqrt(nAntenna);
        conf.DelayProfile = cell2mat(cdlProfile(i));
        
        candSet = 1 : problem.nUsers;
        PSDULENGTH = psduLength * ones(1, problem.nUsers);
        
        [problem,~,~] = f_configuration(conf,problem);
        
        if expID == 1
            % Call Beamforming
            [W_LCMV,W_CBF,arrayHandle_old,estObj_LCMV,estObj_CBF] = ...
                f_conventionalBF(problem,conf,candSet);
            [~,~,~,SNRList_LCMV]  = f_BF_results(W_LCMV,arrayHandle_old,problem,conf,false);
            [~,~,~,SNRList_CBF]  = f_BF_results(W_CBF,arrayHandle_old,problem,conf,false);
            % Select MCS for estimated SNR
            [MCS_LCMV,~] = f_selectMCS(candSet,pow2db(SNRList_LCMV),problem.targetPER,problem.MCSPER,problem.DEBUG);
            [MCS_CBF,~] = f_selectMCS(candSet,pow2db(SNRList_CBF),problem.targetPER,problem.MCSPER,problem.DEBUG);
            % Evaluate PER
            PER_LCMV(i, j) = f_PER_trd(candSet, problem, W_LCMV, PSDULENGTH, MCS_LCMV, problem.fullChannels, arrayHandle_old, totPkt);
            PER_CBF(i, j) = f_PER_trd(candSet, problem, W_CBF, PSDULENGTH, MCS_CBF, problem.fullChannels, arrayHandle_old, totPkt);
            fprintf('Solved for %d user, profile = %s, PER_LCMV = %.3f, PER_CBF = %.3f\n', nUserList(j), cell2mat(cdlProfile(i)), PER_LCMV(i, j), PER_CBF(i, j));
        else
            % Call Beamforming
            [sol_found,W_heu,arrayHandle_heu,estObj] = f_heuristics(problem,conf,candSet);
            [~,~,~,SNRList]  = f_BF_results(W_heu,arrayHandle_heu,problem,conf,false);
            % Select MCS for estimated SNR
            [MCS,~] = f_selectMCS(candSet,SNRList,problem.targetPER,problem.MCSPER,problem.DEBUG);
            if ~sol_found; disp('Solution not found, results may be unstable!'); end
            % Evaluate PER
            PER_HEU(i, j) = f_PER_stats(candSet, problem, W_heu, PSDULENGTH, MCS, problem.fullChannels, arrayHandle_heu, totPkt);
            fprintf('Solved for %d user, PER_HEU = %.3f\n', nUserList(j), PER_HEU(i, j));
        end
    end
    
    if expID == 1
        save('ANTENNA_PER_EXP_TRD.mat');
    else
        save('ANTENNA_PER_EXP_HEU.mat');
    end
    
end

if expID == 1
    save('ANTENNA_PER_EXP_TRD.mat');
else
    save('ANTENNA_PER_EXP_HEU.mat');
end

