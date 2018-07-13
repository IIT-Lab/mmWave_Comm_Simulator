clear; clc; close all; warning ('off','all');
addpath('../mmWave-MU-MIMO');
addpath('../mmWave-MU-MIMO/utilities');
addpath('../mmWave-MU-MIMO/data');
addpath('../mmWave-MU-MIMO/PER_TABLES');
addpath('./Plugins/lte5g/');

problem = o_read_input_problem('../mmWave-MU-MIMO/data/metaproblem_test.dat');
problem.DEBUG = false;
conf = o_read_config('../mmWave-MU-MIMO/data/config_test.dat');
% overwrite configurations...
conf.verbosity = 1;
conf.NumPhaseShifterBits = 0;                                       % Number of bits to control heuristic solution
conf.FunctionTolerance_Data = 1e-10;                                % Heuristics stops when not improving solution by this much
conf.PopSizeList = 30;                                              % Number of genes in heuristics
conf.Maxgenerations_Data = 150;                                     % Maximum generations (iterations)
conf.EliteCount_Data = ceil(conf.PopulationSize_Data/5);            % Number of genes unchanged from one generation to the next
conf.MaxStallgenerations_Data = ceil(conf.Maxgenerations_Data/4);   % Stop if fitness value does not improve

totPkt = 100;
psduLength = 128;
cdlProfile = {'CDL-A', 'CDL-B', 'CDL-C', 'CDL-D', 'CDL-E'};
nAntenna = 64;
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
%             [~,~,~,SNRList_LCMV]  = f_BF_results(W_LCMV,arrayHandle_old,problem,conf,false);
%             [~,~,~,SNRList_CBF]  = f_BF_results(W_CBF,arrayHandle_old,problem,conf,false);
%             % Select MCS for estimated SNR
%             [MCS_LCMV,~] = f_selectMCS(candSet,pow2db(SNRList_LCMV),problem.targetPER,problem.MCSPER,problem.DEBUG);
%             [MCS_CBF,~] = f_selectMCS(candSet,pow2db(SNRList_CBF),problem.targetPER,problem.MCSPER,problem.DEBUG);
            MCS_LCMV = mcsIndex * ones(1, nUserList(j));
            MCS_CBF = mcsIndex * ones(1, nUserList(j));
            % Evaluate PER
            PER_LCMV(i, j) = f_PER_trd(candSet, problem, W_LCMV, PSDULENGTH, MCS_LCMV, problem.fullChannels, arrayHandle_old, totPkt);
            PER_CBF(i, j) = f_PER_trd(candSet, problem, W_CBF, PSDULENGTH, MCS_CBF, problem.fullChannels, arrayHandle_old, totPkt);
            fprintf('Solved for %d user, profile = %s, PER_LCMV = %.3f, PER_CBF = %.3f\n', nUserList(j), cell2mat(cdlProfile(i)), PER_LCMV(i, j), PER_CBF(i, j));
        else
            % Call Beamforming
            perTemp = zeros(1, 5);
            for iter = 1 : 5
                [~,W_heu,arrayHandle_heu,~] = f_heuristics(problem,conf,candSet);
                MCS_heu = mcsIndex * ones(1, nUserList(j));
                perTemp(iter) = f_PER_trd(candSet, problem, W_heu, PSDULENGTH, MCS_heu, problem.fullChannels, arrayHandle_heu, totPkt);
            end
            PER_HEU(i, j) = min(perTemp);
%             [~,~,~,SNRList]  = f_BF_results(W_heu,arrayHandle_heu,problem,conf,false);
%             % Select MCS for estimated SNR
%             [MCS,~] = f_selectMCS(candSet,SNRList,problem.targetPER,problem.MCSPER,problem.DEBUG);
%             if ~sol_found; disp('Solution not found, results may be unstable!'); end
%             % Evaluate PER
%             PER_HEU(i, j) = f_PER_stats(candSet, problem, W_heu, PSDULENGTH, MCS, problem.fullChannels, arrayHandle_heu, totPkt);
            fprintf('Solved for %d user, PER_HEU = %.3f\n', nUserList(j), PER_HEU(i, j));
        end
    end
    
    if expID == 1
        save(['USER_PER_EXP_TRD_' num2str(nAntenna) '_MCS_' num2str(mcsIndex) '.mat']);
    else
        save(['USER_PER_EXP_HEU_' num2str(nAntenna) '_MCS_' num2str(mcsIndex) '.mat']);
    end
    
end

if expID == 1
    save(['USER_PER_EXP_TRD_' num2str(nAntenna) '_MCS_' num2str(mcsIndex) '.mat']);
else
    save(['USER_PER_EXP_HEU_' num2str(nAntenna) '_MCS_' num2str(mcsIndex) '.mat']);
end

