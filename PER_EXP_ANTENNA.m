clear; clc; close all;
addpath('utilities');
addpath('lte5g');
addpath('Codes');

problem = o_read_input_problem('data/metaproblem_test.dat');
conf = o_read_config('data/config_test.dat');
conf.verbosity = 0;

totPkt = 200;
mcs = 12;
psduLength = 128;

problem.nUsers = 4;
cdlProfile = {'CDL-A', 'CDL-B', 'CDL-C', 'CDL-D', 'CDL-E'};
ant = [16 64 144 256];

PER_LCMV = zeros(5, 4);
PER_CBF = zeros(5, 4);
PER_HEU = zeros(5, 4);

expID = 1; % stands for LCMV etc
% expID = 2; % stands for Heuristics

candSet = 1 : problem.nUsers;
PSDULENGTH = psduLength * ones(1, problem.nUsers);
MCS = mcs * ones(1, problem.nUsers);

for i = 1 : length(cdlProfile)
    for j = 1 : length(ant)
        
        problem.N_Antennas = ant(j);
        problem.delayProfile = cell2mat(cdlProfile(i));
        problem.MinObjF = 1.*ones(1,problem.nUsers);
        
        if ant(j) == 32
            problem.NxPatch = 8;
            problem.NyPatch = 4;
        else
            problem.NxPatch = sqrt(ant(j));
            problem.NyPatch = sqrt(ant(j));
        end
        
        [problem,~,~] = f_configuration(conf,problem);
        
        if expID == 1
            [W_LCMV,W_CBF,arrayHandle_old,estObj_LCMV,estObj_CBF] = ...
                f_conventionalBF(problem,conf,candSet);
            PER_LCMV(i, j) = f_PER_trd(candSet, problem, W_LCMV, PSDULENGTH, MCS, problem.fullChannels, arrayHandle_old, totPkt);
            PER_CBF(i, j) = f_PER_trd(candSet, problem, W_CBF, PSDULENGTH, MCS, problem.fullChannels, arrayHandle_old, totPkt);
            fprintf('Solved for %d antenna, profile = %s, PER_LCMV = %.3f, PER_CBF = %.3f\n', ant(j), cell2mat(cdlProfile(i)), PER_LCMV(i, j), PER_CBF(i, j));
        else
            [sol_found,W_heu,arrayHandle_heu,estObj] = f_heuristics(problem,conf,candSet);
            if ~sol_found
                disp('Solution not found, results may be unstable!');
            end
            PER_HEU(i, j) = f_PER_stats(candSet, problem, W_heu, PSDULENGTH, MCS, problem.fullChannels, arrayHandle_heu, totPkt);
            fprintf('Solved for %d antenna, PER_HEU = %.3f\n', ant(j), PER_HEU(i, j));
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

