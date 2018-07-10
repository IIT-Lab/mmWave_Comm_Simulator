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

PER_LCMV = [];
PER_CBF = [];
PER_HEU = [];

% expID = 1; % stands for LCMV etc
expID = 2; % stands for Heuristics

for n = [2 4 8 16]
    problem.nUsers = n;
    
    problem.N_Antennas = 128;
    problem.NxPatch = 16;
    problem.NyPatch = 8;
    
    problem.MinObjF = 1.*ones(1,problem.nUsers);
    [problem,~,~] = f_configuration(conf,problem);
    candSet = 1 : problem.nUsers;
    PSDULENGTH = psduLength * ones(1, problem.nUsers);
    MCS = mcs * ones(1, problem.nUsers);
    
    if expID == 1
        [W_LCMV,W_CBF,arrayHandle_old,estObj_LCMV,estObj_CBF] = ...
            f_conventionalBF(problem,conf,candSet);
        PER_LCMV = [PER_LCMV f_PER_trd(candSet, problem, W_LCMV, PSDULENGTH, MCS, problem.fullChannels, arrayHandle_old, totPkt)];
        PER_CBF = [PER_CBF f_PER_trd(candSet, problem, W_CBF, PSDULENGTH, MCS, problem.fullChannels, arrayHandle_old, totPkt)];
        fprintf('Solved for %d users, PER_LCMV = %.3f, PER_CBF = %.3f\n', n, PER_LCMV(end), PER_CBF(end));
    else
        [~, ~, arrayHandle_heu, ~] = f_heuristics(problem, conf, candSet);
        nGoodPkt = zeros(1, totPkt);
        parfor npkt = 1 : totPkt
            [sol_found, W_heu, ~, ~] = f_heuristics(problem, conf, candSet);
            % If all good, nGoodPkt = nUser;
            [nGoodPkt(npkt)] = f_PER_heuristics(candSet, problem, W_heu, PSDULENGTH, MCS, problem.fullChannels, arrayHandle_heu); 
        end
        PER_HEU = [PER_HEU sum(nGoodPkt) / n / totPkt];
        fprintf('Solved for %d users, PER_HEU = %.3f\n', n, PER_HEU(end));
    end
end

if expID == 1
    save('USER_PER_EXP_TRD.mat');
else
    save('USER_PER_EXP_HEU.mat');
end

