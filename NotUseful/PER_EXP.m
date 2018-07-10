clear; clc; close all;
addpath('utilities','-end');
addpath('./Plugins/lte5g');

problem = o_read_input_problem('data/metaproblem_test.dat');
conf = o_read_config('data/config_test.dat');
conf.verbosity = 0;

totPkt = 200;

PER_LCMV = [];
PER_CBF = [];
PER_HEU = [];

problem.nUsers = 4;

expID = 1; % stands for LCMV etc
% expID = 2; % stands for Heuristics

for antenna = [16 32 64]
    problem.N_Antennas = antenna;
    if antenna == 32
        problem.NxPatch = 8;
        problem.NyPatch = 4;
    else
        problem.NxPatch = sqrt(antenna);
        problem.NyPatch = sqrt(antenna);
    end
    problem.MinObjF = 100.*ones(1,problem.nUsers);
    
    [problem,~,~] = f_configuration(conf,problem);
    candSet = 1 : problem.nUsers;
    
    
    
    
    
    
    %[W_LCMV,W_CBF,arrayHandle_old,estObj_LCMV,estObj_CBF] = ...
    %         f_conventionalBF(problem,conf,candSet);
    [sol_found,W_heu,arrayHandle_heu,estObj] = f_heuristics(problem,conf,candSet);
    
    TXbits = 128 * ones(1, problem.nUsers);
    MCS = 1 * ones(1, problem.nUsers);
    %     PER_LCMV = [PER_LCMV f_PER_stats(candSet, problem, W_LCMV, TXbits, MCS, problem.fullChannels, arrayHandle_old, totPkt)];
    %     PER_CBF = [PER_CBF f_PER_stats(candSet, problem, W_CBF, TXbits, MCS, problem.fullChannels, arrayHandle_old, totPkt)];
    
    PER_HEU = [PER_HEU f_PER_stats(candSet, problem, W_heu, TXbits, MCS, problem.fullChannels, arrayHandle_heu, totPkt)];
    if ~sol_found
        disp('Solution not found, results may be unstable!');
    end
    
    fprintf('Solved for %d antenna, PER_HEU = %.3f\n', antenna, PER_HEU(end));
%     fprintf('Solved for %d antenna, PER_LCMV = %.3f, PER_CBF = %.3f\n', antenna, PER_LCMV(end), PER_CBF(end));
    
end

save('PER_EXP_HEU_ANTENNA.mat');