clear all; clc; close all;

load('testEnvironment.mat');

answers = f_PER(candSet, problem, W, TXbits, MCS, channelHandles, arrayHandle);