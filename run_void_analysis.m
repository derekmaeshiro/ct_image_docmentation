% run_void_analysis.m
clear; clc;
type = 'Loose';      % change to 'OtherType', etc as needed
load([type, '_packing.mat']);   % loads BW variable from the file

void_shape_analysis(BW);  % invokes your analysis on the data