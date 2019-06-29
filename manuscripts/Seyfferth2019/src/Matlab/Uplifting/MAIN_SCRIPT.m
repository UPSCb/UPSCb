% MAIN SCRIPT TO ANALYZE SINGLE EXPERIMENT
% 2018-07-06 
%% result5 struct:
% Contains data for all trees:
% 
%   EXP5_T89_1_1_28_H2O: [1×1 struct]
%   EXP5_T89_2_1_28_H2O: [1×1 struct]
%   EXP5_T89_5_1_28_H2O: [1×1 struct]
%   EXP5_T89_6_1_28_H2O: [1×1 struct]
%   EXP5_T89_7_1_28_H2O: [1×1 struct]
%   EXP5_T89_9_1_28_H2O: [1×1 struct]
%   EXP5_etr1L6_2_1_28_H2O: [1×1 struct]
%   EXP5_etr1L6_3_1_28_H2O: [1×1 struct]
%   EXP5_etr1L6_4_1_28_H2O: [1×1 struct]
%   EXP5_etr1L6_5_1_28_H2O: [1×1 struct]
%   EXP5_etr1L6_6_1_28_H2O: [1×1 struct]
%   EXP5_etr1L6_7_1_28_H2O: [1×1 struct]
% 
% Each tree has 3 fields: output, out, vec, example:
% >> result.EXP4_139L4_3_1_28_H2O
%   struct with fields: 
%     output: [28×4 table]- contains curvature and lift for each day:
%     time        dir        curvature     lift  
%     ____    ___________    _________    _______
% 
%       0     [1x32 char]      11.144     -6.5613
%       1     [1x32 char]      4.8036     -6.7875
%       2     [1x32 char]     -4.3793     -6.6857
%       3     [1x32 char]     -5.0444     -6.3869
%       ...
%        out: {14140×3 cell} - contains XY position for each measuring
%        point (several during one day)
%        vec: [14084×5 table] - contains difference between x,y and angle
%        at each measuring time points.

%% meta struct 
% Contains meta information about each tree line:
%          FolderName             genotype       tree     data     treatment     hook     length
%     ________________________    _____________    ____    ______    _________    ______    ______
% 
%     'EXP5-T89_1_1-28_H2O'       'EXP5-T89'       '1'     '1-28'      'H2O'      69.229    121.24
%     'EXP5-T89_2_1-28_H2O'       'EXP5-T89'       '2'     '1-28'      'H2O'       68.21    113.24
%     'EXP5-T89_5_1-28_H2O'       'EXP5-T89'       '5'     '1-28'      'H2O'      71.679    120.11
%     'EXP5-T89_6_1-28_H2O'       'EXP5-T89'       '6'     '1-28'      'H2O'       69.51    118.03
%     'EXP5-T89_7_1-28_H2O'       'EXP5-T89'       '7'     '1-28'      'H2O'      61.191    112.93
%     'EXP5-T89_9_1-28_H2O'       'EXP5-T89'       '9'     '1-28'      'H2O'      79.294    125.16
%     'EXP5-etr1L6_2_1-28_H2O'    'EXP5-etr1L6'    '2'     '1-28'      'H2O'      64.329    109.69
%     'EXP5-etr1L6_3_1-28_H2O'    'EXP5-etr1L6'    '3'     '1-28'      'H2O'      63.171    122.59
%     'EXP5-etr1L6_4_1-28_H2O'    'EXP5-etr1L6'    '4'     '1-28'      'H2O'      70.974    124.06
%     'EXP5-etr1L6_5_1-28_H2O'    'EXP5-etr1L6'    '5'     '1-28'      'H2O'      69.494    129.24
%     'EXP5-etr1L6_6_1-28_H2O'    'EXP5-etr1L6'    '6'     '1-28'      'H2O'      71.888     144.4
%     'EXP5-etr1L6_7_1-28_H2O'    'EXP5-etr1L6'    '7'     '1-28'      'H2O'      67.168    137.93

%%
clear all
close all
clc
fprintf('Analyze EXP5 \n')
% Files result5 and meta5 are created by StemCurve_EXP5_CreateSU.m
StemCurve_EXP5_CreateSU()
% Load already calculated:
load result5
load meta5

Curv5 = calculatecurv_final(meta, result,40); 
% Save result as Excell File:
writetable(Curv5, 'Results\Curv5.xls')
