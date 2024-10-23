% script runGenROM_DRA.m
%
%   Example script that shows how to call the genROM_DRA.m function.

% Copyright (c) 2015 by Gregory L. Plett of the University of Colorado 
% Colorado Springs (UCCS). This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0.
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L., "Battery
% Management Systems, Volume I, Battery Modeling," Artech House, 2015.

clc
close all
clear all

% Generate reduced-order model files in the folder named "ROMFiles" based
% on cell parameters in "Doyle_parameter_list.xlsx" and on DRA tuning
% parameters in "Doyle_DRA_list.xlsx".
genROM('Doyle_parameter_list.xlsx','Doyle_DRA_list.xlsx','ROMFiles');



