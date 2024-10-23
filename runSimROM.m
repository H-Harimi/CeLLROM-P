% script runSimROM.m
%
%   Example script that shows how to call the simROM.m function.

% Copyright (c) 2015 by Gregory L. Plett of the University of Colorado 
% Colorado Springs (UCCS). This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0.
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L., "Battery
% Management Systems, Volume I, Battery Modeling," Artech House, 2015.

clc
clear all
close all

% set the ROMDIR variable to the directory containing the ROMs created by
% the DRA via genROM.m
ROMDIR = './Doyle_ROM_Files_4/'; % four states plus integrator
cellType = 'doyle';              % Cell type: 'doyle' or ...

% Precompute vectors of 15:00 minutes, and 15:01 minutes
z900 = zeros(900,1); z901 = zeros(901,1);

runSims = 1:8; % set to vector of simulation cases to be executed
for sim = runSims,
  switch sim      
    case 1, % charge-balanced UDDS profile, max 2C rate, 60% SOC
      UDDS = load('UDDS_2C_SOCconst.txt');
      deltat = 1; soc0 = 0.6;
      Tvec = 298.15*ones(size(UDDS(:,1)));
      ivec = [0; UDDS(1:end-1,2)]*20.5; % scale to max 2C rate; 1C=20.5A
      blend = 1; % model blending? 0 = no; 1 = yes
      dTerm = 'uniform'; % "D" computed from Dstr every sample
      [vCell,ROM] = simROM(ROMDIR,cellType,ivec,Tvec,soc0,blend,dTerm);
      save doyle_udds_balanced_ROM.mat ROM
    case 2, % rest, 1C discharge, rest, starting at 100% SOC
      Crate = 1; ik = 20.5*Crate; soc0 = 1;
      simLen = 900+3600/Crate+900; 
      Tvec = 298.15*ones(simLen+1,1);
      blend = 1; % model blending? 0 = no; 1 = yes
      dTerm = 'uniform'; % "D" computed from Dstr every sample
      ivec = [z901; ik*ones(3600/Crate,1); z900];
      [vCell,ROM] = simROM(ROMDIR,cellType,ivec,Tvec,soc0,blend,dTerm);
      save doyle_1C_ROM.mat ROM
    case 3, % rest, C/2 discharge, rest, starting at 100% SOC
      Crate = 0.5; ik = 20.5*Crate; soc0 = 1;
      simLen = 900+3600/Crate+900; 
      Tvec = 298.15*ones(simLen+1,1);
      blend = 1; % model blending? 0 = no; 1 = yes
      dTerm = 'linear'; % "D" computed by interpolating DRA D matrix
      ivec = [z901; ik*ones(3600/Crate,1); z900];
      [vCell,ROM] = simROM(ROMDIR,cellType,ivec,Tvec,soc0,blend,dTerm);
      save doyle_0.5C_ROM.mat ROM
    case 4, % rest, C/4 discharge, rest, starting at 100% SOC
      Crate = 0.25; ik = 20.5*Crate; soc0 = 1;
      simLen = 900+3600/Crate+900; 
      Tvec = 298.15*ones(simLen+1,1);
      blend = 1; % model blending? 0 = no; 1 = yes
      dTerm = 'linear'; % "D" computed by interpolating DRA D matrix
      ivec = [z901; ik*ones(3600/Crate,1); z900];
      [vCell,ROM] = simROM(ROMDIR,cellType,ivec,Tvec,soc0,blend,dTerm);
      save doyle_0.25C_ROM.mat ROM
    case 5,  % rest, C/10 discharge, rest, starting at 100% SOC
      Crate = 0.1; ik = 20.5*Crate; soc0 = 1;
      simLen = 900+3600/Crate+900; 
      Tvec = 298.15*ones(simLen+1,1);
      blend = 1; % model blending? 0 = no; 1 = yes
      dTerm = 'linear'; % "D" computed by interpolating DRA D matrix
      ivec = [z901; ik*ones(3600/Crate,1); z900];
      [vCell,ROM] = simROM(ROMDIR,cellType,ivec,Tvec,soc0,blend,dTerm);
      save doyle_0.1C_ROM.mat ROM         
    case 6, % rest, short 1C discharge, rest, starting at 100% SOC
      ivec = [zeros(10,1); 20.5*ones(10,1); zeros(10,1)];
      simLen = length(ivec)-1;
      soc0 = 1;
      Tvec = 298.15*ones(simLen+1,1);
      blend = 1; % model blending? 0 = no; 1 = yes
      dTerm = 'linear'; % "D" computed by interpolating DRA D matrix
      [vCell,ROM] = simROM(ROMDIR,cellType,ivec,Tvec,soc0,blend,dTerm);
      save doyle_dpulse_ROM.mat ROM   
    case 7,  % rest, short 2C charge, rest, starting at 100% SOC
      ivec = -2*[zeros(10,1); 20.5*ones(10,1); zeros(10,1)];
      simLen = length(ivec)-1;
      soc0 = 1;
      Tvec = 298.15*ones(simLen+1,1);
      blend = 1; % model blending? 0 = no; 1 = yes
      dTerm = 'uniform'; % "D" computed from Dstr every sample
      [vCell,ROM] = simROM(ROMDIR,cellType,ivec,Tvec,soc0,blend,dTerm);
      save doyle_cpulse_ROM.mat ROM   
    case 8, % ten consecutive charge-depleting UDDS cycles, starting at ...
      UDDS = load('UDDS_2C_10cycles.txt');  % 80% SOC
      deltat = 1; soc0 = 0.8;
      Tvec = 298.15*ones(size(UDDS(:,1)));
      ivec = [0; UDDS(1:end-1,2)]*20.5;
      blend = 1; % model blending? 0 = no; 1 = yes
      dTerm = 'uniform'; % "D" computed from Dstr every sample
      [vCell,ROM] = simROM(ROMDIR,cellType,ivec,Tvec,soc0,blend,dTerm);
      save doyle_long_UDDS_ROM.mat ROM
  end
  % Plot cell voltage of most recent simulation
  figure(sim); clf;
  plot(0:length(vCell)-1,vCell,'r'); grid on;
  xlabel('Time (s)');
  ylabel('Cell voltage (V)');
  title(sprintf('Doyle cell voltage for case %d',sim)); 
  drawnow;
end