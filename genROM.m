% function genROM(cellParams,DRAparams,ROMDIR)
%
%   Generates reduced-order model (ROM) files for the cell described by the
%   parameters in the Excel spreadsheet "cellParams" using the DRA and DRA
%   tuning parameters in the "DRAparams" spreadsheet.  The output ROM files
%   are placed in the directory ROMDIR, which is created if it does not
%   already exist.
%
%   Inputs:
%     cellParams: Filename for Excel spreadsheet holding cell parameters
%     DRAparams: Filename for Excel spreadsheet holding DRA tuning values
%     ROMDIR: Path to directory in which to place output ROM files

% Copyright (c) 2015 by Gregory L. Plett of the University of Colorado 
% Colorado Springs (UCCS). This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0.
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L., "Battery
% Management Systems, Volume I, Battery Modeling," Artech House, 2015.

function genROM(cellParams,DRAparams,ROMDIR)
  reg = {'neg','pos'}; % Names for two electrode regions
  R   = 8.3144621;     % Universal gas constant (used in Arrhenius eq.)
  
  % Load cell parameters
  [cell,RevCell] = readParamTable(cellParams,'Parameters','No');

  % Load DRA configuration parameters
  [draData,tflist,RevDRA] = readDRATable(DRAparams,'Parameters');

  % Compatibility check: This version of is compatible with version 
  % 1.0 of Cell Params and version 1.0 of DRA Params
  if RevCell ~= 1.0,
    error('genROM.m isn''t compatible with this cell param version');
  elseif RevDRA ~= 1.0,
    error('genROM.m isn''t compatible with this DRA param version');
  end

  % Check to see if last character is '/' and if the directory exists
  if ROMDIR(end) ~= '/', ROMDIR = [ROMDIR, '/']; end
  if ~exist(ROMDIR,'dir'), mkdir(ROMDIR); end

  % Run DRA for every temperature in draData.T and for every
  % state-of-charge in draData.SOC
  for tt = 1:length(draData.T), 
    for qq = 1:length(draData.SOC),
      T   = draData.T(tt)+273.15; % Temperature in Kelvin (convert from C)
      SOC = draData.SOC(qq)/100;  % SOC as fraction (convert from %)
      cell.const.init_SOC = SOC;  % Store this SOC setpoint in output ROM
      cell.const.T        = T;    % Store this temperature setpoint too

      % Compute values of parameters that vary with temperature for this
      % particular setpoint using Arrhenius equation
      RTfact = (1/cell.const.Tref - 1/T)/R;
      cell.const.De = cell.const.De_ref * exp(cell.const.Eact_De * RTfact);
      % If the kappa_ref field is a string, then it represents a function
      % (with input parameter c_e). So, evaluate the function at c_{e,0} to
      % get reference kappa value.  Otherwise, kappa_ref is a constant
      if ischar(cell.const.kappa_ref),                            
        kappa = eval(cell.const.kappa_ref);
        cell.const.kappa = kappa(cell.const.ce0) ...
                                 * exp(cell.const.Eact_kappa*RTfact);
      else
        cell.const.kappa = cell.const.kappa_ref ...
                                 * exp(cell.const.Eact_kappa*RTfact);
      end
      for ii = 1:length(reg), % Compute Ds,k_norm,sigma and Uocp
        cell.(reg{ii}).Ds = cell.(reg{ii}).Ds_ref...
                                  * exp(cell.(reg{ii}).Eact_Ds * RTfact);
        cell.(reg{ii}).k_norm = cell.(reg{ii}).k_norm_ref...
                                  * exp(cell.(reg{ii}).Eact_k * RTfact);
        cell.(reg{ii}).sigma = cell.(reg{ii}).sigma_ref...
                                * exp(cell.(reg{ii}).Eact_sigma * RTfact);
        % If the Uocp field is a string, then it represents an inline
        % function definition (with input parameter c_{s,e}). This string
        % is converted into a standard inline function via the "eval" cmd.
        if ischar(cell.(reg{ii}).Uocp),
          cell.(reg{ii}).Uocp_str = cell.(reg{ii}).Uocp;
          cell.(reg{ii}).Uocp = eval(cell.(reg{ii}).Uocp_str);
        end
      end

      % Everything is set up, so run the DRA.
      % First, print the date and time so we know computer is not hung,
      % since the operation can take some time
      fprintf('\nStarting: %s...\n',datestr(now));
      fprintf('T: %2.2fK, SOC: %d%%, Fs: %d, Tlen: %d, Order: [',...%s\n',...
               T,SOC*100,draData.Fs,draData.Tlen);
      fprintf(' %g ',draData.order); fprintf(']\n');               
      fprintf('-----------------------------------------------------\n');
      out = dra(draData,tflist,cell); % Call the DRA routine itself

      numROMs = length(out); % Can produce more than one ROM at once
      for ind = 1:numROMs
        ROM = out(ind);
        % Remove @(x) Uocp functions, but keep string versions
        ROM.cellData.pos = rmfield(ROM.cellData.pos,'Uocp');
        ROM.cellData.neg = rmfield(ROM.cellData.neg,'Uocp');

        outfile = sprintf('%scell_%s_ROM_%dSOC_%dC_%dx_%dx%d.mat',...
                    ROMDIR,cell.name,round(100*SOC),...
                    round(T-273.15),ROM.order,draData.Fs,...
                    max(max(draData.hank1loc),max(draData.hank2loc)));
        save(outfile,'ROM');
      end      
    end
  end
end