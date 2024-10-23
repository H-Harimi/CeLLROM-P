% function [gradphis_tf,D,res0,cellData] = tf_gradphis(s,locs,cellData,electrode) 
%
%   Returns the frequency response of the Nabla Phi_s(locs,s)/I_app(s) 
%   transfer function as well as any "D" term and integrator residue.
%   (Used in chapter 7 when computing heat-generation terms.)
%
%   Inputs:
%     s: frequency samples to where TF is to be evaluated
%     locs: non-normalized cellData locations in [m]: 0 = negative-electrode 
%           current collector, and so forth
%     cellData: structure with battery parameters
%     electrode: "neg" for negative and "pos" for positive electrode
%
%   Outputs:
%     gradphis_tf: (complex) freq response at requested locations, freqs
%     D: model "D" term at requested locations
%     res0: integrator residue at requested locations
%     cellData: cellData structure, with added transfer-function terms

% Copyright (c) 2015 by Gregory L. Plett of the University of Colorado 
% Colorado Springs (UCCS). This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0.
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L., "Battery
% Management Systems, Volume I, Battery Modeling," Artech House, 2015.

function [gradphis_tf,Dterm,res0,cellData] = tf_gradphis(s,locs,cellData,electrode)
  % Create cellData structure field names if they do not exist
  if ~isfield(cellData,'tf')
    cellData.tf = [];
    cellData.tf.name = {};
    cellData.tf.val = [];
    cellData.tf.Dstr = {};
  end
  
  % Check parameters passed to the function
  if(max(locs)>1 || min(locs)<0)
    error('ERROR (tf_gradphis): All values of "locs" must be in [0,1]');
  end
  if(strcmpi(electrode,'neg'))
    elec = cellData.neg;
  elseif(strcmpi(electrode,'pos'))
    elec = cellData.pos;
  else
    error('ERROR (tf_gradphis): "electrode" must be "neg" or "pos"');
  end

  % Set up constants and calculations needed by transfer function
  F = 96485.3365;                 % Faraday constant, [Coulomb/mol]
  R = 8.3144621;                  % Gas constant, [J/mol-K]
  T = cellData.const.T;           % Cell temperature, [K]
  Rs = elec.Rs;                   % Particle radius, [m]
  Ds = elec.Ds;                   % Solid diffusivity [m^2/s]
  Acell = cellData.const.Acell;  % Current-collector area [m^2]
  as = 3*elec.eps_s/Rs; % Specific interfacial surface area [m^2/m^3]
  % Effective electrolyte and solid conductivities
  kappa_eff = cellData.const.kappa*(elec.eps_e)^(elec.brug_kappa);
  sigma_eff = elec.sigma*(elec.eps_s)^elec.brug_sigma;
  L = elec.L;                     % Electrode length
  % Electrode stoichiometry at this cell SOC
  theta = cellData.const.init_SOC*(elec.theta100-elec.theta0)+elec.theta0;
  csmax = elec.csmax;             % Maximum lithium concentration
  alpha = elec.alpha;             % Charge-transfer coefficient
  ce0 = cellData.const.ce0;      % Equilibrium electrolyte concentration
  cs0 = csmax * theta;            % Equilibrium solid concentration
  k = elec.k_norm/csmax/ce0^(1-alpha); % Reaction-rate constant
  j0 = k*(ce0*(csmax-cs0))^(1-alpha)*cs0^(alpha); % Exchange flux density
  Rct = R*T/(j0*F^2);             % Charge-transfer resistance
  Rfilm = elec.Rfilm;             % Film resistance
  Rtot = Rct + Rfilm;             % Total interfacial resistance
  % Evaluate derivative of Uocp of electrode at given cell SOC
  dudc = elec.Uocp{2}(theta)/csmax; 
  beta = Rs*sqrt(s./Ds);          % Jacobsen-West beta
  % Compute (scaled) transfer function of C_{s,e}(s)/J(s) [Jacobsen-West]  
  cse_j = Rs/Ds./(1-beta.*coth(beta)); 
  % Compute unitless impedance ratio term (nu_inf as s->infinity)
  nu = L * sqrt((as/sigma_eff + as/kappa_eff)./(Rtot + dudc.*cse_j/F));
  nu_inf = L * sqrt((as/kappa_eff + as/sigma_eff)/(Rct + Rfilm));

  % Now, we are ready to actually calculate the transfer function
  len = length(locs); locs = locs(:);
  gradphis_tf = zeros(len,length(s)); % Initialize output to zero
  % Initialize "D" variables and transfer-function names variable
  Dterm = zeros(len,1); Dstr = cell(len,1); names = cell(len,1);
  for n=1:len, % do for every spatial location
    z = locs(n);
    gradphis_tf(n,:) = (sigma_eff*(sinh(z*nu)-sinh(nu)) + ...
                        kappa_eff*sinh((z-1)*nu)) ./ ...
                       (Acell*sigma_eff*(kappa_eff + sigma_eff)*...
                        sinh(nu));
    Dterm(n) = (sigma_eff*(sinh(z*nu_inf)-sinh(nu_inf)) + ...
                        kappa_eff*sinh((z-1)*nu_inf)) ./ ...
                       (Acell*sigma_eff*(kappa_eff + sigma_eff)*...
                        sinh(nu_inf));
    if(strcmpi(elec.name,'neg'))
      Dstr{n} = sprintf(['(sigma_eff_neg*(sinh(%g*nu_neg)-sinh(nu_neg)) +' ...
                         'kappa_eff_neg*sinh((%g-1)*nu_neg)) ./' ...
                         '(%g*sigma_eff_neg*(kappa_eff_neg + sigma_eff_neg)*'...
                         'sinh(nu_neg))'],z,z,Acell);
    else
      Dstr{n} = sprintf(['-(sigma_eff_pos*(sinh(%g*nu_pos)-sinh(nu_pos)) +' ...
                         'kappa_eff_pos*sinh((%g-1)*nu_pos)) ./' ...
                         '(%g*sigma_eff_pos*(kappa_eff_pos + sigma_eff_pos)*'...
                         'sinh(nu_pos))'],z,z,Acell);
    end    
    names{n} = sprintf('gradphis_%s',electrode);
    
    % value of tf at s->0. Solved with Mathematica.
    tf0 = (locs-1)/(Acell*sigma_eff); 
    gradphis_tf(:,s==0) = tf0(:,ones(size(find(s==0))));

    if any(isnan(gradphis_tf)),
      error('ERROR (tf_gradphis): At least one value computes to NaN');
    end    
  end
  
  % The positive-electrode terms are multiplied by -1.
  if(strcmpi(elec.name,'pos')),
    gradphis_tf = -1*gradphis_tf;
    Dterm = -1*Dterm;
  end
  
  res0 = zeros(len,1); % No integrator residue for this transfer function
  cellData.tf.name = [cellData.tf.name; names];  % Store tf name in output
  cellData.tf.val  = [cellData.tf.val; locs(:)]; % Store tf locations
  cellData.tf.Dstr = [cellData.tf.Dstr; Dstr];   % Store D functions
end % end of function 'tf_gradphis'