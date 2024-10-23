% function [cse_tf,D,res0,cellData] = tf_cse(s,locs,cellData,electrode) 
%
%   Returns the frequency response of the C_{s,e}(locs,s)/I_app(s) transfer 
%   function as well as any "D" term and integrator residue, for the given
%   electrode.
%
%   Inputs:
%     s: frequency samples to where TF is to be evaluated
%     locs: normalized cellData locations in [0..1]: 0 = current collector, 
%           and 1 = separator boundary
%     cellData: structure with battery parameters
%     electrode: "neg" for negative and "pos" for positive electrode
%
%   Outputs:
%     cse_tf: (complex) frequency response at requested locations, freqs
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

function [cse_tf,Dterm,res0,cellData] = tf_cse(s,locs,cellData,electrode)
  % Create (output) cellData structure field names if they do not exist
  if ~isfield(cellData,'tf')
    cellData.tf = [];
    cellData.tf.name = {};
    cellData.tf.val = [];
    cellData.tf.Dstr = {};
  end
  
  % Check parameters passed to the function
  if(max(locs)>1 || min(locs)<0)
    error('ERROR (tf_cse): All values of "locs" must be in [0,1]');
  end
  if(strcmpi(electrode,'neg'))
    elec = cellData.neg;
  elseif(strcmpi(electrode,'pos'))
    elec = cellData.pos;
  else
    error('ERROR (tf_cse): "electrode" must be "neg" or "pos"');
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
  % Compute unitless impedance ratio term
  nu = L * sqrt((as/sigma_eff + as/kappa_eff)./(Rtot + dudc.*cse_j/F));

  % Now, we are ready to actually calculate the transfer function
  len = length(locs); locs = locs(:);
  cse_tf = zeros(len,length(s)); % Initialize output to zero
  % Initialize "D" variables and transfer-function names variable
  Dterm = zeros(len,1); Dstr = cell(len,1); names = cell(len,1);
  res0 = -3/(Acell*as*F*L*Rs); % Residue of integrator pole
  for n=1:len % loop through all requested electrode locations
    z = locs(n); % This is the present electrode location
    % First compute tf_j
    cse_tf(n,:) = nu./(as*F*L*Acell*(kappa_eff+sigma_eff).*sinh(nu)).*...
                (sigma_eff*cosh(nu.*z) + kappa_eff*cosh(nu.*(z-1)));
    cse_tf(n,:) = cse_tf(n,:).*cse_j;    % Convert tf_j to tf_cse
    cse_tf(n,:) = cse_tf(n,:) - res0./s; % Remove pole at origin
    Dstr{n} = '0';                       % D = 0 for this TF
    names{n} = sprintf('cse_%s',electrode); % TF name
    
    % value of tf at s->0. Solved with Mathematica.
    tf0 = (-6*dudc*Rs*kappa_eff*sigma_eff +...
      5*as*Ds*F*L^2*((2-6*locs+3*locs.^2)*kappa_eff + ...
      (3*locs.^2-1)*sigma_eff))...
      ./ (30*Acell*as*Ds*dudc*F*L*kappa_eff*sigma_eff);
    cse_tf(:,s==0) = tf0(:,ones(size(find(s==0))));

    if any(isnan(cse_tf)), % Sometimes have problems when dudc -> 0
      error('ERROR (tf_cse): At least one value computes to NaN');
    end
  end
  
  % The positive-electrode terms are multiplied by -1.
  if(strcmp(elec.name,'pos')),
      cse_tf = -1*cse_tf;
      res0 = -1*res0;
  end
  
  res0 = res0*ones(len,1); % Every location has same res0
  cellData.tf.name = [cellData.tf.name; names];  % Append names, 
  cellData.tf.val  = [cellData.tf.val; locs(:)]; % locations, and "D" 
  cellData.tf.Dstr = [cellData.tf.Dstr; Dstr];   % terms to output
  
  if max(abs(cse_tf(:,1))) > 10*max(abs(cse_tf(:,2))),
    % Sometimes have this problem if dudc->0.  "Solved" by manually
    % changing dudc so it never approaches zero, or by including more
    % frequency points closer to zero frequency (increasing Tlen).
    warning('tf_cse.m: dUocp/dCs is too close to zero for good results');
    cse_tf(:,1)=cse_tf(:,2)*10;
  end
end % end of function 'tf_cse'