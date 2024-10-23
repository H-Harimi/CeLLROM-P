% function [gradphie1_tf,D,res0,cellData] = tf_gradphie1(s,locs,cellData) 
%
%   Returns the frequency response of the Nabla [Phi_e]_1(locs,s)/I_app(s) 
%   transfer function as well as any "D" term and integrator residue.
%   (Used in chapter 7 when computing heat-generation terms.)
%
%   Inputs:
%     s: frequency samples to where TF is to be evaluated
%     locs: non-normalized cellData locations in [m]: 0 = negative-electrode 
%           current collector, and so forth
%     cellData: structure with battery parameters
%
%   Outputs:
%     gradphie1_tf: (complex) freq response at requested locations, freqs
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

function [gradphie1_tf,Dterm,res0,cellData] = tf_gradphie1(s,locs,cellData)
  % Create cellData structure field names if they do not exist
  if ~isfield(cellData,'tf')
    cellData.tf = [];
    cellData.tf.name = {};
    cellData.tf.val = [];
    cellData.tf.Dstr = {};
  end
  
  % Copy some variables from cellData structure into more convenient local
  % variable definitions
  Lneg = cellData.neg.L; Lsep = cellData.sep.L; Lpos = cellData.pos.L;
  Ltot = Lneg + Lsep  + Lpos;

  % Check parameters passed to the function
  if(min(locs) < 0 || max(locs) > Ltot+eps)
    error('ERROR (tf_gradphie1): "x" must be in range 0 to %2.2e',Ltot);
  end

  % Set up constants and calculations needed by transfer function
  F = 96485.3365;                 % Faraday constant, [Coulomb/mol]
  R = 8.3144621;                  % Gas constant, [J/mol-K]
  T = cellData.const.T;           % Cell temperature, [K]
  
  Rs_neg = cellData.neg.Rs;       % Particle radius, [m]
  Rs_pos = cellData.pos.Rs;       % Particle radius, [m]
  Ds_neg = cellData.neg.Ds;       % Solid diffusivity [m^2/s]
  Ds_pos = cellData.pos.Ds;       % Solid diffusivity [m^2/s]
  Acell = cellData.const.Acell;   % Current-collector area [m^2]
  as_neg = 3*cellData.neg.eps_s/Rs_neg; % Specific interfacial surf. area
  as_pos = 3*cellData.pos.eps_s/Rs_pos; % Specific interfacial surf. area

  % Effective conductivities of electrolyte and solid
  kappa_eff_neg = cellData.const.kappa*...
                 (cellData.neg.eps_e)^(cellData.neg.brug_kappa);
  kappa_eff_sep = cellData.const.kappa*...
                 (cellData.sep.eps_e)^(cellData.sep.brug_kappa);
  kappa_eff_pos = cellData.const.kappa*...
                 (cellData.pos.eps_e)^(cellData.pos.brug_kappa);
  sigma_eff_neg = cellData.neg.sigma*...
                 (cellData.neg.eps_s)^cellData.neg.brug_sigma;
  sigma_eff_pos = cellData.pos.sigma*...
                 (cellData.pos.eps_s)^cellData.pos.brug_sigma;

  % Electrode state of charge for each electrode at this cell SOC               
  theta_neg = cellData.const.init_SOC * (cellData.neg.theta100 - ...
              cellData.neg.theta0) + cellData.neg.theta0;
  theta_pos = cellData.const.init_SOC * (cellData.pos.theta100 - ...
              cellData.pos.theta0) + cellData.pos.theta0;

  % Copy some variables from cellData structure into more convenient local
  % variable definitions
  csmax_neg = cellData.neg.csmax;
  csmax_pos = cellData.pos.csmax;
  alpha_neg = cellData.neg.alpha;
  alpha_pos = cellData.pos.alpha;
  ce0 = cellData.const.ce0;
  cs0_neg = csmax_neg * theta_neg; % initial cs in negative
  cs0_pos = csmax_pos * theta_pos; % initial cs in positive
  k_neg = cellData.neg.k_norm/csmax_neg/ce0^(1-alpha_neg); % reaction rate
  k_pos = cellData.pos.k_norm/csmax_pos/ce0^(1-alpha_pos); % ... constants
  % Exchange flux densities at this SOC setpoint
  j0_neg = k_neg * (ce0 * (csmax_neg - cs0_neg))^(1 - alpha_neg) ...
                 * cs0_neg^(alpha_neg);
  j0_pos = k_pos * (ce0 * (csmax_pos - cs0_pos))^(1 - alpha_pos) ...
                 * cs0_pos^(alpha_pos);
  Rct_neg = R*T/(j0_neg*F^2); % Charge-transfer resistance 
  Rct_pos = R*T/(j0_pos*F^2); 
  Rfilm_neg = cellData.neg.Rfilm; % Film resistance at particle surface
  Rfilm_pos = cellData.pos.Rfilm; 
  Rtot_neg = Rct_neg + Rfilm_neg; % Total resistance across interface
  Rtot_pos = Rct_pos + Rfilm_pos; 

  % Compute (scaled) transfer function of C_{s,e}(s)/J(s) [Jacobsen-West]
  beta_neg = Rs_neg*sqrt(s./Ds_neg); 
  beta_pos = Rs_pos*sqrt(s./Ds_pos);
  cse_j_neg = Rs_neg/Ds_neg./(1-beta_neg.*coth(beta_neg));
  cse_j_pos = Rs_pos/Ds_pos./(1-beta_pos.*coth(beta_pos));
  % Evaluate derivative of Uocp of both electrodes at given cell SOC
  dudc_neg = cellData.neg.Uocp{2}(theta_neg)/csmax_neg;
  dudc_pos = cellData.pos.Uocp{2}(theta_pos)/csmax_pos;
  % Compute unitless impedance ratio terms for both electrodes
  nu_neg = Lneg*sqrt((as_neg/sigma_eff_neg + as_neg/kappa_eff_neg)...
                   ./(Rtot_neg + dudc_neg.*cse_j_neg/F));
  nu_pos = Lpos*sqrt((as_pos/sigma_eff_pos + as_pos/kappa_eff_pos)...
                   ./(Rtot_pos + dudc_pos.*cse_j_pos/F));
  nu_neg_inf = Lneg * sqrt((as_neg/kappa_eff_neg + ...
                            as_neg/sigma_eff_neg)/Rtot_neg);
  nu_pos_inf = Lpos * sqrt((as_pos/kappa_eff_pos + ...
                            as_pos/sigma_eff_pos)/Rtot_pos);


  len = length(locs); locs = locs(:); s = s(:);
  gradphie1_tf = zeros(len,length(s)); % Initialize output to zeros
  % Initialize "D" variables and transfer-function names variable
  Dterm = zeros(len,1); Dstr = cell(len,1); names = cell(len,1);
  for n = 1:len, % Do for every requested location
    x = locs(n);
    % First case is for "x" in negative electrode
    if(x <= Lneg + eps),
      if (abs(x) < eps), % Special case, x=0, -- avoid numeric issues
        gradphie1_tf(n,:) = 0*s; % TF = 0 @ x = 0
        Dterm(n) = 0;
        Dstr{n} = '0';
      else
        gradphie1_tf(n,:) = (kappa_eff_neg*(sinh((Lneg-x)*nu_neg/Lneg)...
                    - sinh(nu_neg)) - sigma_eff_neg*sinh(x*nu_neg/...
                    Lneg))./(Acell*kappa_eff_neg*(kappa_eff_neg + ...
                    sigma_eff_neg)*sinh(nu_neg));
        tf0 = -x/(Acell*kappa_eff_neg*Lneg); % value at s == 0
        gradphie1_tf(n,s==0) = tf0(:,ones(size(find(s==0))));
        Dterm(n) = (kappa_eff_neg*(sinh((Lneg-x)*nu_neg_inf/Lneg)...
                    - sinh(nu_neg_inf)) - sigma_eff_neg*sinh(x*nu_neg_inf/...
                    Lneg))/(Acell*kappa_eff_neg*(kappa_eff_neg + ...
                    sigma_eff_neg)*sinh(nu_neg_inf));
        Dstr{n} = sprintf(['(kappa_eff_neg*(sinh(nu_neg*(%g))'...
                    '- sinh(nu_neg)) - sigma_eff_neg*sinh(nu_neg*(%g)))'...
                    '/(%g*kappa_eff_neg*(kappa_eff_neg +' ...
                    'sigma_eff_neg)*sinh(nu_neg))'],...
                    (Lneg-x)/Lneg,x/Lneg,Acell);
      end
    % second case is for "x" in positive electrode
    elseif (x >= Lneg + Lsep - eps)
      gradphie1_tf(n,:) = (kappa_eff_pos*(sinh(nu_pos*(x-Lneg-Lsep)/Lpos) ...
             - sinh(nu_pos)) - sigma_eff_pos*sinh(nu_pos*(Ltot-x)/Lpos)) ...
             ./ (Acell*kappa_eff_pos*(kappa_eff_pos+sigma_eff_pos)*...
             sinh(nu_pos));
      tf0 = (x - Ltot)/(Acell*Lpos*kappa_eff_pos); % value at s == 0
      gradphie1_tf(n,s==0) = tf0(:,ones(size(find(s==0))));
      Dterm(n) = (kappa_eff_pos*(sinh(nu_pos_inf*(x-Lneg-Lsep)/Lpos) ...
             - sinh(nu_pos_inf)) - sigma_eff_pos*sinh(nu_pos_inf*(Ltot-x)/Lpos)) ...
             / (Acell*kappa_eff_pos*(kappa_eff_pos+sigma_eff_pos)*...
             sinh(nu_pos_inf));
      Dstr{n} = sprintf(['(kappa_eff_pos*(sinh(nu_pos*(%g))' ...
             '- sinh(nu_pos)) - sigma_eff_pos*sinh(nu_pos*(%g)))' ...
             '/ (%g*kappa_eff_pos*(kappa_eff_pos+sigma_eff_pos)*'...
             'sinh(nu_pos))'],(x-Lneg-Lsep)/Lpos,(Ltot-x)/Lpos,Acell);
    % third case is for "x" in separator
    else
      gradphie1_tf(n,:) = -1/(Acell*kappa_eff_sep);
      Dterm(n) = -1/(Acell*kappa_eff_sep);
      Dstr{n} = sprintf('-1/(%g*kappa_eff_sep)',Acell);
    end
    names{n} = 'gradphie1';
  end % for loop

  res0 = zeros(len,1); % No integrator residue for this transfer function
  cellData.tf.name = [cellData.tf.name; names];  % Store tf name in output
  cellData.tf.val  = [cellData.tf.val; locs(:)]; % Store tf locations
  cellData.tf.Dstr = [cellData.tf.Dstr; Dstr];   % Store D functions  
end % of function gradphie1_tf.m