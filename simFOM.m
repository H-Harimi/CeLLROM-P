% function [Vcell,FOM,model] = simFOM(cellData,ivec,tvec,initSOC)
%
%   Simulates the full-order model version of the lithium-ion cell defined
%   by the parameters in cellData for the input current in ivec and the
%   time vector in tvec, starting at initial SOC = initSOC.
%
%   Uses COMSOL LiveLink with MATLAB (required add-on software)
%
%   Inputs:
%     cellData: Structure with battery parameters
%     ivec: Vector of cell current versus time [A]
%     tvec: Vector of output sample times [s]
%     initSOC: Cell SOC at beginning of simulation
%
%   Outputs:
%     Vcell: Vector of cell voltage versus time
%     FOM: All the internal cell variables from the simulation
%     model: The COMSOL model, which can be saved to file for more analysis

% Copyright (c) 2015 by Gregory L. Plett of the University of Colorado 
% Colorado Springs (UCCS). This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0.
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L., "Battery
% Management Systems, Volume I, Battery Modeling," Artech House, 2015.

function [Vcell,FOM,model] = simFOM(cellData,ivec,tvec,initSOC)

  % Create a new COMSOL model
  import com.comsol.model.*
  import com.comsol.model.util.*
  try
    model = ModelUtil.create('Model');
  catch
    error('Error creating COMSOL model. Please make sure that "COMSOL with MATLAB" is running.');
  end
  ModelUtil.showProgress(true); % Turn on progress bar
  model.author('Gregory L. Plett');

  ivec = ivec(:); tvec = tvec(:);
  % Calculate initial electrode SOCs from cell SOC
  theta_init_neg = initSOC*(cellData.neg.theta100 - cellData.neg.theta0)+...
                            cellData.neg.theta0;
  theta_init_pos = initSOC*(cellData.pos.theta100 - cellData.pos.theta0)+...
                            cellData.pos.theta0;

  %% Comsol global variables: Filled in with values read from cellData
  model.param.set('UNIVERSAL_PARAMETERS', '0');
  model.param.set('R', '8.314 [J/mol/K]', 'Gas constant');
  model.param.set('F', '96485 [C/mol]', 'Faraday''s constant');
  model.param.set('Tref', '(25 + 273.15) [K]', 'Reference temperature');
  model.param.set('T', ...
                  sprintf('%g [K]',cellData.const.Tref), ...
                  'Simulation temperature');
  model.param.set('NEGATIVE_ELECTRODE_PARAMETERS', '0', '------------');
  model.param.set('L_neg', ...
                  sprintf('%g [m]',cellData.neg.L), ...
                  'Length of negative electrode');
  model.param.set('Rs_neg', ...
                  sprintf('%g [m]',cellData.neg.Rs), ...
                  'Particle radius Negative');
  model.param.set('Eact_Ds_neg', ...
                  sprintf('%g [J/mol]',cellData.neg.Eact_Ds), ...
                  'Activation energy');
  model.param.set('Ds_neg_ref', ...
                  sprintf('%g [m^2/s]',cellData.neg.Ds_ref), ...
                  'Solid phase Li-diffusivity Negative');
  model.param.set('Eact_sigma_neg', ...
                  sprintf('%g [J/mol]',cellData.neg.Eact_sigma), ...
                  'Activation energy');
  model.param.set('sigma_neg_ref', ...
                  sprintf('%g [S/m]',cellData.neg.sigma_ref), ...
                  'Solid phase conductivity Negative');
  model.param.set('Eact_k_neg', ...
                  sprintf('%g [J/mol]',cellData.neg.Eact_k), ...
                  'Activation energy');
  model.param.set('k_neg_norm_ref', ...
                  sprintf('%g [mol/m^2/s]',cellData.neg.k_norm_ref), ...
                  'Reaction rate coefficient Negative');
  model.param.set('as_neg', '3*(eps_s_neg)/Rs_neg', ...
                  'Specific surface area Negative');
  model.param.set('eps_e_neg', cellData.neg.eps_e, ...
                  'Electrolyte phase vol-fraction Negative');
  model.param.set('eps_s_neg', cellData.neg.eps_s, ...
                  'Solid phase vol-fraction Negative');
  model.param.set('csmax_neg', ...
                  sprintf('%g [mol/m^3]',cellData.neg.csmax), ...
                  'Max solid phase concentration Negative');
  model.param.set('theta_init_neg', theta_init_neg, ...
                  'Initial Negative State of Charge');
  model.param.set('cs0_neg', 'csmax_neg*theta_init_neg', ...
                  'Initial solid phase conc Negative');
  model.param.set('alpha_neg', cellData.neg.alpha, ...
                  'Charge transfer coefficient');
  model.param.set('brug_sigma_neg', cellData.neg.brug_sigma);
  model.param.set('brug_De_neg', cellData.neg.brug_De, ...
                  'Bruggeman coefficient');
  model.param.set('brug_kappa_neg', cellData.neg.brug_kappa);
  model.param.set('Rfilm_neg', ...
                  sprintf('%g [V*m^2/A]',cellData.neg.Rfilm), ...
                  'Film resistance');
  model.param.set('SEPARATOR_PARAMETERS', '0', '---------------------');
  model.param.set('L_sep', ...
                  sprintf('%g [m]',cellData.sep.L), ...
                  'Length of separator');
  model.param.set('eps_e_sep', cellData.sep.eps_e, 'Separator porosity');
  model.param.set('brug_De_sep', cellData.sep.brug_De);
  model.param.set('brug_kappa_sep', cellData.sep.brug_kappa);
  model.param.set('POSITIVE_ELECTRODE_PARAMETERS', '0', '------------');
  model.param.set('L_pos', ...
                  sprintf('%g [m]',cellData.pos.L), ...
                  'Length of positive electrode');
  model.param.set('Rs_pos', ...
                  sprintf('%g [m]',cellData.pos.Rs), ...
                  'Particle radius Positive');
  model.param.set('Eact_Ds_pos', ...
                  sprintf('%g [J/mol]',cellData.pos.Eact_Ds), ...
                  'Activation energy');
  model.param.set('Ds_pos_ref', ...
                  sprintf('%g [m^2/s]',cellData.pos.Ds_ref), ...
                  'Solid phase Li-diffusivity Positive');
  model.param.set('Eact_sigma_pos', ...
                  sprintf('%g [J/mol]',cellData.pos.Eact_sigma), ...
                  'Activation energy');
  model.param.set('sigma_pos_ref', ...
                  sprintf('%g [S/m]',cellData.pos.sigma_ref), ...
                  'Solid phase conductivity Positive');
  model.param.set('Eact_k_pos', ...
                  sprintf('%g [J/mol]',cellData.pos.Eact_k), ...
                  'Activation energy');
  model.param.set('k_pos_norm_ref', ...
                  sprintf('%g [mol/m^2/s]',cellData.pos.k_norm_ref), ...
                  'Reaction rate coefficient Positive');
  model.param.set('as_pos', '3*(eps_s_pos)/Rs_pos', ...
                  'Specific surface area Negative');
  model.param.set('eps_e_pos', cellData.pos.eps_e, ...
                  'Electrolyte phase vol-fraction Positive');
  model.param.set('eps_s_pos', cellData.pos.eps_s, ...
                  'Solid phase vol-fraction Positive');
  model.param.set('csmax_pos', ...
                  sprintf('%g [mol/m^3]',cellData.pos.csmax), ...
                  'Max solid phase concentration Positive');
  model.param.set('theta_init_pos', theta_init_pos, ...
                  'Initial Positive State of Charge');
  model.param.set('cs0_pos', 'csmax_pos*theta_init_pos', ...
                  'Initial solid phase conc Positive');
  model.param.set('alpha_pos', cellData.pos.alpha, ...
                  'Charge transfer coefficient Positive');
  model.param.set('brug_sigma_pos', cellData.pos.brug_sigma);
  model.param.set('brug_De_pos', cellData.pos.brug_De);
  model.param.set('brug_kappa_pos', cellData.pos.brug_kappa);
  model.param.set('Rfilm_pos', ...
                  sprintf('%g [V*m^2/A]',cellData.pos.Rfilm), ...
                  'Film resistance');
  model.param.set('ELECTROLYTE_PARAMETERS', '0', '-------------------');
  model.param.set('Eact_De', ...
                  sprintf('%g [J/mol]',cellData.const.Eact_De), ...
                  'Activation energy');
  model.param.set('De_ref', ...
                  sprintf('%g [m^2/s]',cellData.const.De_ref), ...
                  'Salt diffusivity in Electrolyte');
  model.param.set('Eact_kappa', ...
                  sprintf('%g [J/mol]',cellData.const.Eact_kappa), ...
                  'Activation energy');
  model.param.set('kappa_D_fact', '2*R*T/F*(1+dlnfdlnc)*(1-t_plus)', ...
                  'Junction potantial coefficient');
  model.param.set('dlnfdlnc', '0', ...
                  'Activity factor concentration variation');
  model.param.set('t_plus', cellData.const.t_plus, ...
                  'Cationic transport number');
  model.param.set('c_e0', ...
                  sprintf('%g [mol/m^3]',cellData.const.ce0), ...
                  'Initial electrolyte salt concentration');
  model.param.set('SIMULATION_CONTROL', '0', '-----------------------');
  model.param.set('xnorm', '1[m]');
  % model.param.set('i_charge', '-i_1C', 'Charge current');
  % model.param.set('i_disch', 'i_1C', 'Discharge current');
  % model.param.set('i_app_pulse', 'i_1C*input_current(t)');
  % model.param.set('i_app', 'i_app_pulse');
  model.param.set('Acell', ...
                  sprintf('%g [m^2]',cellData.const.Acell), ...
                  'Cross sectional electrode area');
  model.param.set('i_app', 'input_current(t)/(Acell)');  

  %% Comsol global functions: Filled in with data from cellData
  if ischar(cellData.neg.Uocp),
    [UocpNeg,UocpNegVar] = parseString(cellData.neg.Uocp);
    model.func.create('UocpNeg', 'Analytic');
    model.func('UocpNeg').set('funcname', 'UocpNeg');
    model.func('UocpNeg').set('expr', UocpNeg);
    model.func('UocpNeg').set('args', {UocpNegVar});
    model.func('UocpNeg').set('fununit', 'V');
    model.func('UocpNeg').set('plotargs', {'x' '0.015' '0.7247'});
  else
    error('Don''t yet know how to handle non-string Uocp');
  end

  if ischar(cellData.pos.Uocp),
    [UocpPos,UocpPosVar] = parseString(cellData.pos.Uocp);
    model.func.create('UocpPos', 'Analytic');
    model.func('UocpPos').set('funcname', 'UocpPos');
    model.func('UocpPos').set('expr', UocpPos);
    model.func('UocpPos').set('args', {UocpPosVar});
    model.func('UocpPos').set('fununit', 'V');
    model.func('UocpPos').set('plotargs', {'x' '0.443' '0.9515'});
  else
    error('Don''t yet know how to handle non-string Uocp');
  end

  if ischar(cellData.const.kappa_ref),
    [kappa,kappaVar] = parseString(cellData.const.kappa_ref);
    model.func.create('Kappa', 'Analytic');
    model.func('Kappa').set('funcname', 'Kappa_ref');
    model.func('Kappa').set('expr', kappa);
    model.func('Kappa').set('args', {kappaVar});
    model.func('Kappa').set('argunit', 'mol/m^3');
    model.func('Kappa').set('fununit', 'S/m');
    model.func('Kappa').set('plotargs', {'ce' '0' '3000'});
  else
    error('Don''t yet know how to handle non-string kappa_ref');
  end

  if length(tvec) == length(ivec),
    ivec = ivec(1:end-1);
  end
  v1 = tvec(1:end-1); v2 = tvec(2:end); v3 = ivec;
  Iapp = sprintf('''%g'' ''%g'' ''%g'';',[v1(:),v2(:),v3(:)]');
  Iapp = eval(sprintf('{ %s }',Iapp));
  model.func.create('Iapp', 'Piecewise');
  model.func('Iapp').set('funcname', 'input_current');
  model.func('Iapp').set('arg', 't');
  model.func('Iapp').set('extrap', 'periodic');
  model.func('Iapp').set('smooth', 'contd2');
  model.func('Iapp').set('smoothzone', '3E-7');
  model.func('Iapp').set('pieces', Iapp);
  model.func('Iapp').set('argunit', 's');

  model.func.create('NicePow', 'Analytic');
  model.func('NicePow').set('plotargs', {'x' '-1' '1'; 'y' '0' '1'});
  model.func('NicePow').set('funcname', 'NicePow');
  model.func('NicePow').set('args', {'x' 'y'});
  model.func('NicePow').set('expr', '(step(x)*abs(x))^(y)');

  model.func.create('Step', 'Step');
  model.func('Step').set('funcname', 'step');

  model.modelNode.create('mod1D');
  model.modelNode('mod1D').name('1D cellData model');
  model.modelNode('mod1D').identifier('main_dim');

  model.modelNode.create('mod2D');
  model.modelNode('mod2D').name('2D electrode model');
  model.modelNode('mod2D').identifier('pseudo_dim');

  %% Model geometries and meshes: 1D and 2D electrodes/ 1D separator
  %  Should never need to change this section
  model.geom.create('geom1D', 1);
  model.geom('geom1D').model('mod1D');
  model.geom('geom1D').feature.create('intervals1D', 'Interval');
  model.geom('geom1D').feature('intervals1D').set('intervals', 'many');
  model.geom('geom1D').feature('intervals1D').set('p', 'range(0,1,3)');
  model.geom('geom1D').run;

  model.mesh.create('mesh1D', 'geom1D');
  model.mesh('mesh1D').feature.create('E1', 'Edge');
  model.mesh('mesh1D').feature('E1').feature.create('D1', 'Distribution');
  model.mesh('mesh1D').feature('E1').feature('D1').selection.set(1);
  model.mesh('mesh1D').feature('E1').feature('D1').set('type','predefined');
  model.mesh('mesh1D').feature('E1').feature('D1').set('elemcount','20');
  model.mesh('mesh1D').feature('E1').feature('D1').set('elemratio','0.5');
  model.mesh('mesh1D').feature('E1').feature('D1').set('method','geometric');
  model.mesh('mesh1D').feature('E1').feature.create('D2', 'Distribution');
  model.mesh('mesh1D').feature('E1').feature('D2').selection.set(2);
  model.mesh('mesh1D').feature('E1').feature('D2').set('type','predefined');
  model.mesh('mesh1D').feature('E1').feature('D2').set('elemcount','20');
  model.mesh('mesh1D').feature('E1').feature('D2').set('elemratio','0.5');
  model.mesh('mesh1D').feature('E1').feature('D2').set('symmetric',true);
  model.mesh('mesh1D').feature('E1').feature('D2').set('reverse', true);
  model.mesh('mesh1D').feature('E1').feature.create('D3','Distribution');
  model.mesh('mesh1D').feature('E1').feature('D3').selection.set(3);
  model.mesh('mesh1D').feature('E1').feature('D3').set('type','predefined');
  model.mesh('mesh1D').feature('E1').feature('D3').set('elemcount','20');
  model.mesh('mesh1D').feature('E1').feature('D3').set('elemratio','0.5');
  model.mesh('mesh1D').feature('E1').feature('D3').set('reverse', true);
  model.mesh('mesh1D').run;

  model.geom.create('geom2D', 2);
  model.geom('geom2D').feature.create('squareNeg', 'Square');
  model.geom('geom2D').feature.create('squarePos', 'Square');
  model.geom('geom2D').feature('squareNeg').name('Negative electrode');
  model.geom('geom2D').feature('squarePos').name('Positive electrode');
  model.geom('geom2D').feature('squareNeg').set('pos', '0.0,0.0');
  model.geom('geom2D').feature('squarePos').set('pos', '1.5,0.0');
  model.geom('geom2D').run;

  model.mesh.create('mesh2D', 'geom2D');
  model.mesh('mesh2D').feature.create('ftri1', 'FreeTri');
  model.mesh('mesh2D').feature('ftri1').feature.create('D1', 'Distribution');
  model.mesh('mesh2D').feature('ftri1').feature('D1').selection.set([3 7]);
  model.mesh('mesh2D').feature('size').set('hmax', '0.0925');
  model.mesh('mesh2D').feature('size').set('hmin', '3.13E-4');
  model.mesh('mesh2D').feature('size').set('hcurve', '0.25');
  model.mesh('mesh2D').feature('size').set('hgrad', '1.25');
  model.mesh('mesh2D').feature('size').set('hauto', '3');
  model.mesh('mesh2D').feature('ftri1').set('yscale', '0.5');
  model.mesh('mesh2D').feature('ftri1').feature('D1').set('numelem', '20');
  model.mesh('mesh2D').run;

  %% 1D region variables and algebraic equations
  model.variable.create('Neg1DVars');
  model.variable('Neg1DVars').name('Negative electrode variables');
  model.variable('Neg1DVars').model('mod1D');
  model.variable('Neg1DVars').set('L', 'L_neg');
  model.variable('Neg1DVars').set('as', 'as_neg');
  model.variable('Neg1DVars').set('eps_e', 'eps_e_neg');
  model.variable('Neg1DVars').set('cse', ...
                 'pseudo_dim.cse_neg(extrcpl_source_cse)');
  model.variable('Neg1DVars').set('soc', 'cse/csmax_neg');
  model.variable('Neg1DVars').set('sigma_eff', ...
                 ['sigma_neg_ref*exp(Eact_sigma_neg/R*(1/Tref-1/T))*'...
                 'eps_s_neg^brug_sigma_neg']);
  model.variable('Neg1DVars').set('kappa_eff', ...
                 ['Kappa_ref(c_e)*exp(Eact_kappa/R*(1/Tref-1/T))*'...
                 'eps_e_neg^brug_kappa_neg']);
  model.variable('Neg1DVars').set('De_eff', ...
                 ['De_ref*exp(Eact_De/R*(1/Tref-1/T))*'...
                 'eps_e_neg^brug_De_neg']);
  model.variable('Neg1DVars').set('j0', ...
                 ['k_neg_norm_ref*exp(Eact_k_neg/R*(1/Tref-1/T))*'...
                 'NicePow((csmax_neg-cse)/csmax_neg,1-alpha_neg)*'...
                 'NicePow(cse/csmax_neg,alpha_neg)*'...
                 'NicePow(c_e/c_e0,1-alpha_neg)']);
  model.variable('Neg1DVars').set('j', ...
                 ['j0*(exp((1-alpha_neg)*F*eta/(R*T))-'...
                 'exp(-alpha_neg*F*eta/(R*T)))']);
  model.variable('Neg1DVars').set('extrcpl_source_influx_neg', '-j');
  model.variable('Neg1DVars').set('eta_0', ...
                 'eta-phi_s+phi_e+UocpNeg(soc)+j*Rfilm_neg*F');
  model.variable('Neg1DVars').selection.geom('geom1D', 1);
  model.variable('Neg1DVars').selection.set(1);

  model.variable.create('Sep1DVars');
  model.variable('Sep1DVars').name('Separator variables');
  model.variable('Sep1DVars').model('mod1D');
  model.variable('Sep1DVars').set('L', 'L_sep');
  model.variable('Sep1DVars').set('as', '0 [m^-1]');
  model.variable('Sep1DVars').set('eps_e', 'eps_e_sep');
  model.variable('Sep1DVars').set('kappa_eff', ...
                 ['Kappa_ref(c_e)*exp(Eact_kappa/R*(1/Tref-1/T))*'...
                 'eps_e_sep^brug_kappa_sep']);
  model.variable('Sep1DVars').set('De_eff', ...
                 ['De_ref*exp(Eact_De/R*(1/Tref-1/T))*'...
                 'eps_e_sep^brug_De_sep']);
  model.variable('Sep1DVars').set('j', '0 [mol/m^2/s]');
  model.variable('Sep1DVars').selection.geom('geom1D', 1);
  model.variable('Sep1DVars').selection.set(2);

  model.variable.create('Pos1DVars');
  model.variable('Pos1DVars').name('Positive electrode variables');
  model.variable('Pos1DVars').model('mod1D');
  model.variable('Pos1DVars').set('L', 'L_pos');
  model.variable('Pos1DVars').set('as', 'as_pos');
  model.variable('Pos1DVars').set('eps_e', 'eps_e_pos');
  model.variable('Pos1DVars').set('cse', ...
                 'pseudo_dim.cse_pos(extrcpl_source_cse)');
  model.variable('Pos1DVars').set('soc', 'cse/csmax_pos');
  model.variable('Pos1DVars').set('sigma_eff', ...
                 ['sigma_pos_ref*exp(Eact_sigma_pos/R*(1/Tref-1/T))*'...
                 'eps_s_pos^brug_sigma_pos']);
  model.variable('Pos1DVars').set('kappa_eff', ...
                 ['Kappa_ref(c_e)*exp(Eact_kappa/R*(1/Tref-1/T))*'...
                 'eps_e_pos^brug_kappa_pos']);
  model.variable('Pos1DVars').set('De_eff', ...
                 ['De_ref*exp(Eact_De/R*(1/Tref-1/T))*'...
                 'eps_e_pos^brug_De_pos']);
  model.variable('Pos1DVars').set('eta_old','phi_s-phi_e-UocpPos(soc)');
  model.variable('Pos1DVars').set('j0', ...
                 ['k_pos_norm_ref*exp(Eact_k_pos/R*(1/Tref-1/T))*'...
                 'NicePow((csmax_pos-cse)/csmax_pos,1-alpha_pos)*'...
                 'NicePow(cse/csmax_pos,alpha_pos)*'...
                 'NicePow(c_e/c_e0,1-alpha_pos)']);
  model.variable('Pos1DVars').set('j', ...
                 ['j0*(exp((1-alpha_pos)*F*eta/(R*T))-'...
                 'exp(-alpha_pos*F*eta/(R*T)))']);
  model.variable('Pos1DVars').set('extrcpl_source_influx_pos', '-j');
  model.variable('Pos1DVars').set('eta_0', ...
                 'eta-phi_s+phi_e+UocpPos(soc)+j*Rfilm_pos*F');
  model.variable('Pos1DVars').selection.geom('geom1D', 1);
  model.variable('Pos1DVars').selection.set(3);

  %% 2D region variables and algebraic equations
  model.variable.create('Neg2DVars');
  model.variable('Neg2DVars').name('Negative electrode variables');
  model.variable('Neg2DVars').model('mod2D');
  model.variable('Neg2DVars').set('Rs', 'Rs_neg');
  model.variable('Neg2DVars').set('Ds', ...
                 'Ds_neg_ref*exp(Eact_Ds_neg/R*(1/Tref-1/T))');
  model.variable('Neg2DVars').set('c_s0', 'cs0_neg');
  model.variable('Neg2DVars').set('cs_scale_neg', 'c_s*y^2');
  model.variable('Neg2DVars').set('ynorm', '1[m]');
  model.variable('Neg2DVars').selection.geom('geom2D', 2);
  model.variable('Neg2DVars').selection.set(1);

  model.variable.create('Pos2DVars');
  model.variable('Pos2DVars').name('Positive electrode variables');
  model.variable('Pos2DVars').model('mod2D');
  model.variable('Pos2DVars').set('Rs', 'Rs_pos');
  model.variable('Pos2DVars').set('Ds', ...
                 'Ds_pos_ref*exp(Eact_Ds_pos/R*(1/Tref-1/T))');
  model.variable('Pos2DVars').set('c_s0', 'cs0_pos');
  model.variable('Pos2DVars').set('ynorm', '1[m]');
  model.variable('Pos2DVars').selection.geom('geom2D', 2);
  model.variable('Pos2DVars').selection.set(2);

  model.variable.create('Neg2DBoundary');
  model.variable('Neg2DBoundary').name('Negative electrode boundary');
  model.variable('Neg2DBoundary').model('mod2D');
  model.variable('Neg2DBoundary').set('cse_influx', ...
                 'main_dim.cse_neg(extrcpl_source_influx_neg)');
  model.variable('Neg2DBoundary').set('extrcpl_source_cse', 'c_s');
  model.variable('Neg2DBoundary').selection.geom('geom2D', 1);
  model.variable('Neg2DBoundary').selection.set(3);

  model.variable.create('Pos2DBoundary');
  model.variable('Pos2DBoundary').name('Positive electrode boundary');
  model.variable('Pos2DBoundary').model('mod2D');
  model.variable('Pos2DBoundary').set('cse_influx', ...
                 'main_dim.cse_pos(extrcpl_source_influx_pos)');
  model.variable('Pos2DBoundary').set('extrcpl_source_cse', 'c_s');
  model.variable('Pos2DBoundary').selection.geom('geom2D', 1);
  model.variable('Pos2DBoundary').selection.set(7);

  %% Connection between 1D and 2D domains
  %  Should never need to modify this
  model.cpl.create('linExtrudeNeg1D', 'LinearExtrusion', 'geom1D');
  model.cpl('linExtrudeNeg1D').selection.set(1);
  model.cpl('linExtrudeNeg1D').name('Linear Extrusion Neg');
  model.cpl('linExtrudeNeg1D').set('dstgeom', 'geom2D');
  model.cpl('linExtrudeNeg1D').set('dstframe', 'material');
  model.cpl('linExtrudeNeg1D').set('srcframe', 'material');
  model.cpl('linExtrudeNeg1D').set('opname', 'cse_neg');
  model.cpl('linExtrudeNeg1D').selection('dstvertex1').set(2);
  model.cpl('linExtrudeNeg1D').selection('dstvertex2').set(4);
  model.cpl('linExtrudeNeg1D').selection('srcvertex1').geom('geom1D',0);
  model.cpl('linExtrudeNeg1D').selection('srcvertex1').set(1);
  model.cpl('linExtrudeNeg1D').selection('srcvertex2').geom('geom1D',0);
  model.cpl('linExtrudeNeg1D').selection('srcvertex2').set(2);

  model.cpl.create('linExtrudePos1D', 'LinearExtrusion', 'geom1D');
  model.cpl('linExtrudePos1D').selection.set(3);
  model.cpl('linExtrudePos1D').name('Linear Extrusion Pos');
  model.cpl('linExtrudePos1D').set('dstgeom', 'geom2D');
  model.cpl('linExtrudePos1D').set('dstframe', 'material');
  model.cpl('linExtrudePos1D').set('srcframe', 'material');
  model.cpl('linExtrudePos1D').set('opname', 'cse_pos');
  model.cpl('linExtrudePos1D').selection('dstvertex1').set(6);
  model.cpl('linExtrudePos1D').selection('dstvertex2').set(8);
  model.cpl('linExtrudePos1D').selection('srcvertex1').geom('geom1D',0);
  model.cpl('linExtrudePos1D').selection('srcvertex1').set(3);
  model.cpl('linExtrudePos1D').selection('srcvertex2').geom('geom1D',0);
  model.cpl('linExtrudePos1D').selection('srcvertex2').set(4);

  model.cpl.create('linExtrudeNeg2D', 'LinearExtrusion', 'geom2D');
  model.cpl('linExtrudeNeg2D').selection.geom('geom2D', 1);
  model.cpl('linExtrudeNeg2D').selection.set(3);
  model.cpl('linExtrudeNeg2D').name('Linear Extrusion Neg');
  model.cpl('linExtrudeNeg2D').set('dstgeom', 'geom1D');
  model.cpl('linExtrudeNeg2D').set('dstframe', 'material');
  model.cpl('linExtrudeNeg2D').set('srcframe', 'material');
  model.cpl('linExtrudeNeg2D').set('opname', 'cse_neg');
  model.cpl('linExtrudeNeg2D').selection('dstvertex1').set(1);
  model.cpl('linExtrudeNeg2D').selection('dstvertex2').set(2);
  model.cpl('linExtrudeNeg2D').selection('srcvertex1').geom('geom2D',0);
  model.cpl('linExtrudeNeg2D').selection('srcvertex1').set(2);
  model.cpl('linExtrudeNeg2D').selection('srcvertex2').geom('geom2D',0);
  model.cpl('linExtrudeNeg2D').selection('srcvertex2').set(4);
  model.cpl('linExtrudeNeg2D').selection('srcvertex3').geom('geom2D',0);

  model.cpl.create('linExtrudePos2D', 'LinearExtrusion', 'geom2D');
  model.cpl('linExtrudePos2D').selection.geom('geom2D', 1);
  model.cpl('linExtrudePos2D').selection.set(7);
  model.cpl('linExtrudePos2D').name('Linear Extrusion Pos');
  model.cpl('linExtrudePos2D').set('dstgeom', 'geom1D');
  model.cpl('linExtrudePos2D').set('dstframe', 'material');
  model.cpl('linExtrudePos2D').set('srcframe', 'material');
  model.cpl('linExtrudePos2D').set('opname', 'cse_pos');
  model.cpl('linExtrudePos2D').selection('dstvertex1').set(3);
  model.cpl('linExtrudePos2D').selection('dstvertex2').set(4);
  model.cpl('linExtrudePos2D').selection('srcvertex1').geom('geom2D',0);
  model.cpl('linExtrudePos2D').selection('srcvertex1').set(6);
  model.cpl('linExtrudePos2D').selection('srcvertex2').geom('geom2D',0);
  model.cpl('linExtrudePos2D').selection('srcvertex2').set(8);
  model.cpl('linExtrudePos2D').selection('srcvertex3').geom('geom2D',0);

  %% PDEs that COMSOL will solve
  model.physics.create('phi_s_eq', 'GeneralFormPDE', 'geom1D');
  model.physics('phi_s_eq').name('Solid charge conservation');
  model.physics('phi_s_eq').identifier('phi_s');
  model.physics('phi_s_eq').field('dimensionless').field('phi_s');
  model.physics('phi_s_eq').field('dimensionless').component({'phi_s'});
  model.physics('phi_s_eq').selection.set([1 3]);
  model.physics('phi_s_eq').feature.create('init2', 'init', 1);
  model.physics('phi_s_eq').feature('init2').selection.set(3);
  model.physics('phi_s_eq').feature.create('cons1', 'Constraint', 0);
  model.physics('phi_s_eq').feature('cons1').selection.set(1);
  model.physics('phi_s_eq').feature.create('flux1', 'FluxBoundary', 0);
  model.physics('phi_s_eq').feature('flux1').selection.set(4);
  model.physics('phi_s_eq').prop('ShapeProperty').set('boundaryFlux','0');
  model.physics('phi_s_eq').prop('Units').set('SourceTermQuantity',...
                                              'currentsource');
  model.physics('phi_s_eq').prop('Units').set(...
                      'DependentVariableQuantity', 'electricpotential');
  model.physics('phi_s_eq').prop('Units').set('CustomSourceTermUnit','V');
  model.physics('phi_s_eq').feature('gfeq1').set('f','-j*as*(L/xnorm)*F');
  model.physics('phi_s_eq').feature('gfeq1').set('da','0');
  model.physics('phi_s_eq').feature('gfeq1').set('Ga', ...
                                         '-sigma_eff/(L/xnorm)*phi_sx');
  model.physics('phi_s_eq').feature('gfeq1').name('General Form PDE');
  model.physics('phi_s_eq').feature('init2').set('phi_s', ...
                     'UocpPos(theta_init_pos)-UocpNeg(theta_init_neg)');
  model.physics('phi_s_eq').feature('cons1').set('R', 'phi_s');
  model.physics('phi_s_eq').feature('flux1').set('g', '-i_app');

  model.physics.create('phi_e_eq', 'GeneralFormPDE', 'geom1D');
  model.physics('phi_e_eq').name('Electrolyte charge conservation');
  model.physics('phi_e_eq').identifier('phi_e');
  model.physics('phi_e_eq').field('dimensionless').field('phi_e');
  model.physics('phi_e_eq').field('dimensionless').component({'phi_e'});
  model.physics('phi_e_eq').prop('ShapeProperty').set('boundaryFlux','0');
  model.physics('phi_e_eq').prop('Units').set('SourceTermQuantity', ...
                                              'currentsource');
  model.physics('phi_e_eq').prop('Units').set(...
                      'DependentVariableQuantity', 'electricpotential');
  model.physics('phi_e_eq').feature('gfeq1').set('f','j*as*(L/xnorm)*F');
  model.physics('phi_e_eq').feature('gfeq1').set('da', '0');
  model.physics('phi_e_eq').feature('gfeq1').set('Ga', ...
               '-kappa_eff/(L/xnorm)*(phi_ex-kappa_D_fact*1/c_e*c_ex)');
  model.physics('phi_e_eq').feature('gfeq1').name('General Form PDE');
  model.physics('phi_e_eq').feature('init1').set('phi_e', ...
                                            '-UocpNeg(theta_init_neg)');

  model.physics.create('c_s_eq', 'ConvectionDiffusionEquation','geom2D');
  model.physics('c_s_eq').name('Solid diffusion');
  model.physics('c_s_eq').identifier('c_s');
  model.physics('c_s_eq').field('dimensionless').field('c_s');
  model.physics('c_s_eq').feature.create('flux1', 'FluxBoundary', 1);
  model.physics('c_s_eq').feature('flux1').selection.set([3 7]);
  model.physics('c_s_eq').prop('ShapeProperty').set('boundaryFlux','0');
  model.physics('c_s_eq').prop('Units').set('SourceTermQuantity', ...
                                            'reactionrate');
  model.physics('c_s_eq').prop('Units').set(...
                          'DependentVariableQuantity', 'concentration');
  model.physics('c_s_eq').feature('cdeq1').set('c', ...
                                       {'0' '0' '0' 'y^2*Ds/Rs/ynorm'});
  model.physics('c_s_eq').feature('cdeq1').set('f', '0');
  model.physics('c_s_eq').feature('cdeq1').set('da','y^2*Rs/(ynorm^3)');
  model.physics('c_s_eq').feature('cdeq1').name('Solid diffusion');
  model.physics('c_s_eq').feature('init1').set('c_s', 'c_s0');
  model.physics('c_s_eq').feature('flux1').set('g', 'cse_influx');

  model.physics.create('c_e_eq', 'ConvectionDiffusionEquation','geom1D');
  model.physics('c_e_eq').name('Electrolyte diffusion');
  model.physics('c_e_eq').identifier('c_e');
  model.physics('c_e_eq').field('dimensionless').field('c_e');
  model.physics('c_e_eq').feature.create('zflx2', 'ZeroFluxBoundary',0);
  model.physics('c_e_eq').feature('zflx2').selection.set([1 4]);
  model.physics('c_e_eq').prop('ShapeProperty').set('boundaryFlux','0');
  model.physics('c_e_eq').prop('Units').set('SourceTermQuantity', ...
                                            'reactionrate');
  model.physics('c_e_eq').prop('Units').set(...
                          'DependentVariableQuantity', 'concentration');
  model.physics('c_e_eq').feature('cdeq1').set('c', 'De_eff/(L/xnorm)');
  model.physics('c_e_eq').feature('cdeq1').set('f', ...
                                           'j*(L/xnorm)*as*(1-t_plus)');
  model.physics('c_e_eq').feature('cdeq1').set('da', 'eps_e*(L/xnorm)');
  model.physics('c_e_eq').feature('cdeq1').name(...
                                       'Convection-Diffusion Equation');
  model.physics('c_e_eq').feature('init1').set('c_e', 'c_e0');

  model.physics.create('eta_0_eq', 'DomainODE', 'geom1D');
  model.physics('eta_0_eq').name('Intercalation overpotential');
  model.physics('eta_0_eq').identifier('eta_0');
  model.physics('eta_0_eq').field('dimensionless').field('eta');
  model.physics('eta_0_eq').field('dimensionless').component({'eta'});
  model.physics('eta_0_eq').selection.set([1 3]);
  model.physics('eta_0_eq').prop('Units').set('SourceTermQuantity', ...
                                              'electricpotential');
  model.physics('eta_0_eq').prop('Units').set(...
                      'DependentVariableQuantity', 'electricpotential');
  model.physics('eta_0_eq').feature('dode1').set('f', 'eta_0');
  model.physics('eta_0_eq').feature('dode1').set('da', '0');

  %% Set up simulation and run it
  model.study.create('study1');
  model.study('study1').feature.create('time', 'Transient');
  model.study('study1').feature('time').set('initstudyhide', 'on');
  model.study('study1').feature('time').set('initsolhide', 'on');
  model.study('study1').feature('time').set('notstudyhide', 'on');
  model.study('study1').feature('time').set('notsolhide', 'on');
  model.study('study1').name('Study');

  tvar = tvec + 0.01; 
  model.study('study1').feature('time').set('tlist', tvar);
  model.study('study1').feature('time').set('rtolactive', true);
  model.study('study1').feature('time').set('rtol', '0.001');

  model.sol.create('solution1');
  model.sol('solution1').study('study1');
  model.sol('solution1').attach('study1');
  model.sol('solution1').feature.create('st1', 'StudyStep');
  model.sol('solution1').feature.create('v1', 'Variables');
  model.sol('solution1').feature.create('t1', 'Time');
  model.sol('solution1').feature('t1').feature.create('fc1',...
                                                      'FullyCoupled');
  model.sol('solution1').feature('t1').feature.create('d1', 'Direct');
  model.sol('solution1').attach('study1');
  model.sol('solution1').feature('st1').name(...
                                   'Compile Equations: Time Dependent');
  model.sol('solution1').feature('st1').set('studystep', 'time');
  model.sol('solution1').feature('v1').feature('main_dim_phi_e').name(...
                                                      'main_dim.phi_e');
  model.sol('solution1').feature('v1').feature('main_dim_phi_e').set(...
                                              'variables', 'mod1_phi1');
  model.sol('solution1').feature('v1').feature('main_dim_phi_s').name(...
                                                      'main_dim.phi_s');
  model.sol('solution1').feature('v1').feature('main_dim_phi_s').set(...
                                              'variables', 'mod1_phi2');
  model.sol('solution1').feature('v1').feature('pseudo_dim_c_s').name(...
                                                      'pseudo_dim.c_s');
  model.sol('solution1').feature('v1').feature('pseudo_dim_c_s').set(...
                                                 'variables', 'mod2_u');
  model.sol('solution1').feature('v1').feature('main_dim_c_e').name(...
                                                        'main_dim.c_e');
  model.sol('solution1').feature('v1').feature('main_dim_c_e').set(...
                                                 'variables', 'mod1_u');
  model.sol('solution1').feature('v1').feature('main_dim_eta').name(...
                                                        'main_dim.eta');
  model.sol('solution1').feature('v1').feature('main_dim_eta').set(...
                                             'variables', 'main_dim_u');
  model.sol('solution1').feature('t1').set('atoludotactive',...
                      {'main_dim_phi_e' 'off' 'main_dim_phi_s' 'off' ...
     'pseudo_dim_c_s' 'off' 'main_dim_c_e' 'off' 'main_dim_eta' 'off'});
  model.sol('solution1').feature('t1').set('tlist', tvar);
  model.sol('solution1').feature('t1').set('atolglobalmethod','unscaled');
  model.sol('solution1').feature('t1').set('solfile', false);
  model.sol('solution1').feature('t1').set('atoludot', ...
                    {'main_dim_phi_e' '1e-3' 'main_dim_phi_s' '1e-3' ...
                     'pseudo_dim_c_s' '1e-3' 'main_dim_c_e' '1e-3' ...
                     'main_dim_eta' '1e-3'});
  model.sol('solution1').feature('t1').set('atolglobal', '0.001');
  model.sol('solution1').feature('t1').set('maxstepbdfactive', true);
  model.sol('solution1').feature('t1').set('storeudot', false);
  model.sol('solution1').feature('t1').set('rtol', '0.0010');
  model.sol('solution1').feature('t1').set('bwinitstepfrac', '1.0');
  model.sol('solution1').feature('t1').set('atolmethod', ...
                {'main_dim_phi_e' 'global' 'main_dim_phi_s' 'global' ...
                 'pseudo_dim_c_s' 'global' 'main_dim_c_e' 'global' ...
                 'main_dim_eta' 'global'});
  model.sol('solution1').feature('t1').set('control', 'time');
  model.sol('solution1').feature('t1').set('tstepsbdf', 'strict');
  model.sol('solution1').feature('t1').set('atol', ...
                     {'main_dim_phi_e' '1e-3' 'main_dim_phi_s' '1e-3'...
                      'pseudo_dim_c_s' '1e-3' 'main_dim_c_e' '1e-3' ...
                      'main_dim_eta' '1e-3'});
  model.sol('solution1').feature('t1').set('fieldselection', ...
                                           'main_dim_phi_e');
  model.sol('solution1').feature('t1').set('ewtrescale', false);
  model.sol('solution1').feature('t1').feature('fcDef').set(...
                                                  'probesel', 'manual');
  model.sol('solution1').feature('t1').feature('fc1').active(false);
  model.sol('solution1').feature('t1').feature('fc1').set(...
                                               'ratelimitactive', true);
  model.sol('solution1').feature('t1').feature('fc1').set('damp', '1.0');
  model.sol('solution1').feature('t1').feature('fc1').set('probesel',...
                                                              'manual');
  model.sol('solution1').feature('t1').feature('d1').set('errorchk','off');
  model.sol('solution1').feature('t1').feature('d1').set('linsolver','pardiso');
  model.sol('solution1').runAll;

  %% Set up results to be plotted
  model.result.dataset('dset1').name('1D solution pulses');
  model.result.dataset('dset2').name('2D solution pulses');

  model.result.create('volt_plot', 'PlotGroup1D');
  model.result('volt_plot').feature.create('ptgr2', 'PointGraph');
  model.result('volt_plot').feature('ptgr2').selection.set(4);
  model.result('volt_plot').name('Cell voltage');
  model.result('volt_plot').set('ylabel', 'Cell voltage (V)');
  model.result('volt_plot').set('xlabel', 'Time (s)');
  model.result('volt_plot').set('xlabelactive', true);
  model.result('volt_plot').set('title', 'Cell voltage profile');
  model.result('volt_plot').set('titletype', 'manual');
  model.result('volt_plot').set('ylabelactive', true);
  model.result('volt_plot').feature('ptgr2').name('Voltage vs t');
  model.result('volt_plot').feature('ptgr2').set('data', 'dset1');
  model.result('volt_plot').feature('ptgr2').set('titletype', 'custom');
  model.result('volt_plot').feature('ptgr2').set('looplevelinput', ...
                                                 {'all'});
  model.result('volt_plot').feature('ptgr2').set('descr', ...
                                                 'Cell voltage');
  model.result('volt_plot').feature('ptgr2').set('descractive', true);

  model.result.create('phi_s_plot', 'PlotGroup1D');
  model.result('phi_s_plot').feature.create('gr', 'LineGraph');
  model.result('phi_s_plot').feature('gr').selection.all;
  model.result('phi_s_plot').feature('gr').selection.all;
  model.result('phi_s_plot').name('Phi_s vs x');
  model.result('phi_s_plot').set('ylabel', 'Potential (V)');
  model.result('phi_s_plot').set('xlabel', ...
                                 'Normalized x-coordinate (unitless)');
  model.result('phi_s_plot').set('xlabelactive', true);
  model.result('phi_s_plot').set('title', ...
                                 'Solid potential \phi<sub>s</sub>');
  model.result('phi_s_plot').set('titletype', 'manual');
  model.result('phi_s_plot').set('ylabelactive', true);
  model.result('phi_s_plot').feature('gr').set('data', 'dset1');
  model.result('phi_s_plot').feature('gr').set('xdataexpr', 'x');
  model.result('phi_s_plot').feature('gr').set('xdata', 'expr');
  model.result('phi_s_plot').feature('gr').set('xdataunit', 'm');
  model.result('phi_s_plot').feature('gr').set('looplevelinput', ...
                                               {'all'});
  model.result('phi_s_plot').feature('gr').set('smooth', 'none');
  model.result('phi_s_plot').feature('gr').set('xdatadescr', ...
                                               'x-coordinate');
  model.result('phi_s_plot').feature('gr').selection.all;

  model.result.create('phi_e_plot', 'PlotGroup1D');
  model.result('phi_e_plot').feature.create('gr', 'LineGraph');
  model.result('phi_e_plot').feature('gr').selection.all;
  model.result('phi_e_plot').feature('gr').selection.all;
  model.result('phi_e_plot').name('Phi_e vs x');
  model.result('phi_e_plot').set('ylabel', 'Potential (V)');
  model.result('phi_e_plot').set('xlabel', ...
                                 'Normalized x-coordinate (unitless)');
  model.result('phi_e_plot').set('looplevelinput', {'all'});
  model.result('phi_e_plot').set('xlabelactive', true);
  model.result('phi_e_plot').set('title', ...
                              'Electrolyte potential \phi<sub>e</sub>');
  model.result('phi_e_plot').set('titletype', 'manual');
  model.result('phi_e_plot').set('ylabelactive', true);
  model.result('phi_e_plot').feature('gr').set('xdataexpr', 'x');
  model.result('phi_e_plot').feature('gr').set('descr', ...
                                            'Dependent variable phi_e');
  model.result('phi_e_plot').feature('gr').set('xdata', 'expr');
  model.result('phi_e_plot').feature('gr').set('expr', 'phi_e');
  model.result('phi_e_plot').feature('gr').set('xdataunit', 'm');
  model.result('phi_e_plot').feature('gr').set('smooth', 'none');
  model.result('phi_e_plot').feature('gr').set('xdatadescr', ...
                                               'x-coordinate');
  model.result('phi_e_plot').feature('gr').selection.all;

  model.result.create('c_e_plot', 'PlotGroup1D');
  model.result('c_e_plot').feature.create('gr', 'LineGraph');
  model.result('c_e_plot').feature('gr').selection.all;
  model.result('c_e_plot').feature('gr').selection.all;
  model.result('c_e_plot').name('c_e vs x');
  model.result('c_e_plot').set('ylabel', ...
                               'Concentration (mol m<sup>-3</sup>)');
  model.result('c_e_plot').set('xlabel', ...
                               'Normalized x-coordinate (unitless)');
  model.result('c_e_plot').set('looplevelinput', {'all'});
  model.result('c_e_plot').set('xlabelactive', true);
  model.result('c_e_plot').set('title', ...
                             'Electrolyte concentration c<sub>e</sub>');
  model.result('c_e_plot').set('titletype', 'manual');
  model.result('c_e_plot').set('ylabelactive', true);
  model.result('c_e_plot').set('looplevel', {'1'});
  model.result('c_e_plot').feature('gr').set('xdataexpr', 'x');
  model.result('c_e_plot').feature('gr').set('descr', ...
                                             'Dependent variable c_e');
  model.result('c_e_plot').feature('gr').set('unit', 'mol/m^3');
  model.result('c_e_plot').feature('gr').set('xdata', 'expr');
  model.result('c_e_plot').feature('gr').set('expr', 'c_e');
  model.result('c_e_plot').feature('gr').set('xdataunit', 'm');
  model.result('c_e_plot').feature('gr').set('smooth', 'none');
  model.result('c_e_plot').feature('gr').set('xdatadescr', ...
                                             'x-coordinate');
  model.result('c_e_plot').feature('gr').selection.all;

  model.result.create('c_se_plot', 'PlotGroup1D');
  model.result('c_se_plot').feature.create('gr', 'LineGraph');
  model.result('c_se_plot').feature('gr').selection.all;
  model.result('c_se_plot').feature('gr').selection.all;
  model.result('c_se_plot').name('c_{s,e} vs x');
  model.result('c_se_plot').set('ylabel', ...
                                'Concentration (mol m<sup>-3</sup>)');
  model.result('c_se_plot').set('xlabel', ...
                                'Normalized x-coordinate (unitless)');
  model.result('c_se_plot').set('xlabelactive', true);
  model.result('c_se_plot').set('title', ...
                         'Solid surface concentration c<sub>s,e</sub>');
  model.result('c_se_plot').set('titletype', 'manual');
  model.result('c_se_plot').set('ylabelactive', true);
  model.result('c_se_plot').feature('gr').set('xdataexpr','x/xnorm');
  model.result('c_se_plot').feature('gr').set('descr', '');
  model.result('c_se_plot').feature('gr').set('unit', 'mol/m^3');
  model.result('c_se_plot').feature('gr').set('title', ...
                               'Solid surface concentration (mol/m^3)');
  model.result('c_se_plot').feature('gr').set('titletype', 'manual');
  model.result('c_se_plot').feature('gr').set('xdata', 'expr');
  model.result('c_se_plot').feature('gr').set('expr', 'cse');
  model.result('c_se_plot').feature('gr').set('xdataunit', '1');
  model.result('c_se_plot').feature('gr').set('xdatadescractive', true);
  model.result('c_se_plot').feature('gr').set('smooth', 'none');
  model.result('c_se_plot').feature('gr').set('xdatadescr', ...
                                             'x-coordinate (unitless)');
  model.result('c_se_plot').feature('gr').selection.all;

  model.result.create('j_plot', 'PlotGroup1D');
  model.result('j_plot').feature.create('gr', 'LineGraph');
  model.result('j_plot').feature('gr').selection.all;
  model.result('j_plot').feature('gr').selection.all;
  model.result('j_plot').name('j vs x');
  model.result('j_plot').set('ylabel', ...
                            'Flux (mol m<sup>-2</sup> s<sup>-1</sup>)');
  model.result('j_plot').set('xlabel', ...
                             'Normalized x-coordinate (unitless)');
  model.result('j_plot').set('xlabelactive', true);
  model.result('j_plot').set('title', 'Intercalation flux j');
  model.result('j_plot').set('titletype', 'manual');
  model.result('j_plot').set('ylabelactive', true);
  model.result('j_plot').feature('gr').set('xdataexpr', 'x');
  model.result('j_plot').feature('gr').set('descr', '');
  model.result('j_plot').feature('gr').set('unit', 'mol/m^2/s');
  model.result('j_plot').feature('gr').set('xdata', 'expr');
  model.result('j_plot').feature('gr').set('expr', 'j');
  model.result('j_plot').feature('gr').set('xdataunit', 'm');
  model.result('j_plot').feature('gr').set('smooth', 'none');
  model.result('j_plot').feature('gr').set('xdatadescr', ...
                                              'x-coordinate');
  model.result('j_plot').feature('gr').selection.all;

  model.result.create('c_s_plot', 'PlotGroup2D');
  model.result('c_s_plot').feature.create('surf1', 'Surface');
  model.result('c_s_plot').name('c_s vs x');
  model.result('c_s_plot').set('xlabelactive', true);
  model.result('c_s_plot').set('ylabelactive', true);
  model.result('c_s_plot').set('titletype', 'custom');
  model.result('c_s_plot').set('typeintitle', false);
  model.result('c_s_plot').set('descriptionintitle', false);
  model.result('c_s_plot').set('ylabel', ...
                               'Normalized y-coordinate (unitless)');
  model.result('c_s_plot').set('looplevel', {'1'});
  model.result('c_s_plot').set('xlabel', ...
                               'Normalized x-coordinate (unitless)');
  model.result('c_s_plot').set('suffixintitle', ...
                               ': Solid concentration c<sub>s</sub>');
  model.result('c_s_plot').set('legendpos', 'rightdouble');
  model.result('c_s_plot').feature('surf1').set('colortable', ...
                                                'DiscoLight');

  %% Set up results to be exported
  model.result.export.create('plot1', 'Plot');
  model.result.export('plot1').name('Cell voltage versus t');
  model.result.export('plot1').set('filename', 'voltage.txt');

  model.result.export.create('data1', 'Data');
  model.result.export('data1').name('phi_s vs x');
  model.result.export('data1').set('expr', {'phi_s'});
  model.result.export('data1').set('solnum', {'1'});
  model.result.export('data1').set('filename', 'phi_s.txt');
  model.result.export('data1').set('descr', {'phi_s'});
  model.result.export('data1').set('unit', {''});

  model.result.export.create('data2', 'Data');
  model.result.export('data2').name('phi_e vs x');
  model.result.export('data2').set('expr', {'phi_e'});
  model.result.export('data2').set('solnum', {'1'});
  model.result.export('data2').set('filename', 'phi_e.txt');
  model.result.export('data2').set('descr', {'phi_e'});
  model.result.export('data2').set('unit', {''});

  model.result.export.create('data3', 'Data');
  model.result.export('data3').name('c_e vs x');
  model.result.export('data3').set('expr', {'c_e'});
  model.result.export('data3').set('solnum', {'1'});
  model.result.export('data3').set('filename', 'ce.txt');
  model.result.export('data3').set('descr', {'c_e'});
  model.result.export('data3').set('unit', {''});

  model.result.export.create('data4', 'Data');
  model.result.export('data4').name('c_{s,e} vs x');
  model.result.export('data4').set('expr', {'cse'});
  model.result.export('data4').set('solnum', {'1'});
  model.result.export('data4').set('filename', 'cse.txt');
  model.result.export('data4').set('descr', {'c_{s,e}'});
  model.result.export('data4').set('unit', {''});

  model.result.export.create('data5', 'Data');
  model.result.export('data5').name('j vs x');
  model.result.export('data5').set('expr', {'j'});
  model.result.export('data5').set('solnum', {'1'});
  model.result.export('data5').set('filename', 'j.txt');
  model.result.export('data5').set('descr', {'j'});
  model.result.export('data5').set('unit', {''});

  % model.sol('solution1').runAll;
  model.result('volt_plot').run;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % End of Code from COMSOL Gui Model
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Save Data from COMSOL Simulation
  input = mpheval(model,'i_app');
  data_neg = mpheval(model,{'soc','cse','j','c_e','phi_s','phi_e',...
                     'eta'},'Edim',1,'Selection',1);
  neg_locs = data_neg.p;
  data_sep = mpheval(model,{'c_e','phi_e'},'Edim',1,'Selection',2);
  sep_locs = data_sep.p;
  data_pos = mpheval(model,{'soc','cse','j','c_e','phi_s','phi_e',...
                     'eta'},'Edim',1,'Selection',3);
  pos_locs = data_pos.p;

  FOM.cellData = cellData;
  FOM.time = tvar-0.01;
  FOM.Iapp = input.d1(:,end);
  FOM.Vcell = data_pos.d5(:,end); Vcell = FOM.Vcell;

  % Calculate cell capacity, and with it calculate cell SOC
  Qneg = cellData.const.Acell*cellData.neg.L*cellData.neg.eps_s*...
         cellData.neg.csmax*(cellData.neg.theta100-cellData.neg.theta0)*...
         96485.3365/3600;
  Qpos = cellData.const.Acell*cellData.pos.L*cellData.pos.eps_s*...
         cellData.pos.csmax*(cellData.pos.theta0-cellData.pos.theta100)*...
         96485.3365/3600;
  Q = min(Qneg,Qpos);
  FOM.SOC = initSOC - [0; cumsum(ivec(1:end-1))]/Q/3600; 
  FOM.SOC = FOM.SOC(:);

  % Negative Electrode Saving
  FOM.j_neg = data_neg.d3;
  FOM.eta_neg = data_neg.d7;
  FOM.cse_neg = data_neg.d2;
  FOM.phis_neg = data_neg.d5;
  FOM.phise_neg = data_neg.d5 - data_neg.d6;

  % Positive Electrode Saving
  FOM.j_pos = data_pos.d3;
  FOM.eta_pos = data_pos.d7;
  FOM.cse_pos = data_pos.d2;
  FOM.phis_pos = data_pos.d5;
  FOM.phise_pos = data_pos.d5 - data_pos.d6;

  % Electrolyte saving
  FOM.phie = [data_neg.d6, data_sep.d2, data_pos.d6];
  FOM.ce = [data_neg.d4, data_sep.d1, data_pos.d4];

  % Saving locations of sampled variables in z/x coordinate systems
  locs.j_neg_locs = neg_locs(:);
  locs.eta_neg_locs = neg_locs(:);
  locs.cse_neg_locs = neg_locs(:);
  locs.phis_neg_locs = neg_locs(:);
  locs.phise_neg_locs = neg_locs(:);
  locs.j_pos_locs = 3 - pos_locs(:);
  locs.eta_pos_locs = 3 - pos_locs(:);
  locs.cse_pos_locs = 3 - pos_locs(:);
  locs.phis_pos_locs = 3 - pos_locs(:);
  locs.phise_pos_locs = 3 - pos_locs(:);

  locs.phie_locs = [cellData.neg.L*neg_locs(:); ...
                    cellData.neg.L + cellData.sep.L*(sep_locs(:)-1); ...
                    cellData.neg.L + cellData.sep.L + ...
                    cellData.pos.L*(pos_locs(:)-2)];
  locs.phie_neg_locs = cellData.neg.L*neg_locs(:);
  locs.phie_sep_locs = cellData.neg.L + cellData.sep.L*(sep_locs(:)-1); 
  locs.phie_pos_locs = cellData.neg.L + cellData.sep.L + ...
                       cellData.pos.L*(pos_locs(:)-2);
  locs.ce_locs = locs.phie_locs;
  locs.ce_neg_locs = locs.phie_neg_locs;
  locs.ce_sep_locs = locs.phie_sep_locs;
  locs.ce_pos_locs = locs.phie_pos_locs;
  FOM.locs = locs;

  function [str,strVar] = parseString(inStr)
    if inStr(1) == '{', % function plus derivative... keep function
      ind = find(inStr == ',',1,'first');
      inStr = inStr(2:ind-1);
    end
    ind1 = find(inStr == '@',1,'first');
    ind2 = find(inStr == ')',1,'first');
    varStr = inStr(ind1+1:ind2);
    varStr(varStr == '(')=[];
    varStr(varStr == ')')=[];
    inStr(ind1:ind2) = [];
    
    dotText = {'.*','./','.^','.+','.-'};
    for kk = 1:length(dotText)
      while(true)
        ind = strfind(inStr,dotText{kk});
        if isempty(ind), break; end;
        inStr(ind(1)) = [];
      end
    end
    str = inStr;
    strVar = varStr;   
  end
end
