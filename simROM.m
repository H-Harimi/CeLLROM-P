% function [vCell,ROMout] = simROM(ROMDIR,cellType,ik,Tk,SOC0,blend,dTerm)
%
%   Simulates the physics-based ROM in MATLAB and generates an output
%   structure. Given SOC0, i0...iN, T0...TN, predicts x0...xN, y0...yN, and
%   nonlinear outputs at time steps 0...N.
%
%   Inputs:
%     ROMDIR: path to ROM files produced by DRA
%     cellType: cell name/identifier (e.g., 'doyle')
%     SOC0: initial SOC (0..1, must presently correspond to a ROM file)
%     ik: N+1 vector of current (A)
%     Tk: N+1 vector of temperature (K)
%     blend: Use model blending if 1, or a single model only if 0
%     dTerm: 'uniform' or 'linear' -- how to compute "D" term in model
%
%   Outputs:
%     vCell: Cell voltage versus time
%     ROMout: All other internal cell variables

% Copyright (c) 2015 by Gregory L. Plett of the University of Colorado 
% Colorado Springs (UCCS). This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0.
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L., "Battery
% Management Systems, Volume I, Battery Modeling," Artech House, 2015.

function [vCell,ROMout] = simROM(ROMDIR,cellType,ik,Tk,SOC0,blend,dTerm)
  duration = size(ik,1);

  % Verify that input stimuli are equal length
  if size(ik,1) ~= size(Tk,1),
    error('Current and temperature profiles must be same length');
  end

  % Load base ROM parameter file from ROMDIR
  if blend, % If blend == 1, load a representative ROM for general params
    fileMask = sprintf('%scell*%s*.mat',ROMDIR,cellType);
    params = dir(fileMask);
    if isempty(params),
      error( 'No relevant ROMs exist in specified ROMDIR.' );
    else
      load( sprintf('%s%s',ROMDIR,params(1).name) );
      fprintf('Initial ROM file: %s loaded\n',params(1).name);
    end
  else % else load the specific needed ROM for all the cell params
    fileMask = sprintf('%scell*%s_ROM_%dSOC*.mat',ROMDIR,cellType,round(100*SOC0));
    params = dir(fileMask);
    if isempty(params),
      error( 'ROM file matching pattern "%s" does not exist.',fileMask );
    else
      load( sprintf('%s%s',ROMDIR,params(1).name) );
      fprintf('Initial ROM file: %s loaded\n',params(1).name);
    end
  end

  %% Load params into simple local variables
  % Convert Uocp_str to Uocp anonymous functions
  % Convert kappa_func_str to kappa_ref anonymous function
  cellData = ROM.cellData;
  cellData.neg.Uocp = eval(cellData.neg.Uocp_str);
  cellData.pos.Uocp = eval(cellData.pos.Uocp_str);
  if ~ischar(cellData.const.kappa_ref),
    cellData.const.kappa_ref = sprintf('(@(x)(%g))',...
      cellData.const.kappa_ref);
  end
  kappa_ref = eval(cellData.const.kappa_ref);

  % Define universal constants and copy param values from cellData
  % structure into variables for easier access
  F = 96485.3365; % Faraday's constant
  R = 8.3144621;  % Universal gas constant
  Tref = cellData.const.Tref; % Temperature at which this model was created
  t_plus = cellData.const.t_plus;
  Acell = cellData.const.Acell;
  Lneg =  cellData.neg.L; Lsep = cellData.sep.L; Lpos = cellData.pos.L;
  Ltot =  Lneg + Lsep + Lpos;
  csmaxNeg = cellData.neg.csmax; csmaxPos = cellData.pos.csmax;
  alphaNeg = cellData.neg.alpha; alphaPos = cellData.pos.alpha;
  ce0 = cellData.const.ce0;
  kNormNeg = cellData.neg.k_norm; kNeg = kNormNeg/csmaxNeg/ce0^(1-alphaNeg);
  kNormPos = cellData.pos.k_norm; kPos = kNormPos/csmaxPos/ce0^(1-alphaPos);
  asNeg = 3*cellData.neg.eps_s/cellData.neg.Rs;
  asPos = 3*cellData.pos.eps_s/cellData.pos.Rs;
  RfilmNeg = cellData.neg.Rfilm; RfilmPos = cellData.pos.Rfilm;
  EactKappa = cellData.const.Eact_kappa;
  EactSigma = cellData.neg.Eact_sigma;
  brugKappaNeg = cellData.neg.brug_kappa; brugKappaSep = cellData.sep.brug_kappa;
  brugKappaPos = cellData.pos.brug_kappa;
  brugSigmaNeg = cellData.neg.brug_sigma; brugSigmaPos = cellData.pos.brug_sigma;
    
  %% The DRA can produce models with transfer functions arranged in 
  % arbitrary order. We need to be able to find specific entries in the
  % model output vector in order to do nonlinear corrections and to compute
  % cell voltage, so here we define index variables that are used to
  % locate these particular model outputs.
  
  % First, where to find "j" in the negative and positive electrodes, and
  % at the corresponding current collectors x=0 and x=L.
  jNegInd      = find(strcmp(cellData.tf.name,'j_neg') == 1);  
  j0Ind        = find(strcmp(cellData.tf.name,'j_neg') == 1 & ...
                     cellData.tf.val == 0);
  jNegLocs     = cellData.tf.val(jNegInd); 
  jPosInd      = find(strcmp(cellData.tf.name,'j_pos') == 1);
  jLInd        = find(strcmp(cellData.tf.name,'j_pos') == 1 & ...
                     cellData.tf.val == 0);
  jPosLocs     = cellData.tf.val(jPosInd); 
  % Next, where to find "c_{s,e}"
  cseNegInd    = find(strcmp(cellData.tf.name,'cse_neg') == 1);
  cse0Ind      = find(strcmp(cellData.tf.name,'cse_neg') == 1 & ...
                     cellData.tf.val == 0);
  cseNegLocs   = cellData.tf.val(cseNegInd); 
  csePosInd    = find(strcmp(cellData.tf.name,'cse_pos') == 1);
  cseLInd      = find(strcmp(cellData.tf.name,'cse_pos') == 1 & ...
                     cellData.tf.val == 0);
  csePosLocs   = cellData.tf.val(csePosInd);
  % Where to find "\phi_s" and the gradient thereof
  phisNegInd   = find(strcmp(cellData.tf.name,'phis_neg') == 1);
  phisPosInd   = find(strcmp(cellData.tf.name,'phis_pos') == 1);  
  gradphisNegInd = find(strcmp(cellData.tf.name,'gradphis_neg') == 1);
  gradphisPosInd = find(strcmp(cellData.tf.name,'gradphis_pos') == 1);
  % Where to find "[\tilde\phi_e]_1" and the gradient thereof
  phie1Ind     = find(strcmp(cellData.tf.name,'phie1') == 1); 
  gradphie1Ind = find(strcmp(cellData.tf.name,'gradphie1') == 1); 
  % Where to find "\tilde c_e" and the gradient thereof
  ceInd        = find(strcmp(cellData.tf.name,'ce') == 1);
  gradceInd    = find(strcmp(cellData.tf.name,'gradce') == 1);
  ceLocs       = cellData.tf.val(ceInd); 
  % Where to find \phi_{s,e} at the negative-electrode current collector
  phiseNeg0Ind = find(strcmp(cellData.tf.name,'phise_neg') == 1 ...
                         & cellData.tf.val == 0);

  %% We store the electrochemical variables in the output structure 
  % returned to the calling program. Here, we initialize them to zero to
  % reserve memory for the results.
  jNeg = zeros(duration,length(jNegInd));
  j0 = zeros(duration,length(j0Ind));
  etaNeg = zeros(duration,length(jNegInd));
  eta0 = zeros(duration,length(j0Ind));
  cseNeg = zeros(duration,length(cseNegInd));
  cse0 = zeros(duration,length(cse0Ind));
  phisNeg = zeros(duration,length(phisNegInd));
  gradphisNeg = zeros(duration,length(gradphisNegInd));
  jPos = zeros(duration,length(jPosInd));
  jL = zeros(duration,length(jLInd));
  etaPos = zeros(duration,length(jPosInd));
  etaL = zeros(duration,length(jLInd));
  csePos = zeros(duration,length(csePosInd));
  cseL = zeros(duration,length(cseLInd));
  phisPos = zeros(duration,length(phisPosInd));
  gradphisPos = zeros(duration,length(gradphisPosInd));
  phie = zeros(duration,length(phie1Ind)+1);
  phietilde1 = zeros(duration,length(phie1Ind));
  gradphietilde1 = zeros(duration,length(gradphie1Ind)+2);
  ce = zeros(duration,length(ceInd));
  gradce = zeros(duration,length(gradceInd)+2);
  phiseNeg0 = zeros(duration,length(phiseNeg0Ind));
  vCell = zeros(duration,1);
    
  %% Verify that signals required to compute cell voltage exist
  % Need:
  %  1. cse_neg at the current collector (cse_0 not empty)
  %  2. cse_pos at the current collector (cse_L not empty)
  %  3. j at both current collectors (j_0 and j_L not empty)
  %  4. phi_e_tilde_1 at the L location
  %  5. ce at location 0 and L
  %  6. phise at the neg current collector (phise_neg_0 not empty)
  if isempty(cse0Ind), warning('cse(0) not found!'); end
  if isempty(cseLInd), warning('cse(L) not found!'); end
  if isempty(j0Ind),   warning('j(0) not found!');  end
  if isempty(jLInd),   warning('j(L) not found!'); end
  xphie = cellData.tf.val(phie1Ind);
  if xphie(1) == 0,
    warning('First phi_e xlocation should not be zero. Ignoring');
    phie1Ind = phie1Ind(2:end);
  end 
  if (xphie(end) > Ltot+eps) || (xphie(end) < Ltot-eps), 
    warning(['Assuming that last phi_e xlocation is at pos'...
             ' current collector but is at: %g'],xphie(end));
  end 
  xce = cellData.tf.val(ceInd);
  if xce(1) > 0,
    warning(['Assuming that first c_e xlocation is at neg'...
             ' current collector but is at: %g'],xce(1)); 
  end
  if (xce(end) > Ltot+eps) || (xce(end) < Ltot-eps),
    warning(['Assuming that last c_e xlocation is at pos'...
             ' current collector but is at: %g'],xce(end));
  end    
  if isempty(phiseNeg0Ind), warning('phise_neg(0) not found!'); end

  %% Verify that all signals needed to compute eta exist
  [C,~,IB] = intersect(jNegLocs,cseNegLocs);
  if length(C) < length(jNegLocs),
    warning('Cannot compute eta in negative electrode');
  else
    etaNegCseInd = IB;
    ceNegLocs = ceLocs(ceLocs < Lneg+eps)/Lneg;
    if length(intersect(jNegLocs,ceNegLocs)) < length(jNegLocs),
      warning('Cannot compute eta in negative electrode');
    else
      etaNegCeInd = find(ceLocs < Lneg+eps);
    end
  end
  [C,~,IB] = intersect(jPosLocs,csePosLocs);
  if length(C) < length(jPosLocs),
    warning('Cannot compute eta in positive electrode');
  else
    etaPosCseInd = IB;
    cePosLocs = 1 - (ceLocs(ceLocs > Lneg+Lsep-eps) - Lneg - Lsep)/Lpos;
    if length(intersect(round(1e6*jPosLocs),round(1e6*cePosLocs))) < length(jPosLocs),
      warning('Cannot compute eta in positive electrode');
    else
      etaPosCeInd = flipud(find(ceLocs > Lneg+Lsep-eps));
    end
  end
  
  %% Find intial electrode SOC 
  neg0 = cellData.neg.theta0; neg100 = cellData.neg.theta100;
  SOCNeg = SOC0*(neg100-neg0)+neg0;
  pos0 = cellData.pos.theta0; pos100 = cellData.pos.theta100;
  SOCPos = SOC0*(pos100-pos0)+pos0; 

  %%  Set up simulation
  A = ROM.A; B = ROM.B; C = ROM.C; D = ROM.D;

  numstates = size(A,1);
  numoutstates = size(C,1);
  x = zeros(duration,numstates); % from x0...xN
  thetaNeg = zeros(1,duration);  % from x0...xN
  thetaPos = zeros(1,duration);  % from x0...xN
  cellSOC = zeros(duration,1);   % from x0...xN

  out = zeros(duration,numoutstates); % from y1...yN
  x(1,:) = zeros(1,numstates); % initialize state
  thetaNeg(1) = SOCNeg; thetaPos(1) = SOCPos;
  cellSOC(1) = (SOCNeg - neg0)/(neg100 - neg0);
  
  % The following two gains are constant even if model blending
  integralCsNegGain = C(cseNegInd(1),end); 
  integralCsPosGain = C(csePosInd(1),end); 

  % If model blending for this sim, then load the appropriate model files
  if blend,
    [A2D,C2D,D2D,Temp,Z] = loadROMFiles( ROMDIR );
  end
    
  %% -----------------------------------------------------------------
  %  Main simulation loop
  %  -----------------------------------------------------------------
  for k = 0:duration-1, % k = time / deltaTime
    % Step 1: Update SOC values at this time step
    csNegAvg = integralCsNegGain*x(k+1,end) + csmaxNeg*SOCNeg;
    if csNegAvg < 0,
      warning('csNegAvg < 0'); csNegAvg = 0;
    end
    if csNegAvg > csmaxNeg,
      warning('cs_neg_avg > csmaxNeg'); csNegAvg = csmaxNeg;
    end
    csPosAvg = integralCsPosGain*x(k+1,end) + csmaxPos*SOCPos;
    if csPosAvg < 0,
      warning('csPosAvg < 0'); csPosAvg = 0;
    end
    if csPosAvg > csmaxPos,
      warning('csPosAvg > csmaxPos'); csPosAvg = csmaxPos;
    end
    thetaNeg(k+1) = csNegAvg/csmaxNeg; 
    thetaPos(k+1) = csPosAvg/csmaxPos;
    cellSOC(k+1)  = (thetaNeg(k+1) - neg0)/(neg100 - neg0);
    % Every 100 iterations, display something on output as progress update
    if (rem(k,100) == 0), fprintf('Iteration = %d, Cell SOC = %2.2f%%\n',...
        k,cellSOC(k+1)*100); end
    
    % Step 2: Get "C" and "D" terms to be used in output equation
    % If "dTerm" is "uniform", then we use closed-form algebraic equations
    % to compute D at this SOC and temperature using the dStr variables
    % stored in the ROM.  We need to first compute the inputs to the
    % functions stored in the strings in dStr.
    if strcmpi(dTerm,'uniform');
      Tfact       = (1/Tref - 1/Tk(k+1))/R;
      kappa       = kappa_ref(ce0) * exp(EactKappa * Tfact);
      sigma_neg   = cellData.neg.sigma_ref * exp(EactSigma * Tfact);
      sigma_pos   = cellData.pos.sigma_ref * exp(EactSigma * Tfact);
      kappa_eff_neg = kappa*(cellData.neg.eps_e)^brugKappaNeg;
      kappa_eff_sep = kappa*(cellData.sep.eps_e)^brugKappaSep; %#ok<NASGU>
      kappa_eff_pos = kappa*(cellData.pos.eps_e)^brugKappaPos;
      sigma_eff_neg = sigma_neg*(cellData.neg.eps_s)^brugSigmaNeg;
      sigma_eff_pos = sigma_pos*(cellData.pos.eps_s)^brugSigmaPos;

      jeq_neg     = kNeg*sqrt(ce0*(csmaxNeg-csNegAvg)*csNegAvg);
      jeq_pos     = kPos*sqrt(ce0*(csmaxPos-csPosAvg)*csPosAvg);

      % Linearizing the "D" matrix around a nonzero j(z,t)
      % See James Lee's PhD dissertation, chapter 6 for deviation from book
      javg_neg    = ik(k+1)/(asNeg*F*Lneg*Acell);
      javg_pos    = -ik(k+1)/(asPos*F*Lpos*Acell); 
      Rct_neg     = R*Tk(k+1)/(F^2*sqrt(jeq_neg^2+javg_neg^2/4));
      Rct_pos     = R*Tk(k+1)/(F^2*sqrt(jeq_pos^2+javg_pos^2/4));
      % Linearizing the "D" matrix around a zero j(z,t)
      nu_neg      = Lneg*sqrt(asNeg*(1/sigma_eff_neg+1/kappa_eff_neg)/...
                  (Rct_neg+RfilmNeg));    %#ok<NASGU>
      nu_pos      = Lpos*sqrt(asPos*(1/sigma_eff_pos+1/kappa_eff_pos)/...
                  (Rct_pos+RfilmPos));    %#ok<NASGU>
      D = [cellfun(@eval,cellData.tf.Dstr)]; %#ok<NBRAK>
    elseif blend, % else we are model blending pre-computed D matrices
      D = interpolate2D(Temp,Z,D2D,Tk(k+1),cellSOC(k+1));
    end % else we use fixed D from ROM.D
    if blend,
      C = interpolate2D(Temp,Z,C2D,Tk(k+1),cellSOC(k+1)); 
    end % else, we use fixed C from ROM.C

    % Step 3: Compute linear outputs
    %         y[k] = C[k]*x[k] + D[k]*u[k]
    out(k+1,:) = C*x(k+1,:)' + D*ik(k+1);

    % Step 4: Parse outputs into output variables, and add affine and
    %         nonlinear corrections, as appropriate
    outk = out(k+1,:);
    % Step 4a: Compute solid surface concentration 
    cseNeg(k+1,:) = outk(cseNegInd) + csmaxNeg*SOCNeg; 
    if any(cseNeg(k+1,:) < 0),
      warning('cseNeg < 0'); cseNeg(cseNeg < 0) = 0;
    end
    if any(cseNeg(k+1,:) > csmaxNeg),
      warning('cseNeg > csmax_neg'); cseNeg(cseNeg > csmaxNeg) = csmaxNeg;
    end
    % Technically part of cseNeg, but computed additionally for voltage eq.
    cse0(k+1,:) = outk(cse0Ind) + csmaxNeg*SOCNeg; % Don't duplicate 
    if any(cse0(k+1,:) < 0), cse0(cse0 < 0) = 0; end % warnings, though
    if any(cse0(k+1,:) > csmaxNeg), cse0(cse0 > csmaxNeg) = csmaxNeg; end
    csePos(k+1,:) = outk(csePosInd) + csmaxPos*SOCPos;
    if any(csePos(k+1,:) < 0),
      warning('cse_pos < 0'); csePos(csePos < 0) = 0;
    end
    if any(csePos(k+1,:) > csmaxPos),
      warning('cse_pos > csmax_pos'); csePos(csePos > csmaxPos) = csmaxPos;
    end
    % Technically part of csePos, but computed additionally for voltage eq.
    cseL(k+1,:) = outk(cseLInd) + csmaxPos*SOCPos; % Don't duplicate
    if any(cseL(k+1,:) < 0), cseL(cseL < 0) = 0; end % warnings, though
    if any(cseL(k+1,:) > csmaxPos), cseL(cseL > csmaxPos) = csmaxPos; end

    % Step 4b: Compute electrolyte concentration
    ce(k+1,:) = outk(ceInd) + ce0;
    gradce(k+1,:) = [0, outk(gradceInd), 0];
    if any(ce(k+1,:) < 0), warning('ce < 0'); ce(ce < 0) = eps; end

    % Step 4c: Compute solid-electrolyte potential difference
    uocpNegAvg = cellData.neg.Uocp{1}(thetaNeg(k+1));
    phiseNeg0(k+1,:) = outk(phiseNeg0Ind) + uocpNegAvg;
    
    % Step 4d: Compute phi_e
    phietilde1(k+1,:) = outk(phie1Ind);
    gradphietilde1(k+1,:) = [0, outk(gradphie1Ind), 0];
    phi_e_tilde_2 = 2*R*Tk(k+1)/F*(1-t_plus)*log(ce(k+1,:)/ce(k+1,1));    
    phie(k+1,:) = [0 phietilde1(k+1,:)] + phi_e_tilde_2 - ...
                         phiseNeg0(k+1,:);

    % Step 4e: Compute j
    jNeg(k+1,:) = outk(jNegInd); j0(k+1,:) = outk(j0Ind);
    jPos(k+1,:) = outk(jPosInd); jL(k+1,:) = outk(jLInd);
    
    % correct eta variable via asinh method
    j0NegCCActual = max(eps,kNeg*(csmaxNeg-cse0(k+1))^(1-alphaNeg)*...
                       cse0(k+1)^alphaNeg*ce(k+1,1)^(1-alphaNeg));
    j0NegActual = max(eps,kNeg*(csmaxNeg-cseNeg(k+1,etaNegCseInd)).^(1-alphaNeg).* ... 
                       cseNeg(k+1,etaNegCseInd).^alphaNeg.*ce(k+1,etaNegCeInd).^(1-alphaNeg));                     %#ok<FNDSB>
    etaNeg(k+1,:) = 2*R*Tk(k+1)/F*asinh(jNeg(k+1,:)./(2*j0NegActual)); 
    eta0(k+1,:)   = 2*R*Tk(k+1)/F*asinh(j0(k+1,:)/(2*j0NegCCActual));
    j0PosCCActual = max(eps,kPos*(csmaxPos-cseL(k+1))^(1-alphaPos).*...
                       cseL(k+1)^alphaPos*ce(k+1,end)^(1-alphaPos));
    j0PosActual = max(eps,kPos*(csmaxPos-csePos(k+1,etaPosCseInd)).^(1-alphaPos).* ... 
                       csePos(k+1,etaPosCseInd).^alphaPos.*ce(k+1,etaPosCeInd).^(1-alphaPos));                    
    etaPos(k+1,:) = 2*R*Tk(k+1)/F*asinh(jPos(k+1,:)./(2*j0PosActual)); 
    etaL(k+1,:)   = 2*R*Tk(k+1)/F*asinh(jL(k+1,:)/(2*j0PosCCActual));
        
    % Step 4g: Compute cell voltage
    Uocp0 = cellData.neg.Uocp{1}(cse0(k+1,:)/csmaxNeg);
    UocpL = cellData.pos.Uocp{1}(cseL(k+1,:)/csmaxPos); 
    vCell(k+1) = etaL(k+1,:) - eta0(k+1,:) + UocpL - Uocp0 + ...
               phietilde1(k+1,end) + phi_e_tilde_2(end) + ...
               F*(RfilmPos*jL(k+1,:) - RfilmNeg*j0(k+1,:));
             
    % Step 4h: Compute phi_s
    phisNeg(k+1,:) = outk(phisNegInd);
    phisPos(k+1,:) = outk(phisPosInd) + vCell(k+1);
    gradphisNeg(k+1,:) = outk(gradphisNegInd);
    gradphisPos(k+1,:) = outk(gradphisPosInd);
    
    % Step 5: Update state equation
    %         x[k] = A[k-1]*x[k-1] + B[k-1]*u[k-1]
    if blend,
      A = diag(interpolate2D(Temp,Z,A2D,Tk(k+1),cellSOC(k+1))); 
    end
    x(k+2,:) = A*x(k+1,:)' + B*ik(k+1); 
  end % run simulation loop
  fprintf('Complete: SOC = %2.2f%% and voltage = %2.4f V\n\n',...
    100*cellSOC(duration),vCell(duration));
      
  %% Return values as a structure
  cellData.neg = rmfield(cellData.neg,'Uocp');
  cellData.pos = rmfield(cellData.pos,'Uocp');
  
  ROMout.cellData = cellData;
  ROMout.time = 0:length(vCell)-1; 
  ROMout.time = ROM.draData.Tsamp*ROMout.time(:);
  ROMout.Iapp = ik;
  ROMout.Vcell = vCell;
  ROMout.SOC = cellSOC;
  
  ROMout.j_neg = jNeg;
  ROMout.eta_neg = etaNeg;
  ROMout.cse_neg = cseNeg;
  ROMout.phis_neg = phisNeg;
  ROMout.gradphis_neg = gradphisNeg;

  ROMout.j_pos = jPos;
  ROMout.eta_pos = etaPos;
  ROMout.cse_pos = csePos;
  ROMout.phis_pos = phisPos;
  ROMout.gradphis_pos = gradphisPos;
  
  ROMout.phie = phie;
  ROMout.gradphie1 = gradphietilde1;
  ROMout.ce = ce;
  ROMout.gradce = gradce;

  ROMout.locs.j_neg_locs = cellData.tf.val(jNegInd);
  ROMout.locs.eta_neg_locs = ROMout.locs.j_neg_locs;
  ROMout.locs.cse_neg_locs = cellData.tf.val(cseNegInd);
  ROMout.locs.phis_neg_locs = cellData.tf.val(phisNegInd);
  
  ROMout.locs.j_pos_locs = cellData.tf.val(jPosInd);
  ROMout.locs.eta_pos_locs = ROMout.locs.j_pos_locs;
  ROMout.locs.cse_pos_locs = cellData.tf.val(csePosInd);
  ROMout.locs.phis_pos_locs = cellData.tf.val(phisPosInd);

  ROMout.locs.phie_locs = [0; cellData.tf.val(phie1Ind)];
  ROMout.locs.phie_neg_locs = ROMout.locs.phie_locs;
  ROMout.locs.phie_neg_locs(ROMout.locs.phie_locs > ...
    cellData.neg.L + eps) = [];
  ROMout.locs.phie_sep_locs = ROMout.locs.phie_locs;
  ROMout.locs.phie_sep_locs(ROMout.locs.phie_sep_locs < ...
    cellData.neg.L - eps) = [];
  ROMout.locs.phie_sep_locs(ROMout.locs.phie_sep_locs > ...
    cellData.neg.L + cellData.sep.L + eps) = [];  
  ROMout.locs.phie_pos_locs = ROMout.locs.phie_locs;
  ROMout.locs.phie_pos_locs(ROMout.locs.phie_locs < ...
    cellData.neg.L + cellData.sep.L - eps) = [];  
  ROMout.locs.ce_locs = cellData.tf.val(ceInd);
  ROMout.locs.ce_neg_locs = ROMout.locs.ce_locs;
  ROMout.locs.ce_neg_locs(ROMout.locs.ce_locs > ...
    cellData.neg.L + eps) = [];
  ROMout.locs.ce_sep_locs = ROMout.locs.ce_locs;
  ROMout.locs.ce_sep_locs(ROMout.locs.ce_sep_locs < ...
    cellData.neg.L - eps) = [];
  ROMout.locs.ce_sep_locs(ROMout.locs.ce_sep_locs > ...
    cellData.neg.L + cellData.sep.L + eps) = [];      
  ROMout.locs.ce_pos_locs = ROMout.locs.ce_locs;
  ROMout.locs.ce_pos_locs(ROMout.locs.ce_locs < ...
    cellData.neg.L + cellData.sep.L - eps) = [];
  ROMout.cellData.const.init_SOC = SOC0;
end

%% Load DRA matrices for model-blended simulations
function [bigA,bigC,bigD,T,Z] = loadROMFiles( ROMDIR )
  files = dir( sprintf('%scell_*_ROM_*.mat',ROMDIR) );

  % Scan Directory for ROM files
  tAxis = []; socAxis = [];
  for ii = 1:length(files),
    ind  = find(files(ii).name == '_');
    bStr = files(ii).name(1:ind(3));
    eStr = files(ii).name(ind(5):end);

    var  = sscanf(files(ii).name,sprintf('%s%%dSOC_%%dC%s',bStr,eStr));
    SOC  = var(1);
    T    = var(2);

    tAxis   = [tAxis, T + 273.15]; %#ok<AGROW>
    socAxis = [socAxis, SOC];      %#ok<AGROW>
  end

  % Remove all of the duplicate data and sort
  T = unique(tAxis);
  Z = unique(socAxis);

  % Assemble big matrices
  for ii = 1:length( Z )
    for jj = 1:length( T )
      romName = sprintf('%s%s%dSOC_%dC%s',ROMDIR,bStr,...
                        Z(ii),T(jj)-273.15,eStr);
      if exist(romName,'file'),
        romResult = load(romName);
        fprintf('Loaded %s\n',romName);
        bigA(jj,ii,1,:) = diag(romResult.ROM.A); %#ok<AGROW>
        if max(abs(eig(romResult.ROM.A))) > 1,
          fprintf('Unstable ROM in %s\n',romName);
        end
        if min(eig(romResult.ROM.A)) < 0,
          fprintf('Oscilating ROM in %s\n',romName);
        end
        bigC(jj,ii,:,:) = romResult.ROM.C;       %#ok<AGROW>
        bigD(jj,ii,1,:) = romResult.ROM.D;       %#ok<AGROW>
      else
        error('Error: Missing %s in ROM directory',romName);
      end
    end
  end    
  Z = Z/100;
end

% Interpolate in two dimensions to find A,C matrices
function M = interpolate2D(T,Z,inMatrix,t,z)
  numT = length(T); numZ = length(Z); dt = 0; dz = 0;
  % Find lower and upper indices for interpolation in both dimensions
  Tlower_ind=sum(T<t);     Tlower_ind(Tlower_ind<1)=1;
  Tupper_ind=Tlower_ind+1; Tupper_ind(Tupper_ind>numT)=numT;
  Zlower_ind=sum(Z<z);     Zlower_ind(Zlower_ind<1)=1;
  Zupper_ind=Zlower_ind+1; Zupper_ind(Zupper_ind>numZ)=numZ;

  % Find fraction between indices for interpolation in both dimensions
  if numT>1 && Tupper_ind~=Tlower_ind
    deltaT = T(Tupper_ind)-T(Tlower_ind); dt=(t-T(Tlower_ind))/deltaT;
  end
  if numZ>1 && Zupper_ind~=Zlower_ind
    deltaZ = Z(Zupper_ind)-Z(Zlower_ind); dz=(z-Z(Zlower_ind))/deltaZ;
  end

  % Extract four matrices and perform the bilinear interpolation. Note,
  % using one "squeeze" instead of four is much faster in MATLAB
  M00 = inMatrix(Tlower_ind,Zlower_ind,:,:);
  M01 = inMatrix(Tupper_ind,Zlower_ind,:,:);
  MdT = inMatrix(Tlower_ind,Zupper_ind,:,:);
  M11 = inMatrix(Tupper_ind,Zupper_ind,:,:);
  M=squeeze((1-dz)*((1-dt)*M00+dt*M01)+dz*((1-dt)*MdT+dt*M11));
end