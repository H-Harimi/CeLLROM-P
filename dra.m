% function ROMs = dra(draData,tflist,cellData)
%
%   Returns reduced-order models created by the DRA corresponding to
%   tuning parameters in draData from transfer functions in tflist
%   populated with data from cellData
%
%   Inputs:
%     draData: DRA tuning parameters, loaded via readDRATable.m
%     tflist: Transfer functions to include, loaded via readDRATable.m
%     cellData: Cell parameter values, loaded via readParamTable.m
%
%   Output:
%     ROMs: Reduced-order models created via DRA (one ROM for every value
%           in the vector draData.order.)

% Copyright (c) 2015 by Gregory L. Plett of the University of Colorado 
% Colorado Springs (UCCS). This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0.
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L., "Battery
% Management Systems, Volume I, Battery Modeling," Artech House, 2015.

function ROMs = dra(draData,tflist,cellData)
  % Error checking:
  if (max(draData.hank1loc) + max(draData.hank2loc)) > ...
     (draData.Tlen/draData.Tsamp)-2,
    error('dra: hankLen must be <= Tlen/Tsamp - 2');
  end

  % Function flow: 
  % 1. Generate "s" vector of freqs
  % 2. Call transfer functions 
  % 3. Compute unit-pulse responses from frequency response
  % 4. Create combined Hankel matrix with these unit-pulse responses
  % 5. Use ERA to get SS model 
  % 6. Augment the SS model for integrator/sort
  
  % 1. Create "s" vector of frequencies
  Ts = 1/draData.Fs;
  Nfft = 2^(ceil(log2(draData.Fs*draData.Tlen)));
  f = 0:Nfft-1;
  s = (2j*draData.Fs)*tan(pi*f/Nfft); %#ok<NASGU>

  % 2. Call all transfer functions in tflist. For each one,
  %    a. Calculate IFFT
  %    b. Calculate downsampled unit-pulse response (save in pulseResp)
  %    c. save D term from transfer function in D
  %    d. save residue in Caug matrix
  pulseResp = []; D = []; Caug = []; dcGain = [];
  tFast = Ts*(0:Nfft-1); % time vector for IFFT
  tFinal = draData.Tsamp*(0:floor(tFast(end)/draData.Tsamp)); % final tvec
  lenTFinal = length(tFinal); % used frequently so calc once
  for theTF = 1:length(tflist)
    fprintf(tflist{theTF},'cellData'); fprintf('\n');
    % Call the transfer functions using "eval".
    % A "%s" in the transfer-function string is replaced with "cellData".
    [tf,Dterm,res0,cellData] = eval(sprintf(tflist{theTF},'cellData'));

    ifftTF = draData.Fs*real(ifft(tf.')); % approx impulse response
    stepTF = (cumsum(ifftTF)*Ts).';       % approx step response
    % Downsample step response by using interpolation at final points
    numRows = size(stepTF,1);
    sampTF = zeros(numRows,lenTFinal);
    discTF = zeros(numRows,lenTFinal);
    for theOut = 1:numRows
      sampTF(theOut,:) = interp1(tFast,stepTF(theOut,:),tFinal,'spline');
      discTF(theOut,:) = [0 diff(sampTF(theOut,:))]; % unit pulse resp
    end 
    % Build unit-pulse response matrix, D matrix and Caug matrix
    pulseResp = [pulseResp; discTF(:,2:end)]; %#ok<AGROW>
    D = [D; Dterm];                           %#ok<AGROW>
    Caug = [Caug; res0];                      %#ok<AGROW>
    dcGain = [dcGain; tf(:,1)];               %#ok<AGROW>
  end  
  
  % Here, we scale all unit-pulse responses to have the same "norm" before
  % we perform the "svd" operation.  Afterward, we then unscale the result
  % so that the final answer has the correct magnitude.  (Scaling helps
  % small-valued transfer functions compete well with large-valued transfer
  % functions in the optimization intrinsically performed by the svd.)
  normFact = sqrt((sum(pulseResp.^2,2))); % scaling factor...
  pulseResp = pulseResp./normFact(:,ones([1,size(pulseResp,2)]));
  pulseResp(isnan(pulseResp))=0;

  %  3. Generate the Hankel matrices from pulse response data
  fprintf('Building Hankel matrices: %s\n',datestr(now));
  sizePulseResp = size(pulseResp,1);
  j = draData.hank1loc; j(1)=0; t = draData.hank2loc; t(1)=0;
  hank0 = zeros(length(j)*sizePulseResp,length(t)); hank1 = hank0;
  for ind1 = 1:length(t)
    for ind2 = 1:length(j)
      hank0(sizePulseResp*(ind2-1)+1:sizePulseResp*ind2,ind1) = ...
            pulseResp(:,t(ind1)+j(ind2)+1);
      hank1(sizePulseResp*(ind2-1)+1:sizePulseResp*ind2,ind1) = ...
            pulseResp(:,t(ind1)+j(ind2)+2);
    end
  end

  %  4. Use ERA to get state-space realization (A,B,C) 
  m = 1; % dim of u, number of inputs
  p = size(pulseResp,1); % number of linear outputs
  fprintf('Starting SVDS: %s\n',datestr(now));
  [U,S,V] = svds(hank0,max(draData.order));
  fprintf('Complete SVDS: %s\n',datestr(now)');

  for theOrder = 1:length(draData.order)
    n = draData.order(theOrder); % size of state
    S1 = sqrtm(S(1:n,1:n)); Obs = U(:,1:n)*S1; Con = S1*V(:,1:n).';
    A = Obs\hank1/Con; eigA = eig(A);
    B = Con(:,1:m);
    C = Obs(1:p,:);    
    if eigA ~= conj(eigA),
      warning('Not all eigenvalues of A are real!\n');
    end
    if max(abs(eigA)) > 1,
      warning('DRA produced unstable ROM!\n');
    end
    if any(real(eigA)<0),
      warning('DRA produced oscillating ROM!\n');
    end
  
    %   5. Create final state space model by adding in an integrator 
    %      term. Do this only if Caug has a nonzero value. Otherwise, 
    %      none of the transfer functions had a pole at the origin 
    %      removed.
    C = C.*normFact(:,ones([1 size(C,2)])); % undo scaling factor
    if (Caug == 0)
      Afin = A; Bfin = B; Cfin = C;
    else 
      Afin = zeros(n+1,n+1); Afin(1,1) = 1; Afin(2:end,2:end) = A; 
      Bfin = [draData.Tsamp; B]; % Input factor for integrator state.
      Cfin = [Caug C];
    end
    
    % Final state-space realization: Convert to diagonal form, sort
    % eigenvalues in order, rescale the B and C matrices so that B contains
    % only units values.
    dra_sys = ss(Afin,Bfin,Cfin,D,draData.Tsamp); 
    sys_diag = canon(dra_sys,'modal',Inf);

    [pp,ind] = sort(diag(sys_diag.a));
    rom.A = diag(pp);
    rom.B = sys_diag.b(ind);
    rom.C = sys_diag.c(:,ind);
    rom.D = D;

    % Force B vector to be all ones; re-scale C matrix...
    for theRow = 1:size(sys_diag.c,1),
      rom.C(theRow,:) = rom.C(theRow,:).*rom.B';
    end
    rom.B = ones(size(rom.B));
    
    % It may make sense in some situations to uncomment the following
    % lines, which force the state-space model's dc gain to match the
    % desired dc gain from the transfer functions.  Your results may vary.
    % n = size(C,2);
    % actualGain = rom.C(:,1:n)*((eye(n)-rom.A(1:n,1:n))\rom.B(1:n,:));
    % wantedGain = dcGain - rom.D;
    % scaleFact  = wantedGain./actualGain;
    % scaleFact(wantedGain == 0) = 0;
    % rom.C = rom.C.*[scaleFact(:,ones(1,n)) ones(size(C,1),1)];

    % Save tuning parameters and computed ROM in structure of final results
    rom.Ts = draData.Tsamp;
    rom.order = n;
    rom.draData = draData;
    rom.cellData = cellData;
    rom.tflist = tflist;
    rom.sing_vals = diag(S(1:n,1:n));
    ROMs(theOrder) = rom; %#ok<AGROW>
  end
end