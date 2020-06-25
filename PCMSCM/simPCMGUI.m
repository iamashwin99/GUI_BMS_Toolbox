function [outputArg1] = simPCMGUI(Ns,Np,maxtime,t0,rt,randsoc,randcap,randres)
%SIMPCMGUI Summary of this function goes here
% simPCM: Simulate parallel-connected-module packs (cells are connected in
% parallel to make modules; these modules are connected in series to make
% packs). There are no inputs or direct outputs from this script.
%
% The parameters for each cell may be different (e.g., capacity,
% resistance, etc.) 

% Copyright (c) 2016 by Gregory L. Plett of 
% University of Colorado Colorado Springs (UCCS). 
%
% This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0
%
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.
% This file was modified by Ashwin Kumar for adapting to a GUI.
% This file is provided as a supplement to: Plett, Gregory L., "Battery
% Management Systems, Volume II, Equivalent-Circuit Methods," Artech House, 
% 2015.
 close all; clc;
outputArg1 = 0;
% Initialize some pack configuration parameters...
load E2model; % creates var. "model" with E2 cell parameter values
%Ns = 3;       % Number of modules connected in series to make a pack
%Np = 3;       % Number of cells connected in parallel in each module

% Initialize some simulation configuration parameters...
%maxtime = 4800; % Simulation run time in simulated seconds
%t0 = 2700; % Pack rests after time t0
storez = zeros([maxtime Ns Np]);  % create storage for SOC
storei = zeros([maxtime Ns Np]);  % create storage for current

% Initialize states for ESC cell model
z  = 0.25*ones(Ns,Np); %#ok<NASGU>
irc = zeros(Ns,Np);
h  = zeros(Ns,Np);

% Default initialization for cells within the pack
kvec = [0; maxtime+1]; % Iteration (time) vector for temp. profile
tvec = [25; 25]; % Default temperature profile
q    = getParamESC('QParam',25,model)*ones(Ns,Np); %#ok<NASGU>
rc   = exp(-1./abs(getParamESC('RCParam',25,model)))'*ones(Ns,Np);
r    = (getParamESC('RParam',25,model))';
m    = getParamESC('MParam',25,model)*ones(Ns,Np);
g    = getParamESC('GParam',25,model)*ones(Ns,Np);
r0   = getParamESC('R0Param',25,model)*ones(Ns,Np);  
%rt   = 0.000125; % 125 microOhm resistance for each tab

m = 0*m; % Eliminate model hysteresis for rough simulation: makes results 
         % easier to interpret. Then, can put hysteresis back via "m = m" 
         % for more realistic results.

% Modified initialization for cell variability...
% ------------------------------------------------------------------------
% Set individual random "initial SOC" values
if randsoc, % set to "if 1," to execute, or "if 0," to skip this code
  z=0.30+0.40*rand([Ns Np]); % rand. init. SOC for ea. cell
end

% Set individual random cell-capacity values
if randcap, % set to "if 1," to execute, or "if 0," to skip this code
  q=4.5+rand([Ns Np]);      % random capacity for ea. cell
end

% Set individual random cell-resistance relationships
if randres, % set to "if 1," to execute, or "if 0," to skip this code
  r0 = 0.005+0.020*rand(Ns,Np);
end
r0 = r0 + 2*rt; % add tab resistance to cell resistance

% Add faults to pack: cells faulted open- and short-circuit
% ------------------------------------------------------------------------
% To delete a PCM (open-circuit fault), set a resistance to Inf
%r0(1,1) = Inf; % for example...

% To delete a cell from a PCM (short-circuit fault), set its SOC to NaN
%z(1,2) = NaN; % for example, delete cell 2 in PCM 1 
Rsc = 0.0025; % Resistance value to use for cell whose SOC < 0%

% Get ready to simulate... first compute pack capacity in Ah
totalCap = min(sum(q,2)); % pack capacity = minimum module capacity
I = 10*totalCap; % cycle at 10C... not realistic, faster simulation


% Okay... now to simulate pack performance using ESC cell model.
for k = 1:maxtime,
  T = interp1(kvec,tvec,k); % cell temperature
  v = OCVfromSOCtemp(z,T,model); % get OCV for each cell: Ns * Np matrix
  v = v + m.*h - r.*irc; % add in capacitor voltages and hysteresis
  r0(isnan(z)) = Rsc; % short-circuit fault has "short-circuit" resistance
  V = (sum(v./r0,2) - I)./sum(1./r0,2);
  ik = (v-repmat(V,1,Np))./r0;
  z = z - (1/3600)*ik./q;  % Update each cell SOC
  z(z<0) = NaN; % set over-discharged cells to short-circuit fault
  irc = rc.*irc + (1-rc).*ik; % Update capacitor voltages
  fac = exp(-abs(g.*ik)./(3600*q));
  h = fac.*h + (1-fac).*sign(ik); % Update hysteresis voltages
  minz = min(z(:)); maxz = max(z(:)); % Check to see if SOC limit hit
  if minz < 0.05, I = -abs(I); end % stop discharging
  if maxz > 0.95, I = abs(I);  end % stop charging
  if k>t0, I = 0; end % rest 
  storez(k,:,:) = z; % Store SOC for later plotting
  storei(k,:,:) = ik; % store current for later plotting
end % for k

% Plot some sample results. In Figure 1, plot the individual SOC vs. time
% for all cells in all series PCMs. There is one subplot for each PCM.
nonocPCMs = ~isinf(sum(r0,2)); % modules that are not open-circuit faulted
figure(1); clf; t = (0:(length(storez(:,:,1))-1))/60; 
xplots = round(1.0*ceil(sqrt(Ns))); yplots = ceil(Ns/xplots); means = [];
for k = 1:Ns,
  zr=squeeze(100*storez(:,k,:));
  subplot(yplots,xplots,k); plot(t,zr); axis([0 ceil(maxtime/60) 0 100]);
  title(sprintf('Cells in PCM %d',k)); 
  ylabel('SOC (%)'); xlabel('Time (min)'); grid on
  zr(isnan(zr))=0; % exclude dead cells (failed short) from mean
  if nonocPCMs(k), % exclude from average if open-circuit!
    means = [means; mean(zr,2)']; %#ok<AGROW>
  end
end

% In Figure 2, plot the average SOC vs. time for all PCMs
figure(2); clf; plot(t,means'); grid on
xlabel('Time (min)'); ylabel('SOC (%)');
title('Average SOC for each PCM');
legendstrings = [];
for k=1:Ns,
  if nonocPCMs(k),
    legendstrings = [legendstrings; { sprintf('PCM %d',k) }]; %#ok<AGROW>
  end
end
legend(legendstrings);

% In Figure 3, plot the (maximum average SOC) minus (minimum average SOC)
figure(3); plot(t,(max(means)-min(means))); grid on
xlabel('Time (min)'); ylabel('Difference in SOC (%)');
title('Maximum-average SOC minus minimum-average SOC');

% In Figure 4, plot the individual cell current vs. time for all cells in 
% all series PCMs. There is one subplot for each PCM.
figure(4); clf; t = (0:(length(storei(:,:,1))-1))/60; 
for k = 1:Ns,
  zr=squeeze(storei(:,k,:));
  subplot(yplots,xplots,k); plot(t,zr); 
  axis([0 ceil(maxtime/60) -101 101]);
  title(sprintf('Cells in PCM %d',k)); 
  ylabel('Current (A)'); xlabel('Time (min)'); grid on 
end


end

