function [outputArg1] = simSCMGUI(Ns,Np,maxtime,t0,rt,randsoc,randcap,randres)
%SIMSCMGUI Summary of this function goes here
%   Detailed explanation goes here
outputArg1 = 0;

% simSCM: Simulate series-connected-module packs (cells are connected in
% series to make modules; these modules are connected in parallel to make
% packs). There are no inputs or direct outputs from this script.
%
% The parameters for each cell may be different (e.g., capacity,
% resistance, etc.) 

% Copyright (c) 2016 by Gregory L. Plett of 
% University of Colorado Colorado Springs (UCCS). 
% This file was modified by Ashwin Kumar for adapting to a GUI.
% This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0
%
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L., "Battery
% Management Systems, Volume II, Equivalent-Circuit Methods," Artech House, 
% 2015.

close all; clc;

% Initialize some pack configuration parameters...
load E2model; % creates var. "model" with E2 cell parameter values
%Ns = 8;       % Number of cells in series to make a module
%Np = 3;       % Number of modules connected in parallel to make pack

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
r0   = getParamESC('R0Param',25,model)*ones(Ns,Np); %#ok<NASGU>
%rt   = 0.000125; % 125 microOhm resistance for each tab

m = 0*m; % Eliminate model hysteresis for rough simulation: makes results 
         % easier to interpret. Then, can put hysteresis back via "m = m" 
         % for more realistic results.

% Modified initialization for cell variability
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
% To delete a SCM (open-circuit fault), set a resistance to Inf
%r0(1,1) = Inf; % for example...

% To delete a cell from a SCM (short-circuit fault), set its SOC to NaN
%z(1,2) = NaN; % for example, delete cell 2 in SCM 1 
Rsc = 0.0025; % Resistance value to use for cell whose SOC < 0%

% Get ready to simulate... first compute pack capacity in Ah
totalCap = min(sum(q,2)); % pack capacity = minimum module capacity
I = 10*totalCap; % cycle at 10C... not realistic, faster simulation

for k = 1:maxtime,
  v = OCVfromSOCtemp(z,25,model); % get OCV for each cell: Ns * Np matrix
  v = v + m.*h - r.*irc; % add in hysteresis and diffusion voltages
  r0(isnan(z)) = Rsc; % s-c fault has "short-circuit" resistance
  V = (sum(sum(v,1)./sum(r0,1),2)-I)./sum(1./sum(r0,1),2); % Bus V
  ik = (sum(v,1)-repmat(V,1,Np))./sum(r0,1); % 1 * Np cell currents
  ik = repmat(ik,Ns,1);
  z = z - (1/3600)*ik./q;  % Update each cell SOC
  z(z<0) = NaN; % set over-discharged cells to short-circuit fault
  irc = rc.*irc + (1-rc).*ik; % Update diffusion currents
  Ah = exp(-abs(g.*ik)./(3600*q));
  h = Ah.*h + (1-Ah).*sign(ik); % Update hysteresis voltages
  if min(z(:)) < 0.05, I = -abs(I); end % stop discharging
  if max(z(:)) > 0.95, I = abs(I);  end % stop charging
  if k>t0, I = 0; end % rest 
  storez(k,:,:) = z; % Store SOC for later plotting
  storei(k,:,:) = ik; % store current for later plotting
end % for k

% Plot some sample results. In Figure 1, plot the individual SOC vs. time
% for all cells in all series PCMs. There is one subplot for each SCM.
t = (0:(length(storez(:,:,1))-1))/60; 
xplots = round(1.0*ceil(sqrt(Np))); yplots = ceil(Np/xplots); 
figure(1); clf;
for k = 1:Np,
  subplot(yplots,xplots,k);
  zr=squeeze(100*storez(:,:,k));
  plot(t,zr); grid on;
  title(sprintf('SOC for cells in SCM %d',k)); 
  axis([0 ceil(maxtime/60) 0 100]);
  ylabel('SOC (%)'); xlabel('Time (min)'); 
  grid on
end

% In Figure 2, plot the individual cell current vs. time for all 
% cells in all series PCMs. There is one subplot for each SCM.
figure(2); clf;
for k = 1:Np,
  subplot(yplots,xplots,k);
  zr=squeeze(storei(:,:,k));
  plot(t,zr); grid on
  axis([0 ceil(maxtime/60) -60.1 60.1]);
  title(sprintf('Current for cells in SCM %d',k)); 
  ylabel('Current (A)'); xlabel('Time (min)'); 
  grid on
end

end

