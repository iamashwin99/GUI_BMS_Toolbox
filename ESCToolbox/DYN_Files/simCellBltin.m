
function [vk,irk,hk,zk,sik,OCV] = simCellBltin(T,deltaT,model,z0,iR0,h0)
  % Force data to be column vector(s)
  ik = ik(:); iR0 = iR0(:);
  
  % Get model parameters from model structure
  RCfact = exp(-deltaT./abs(getParamESC('RCParam',T,model)))';
  G = getParamESC('GParam',T,model);
  Q = getParamESC('QParam',T,model);
  M = getParamESC('MParam',T,model);
  M0 = getParamESC('M0Param',T,model);
  RParam = getParamESC('RParam',T,model);
  R0Param = getParamESC('R0Param',T,model);
  etaParam = getParamESC('etaParam',T,model);
  
  etaik = ik; 
  etaik(ik<0) = etaParam*ik(ik<0); % compensate for coulombic efficiency
  
  % Simulate the dynamic states of the model
  if exist('ss','file'), % use control-system-toolbox method, if available
    sysd= ss(diag(RCfact),1-RCfact,eye(length(RCfact)),0,-1);
    irk = lsim(sysd,etaik,[],iR0);
  else
    irk=zeros([length(ik) length(iR0)]); irk(1,:) = iR0;
    for k = 2:length(ik),
      irk(k) = RCfact.*irk(k-1) + (1-RCfact)*etaik(k-1);
    end
  end
  zk = z0-cumsum([0;etaik(1:end-1)])*deltaT/(Q*3600); 
  if any(zk>1.1),
    warning('Current may have wrong sign as SOC > 110%');
  end
  
  % Hysteresis stuff
  hk=zeros([length(ik) 1]); hk(1) = h0; sik = 0*hk;
  fac=exp(-abs(G*etaik*deltaT/(3600*Q)));
  for k=2:length(ik),
    hk(k)=fac(k-1)*hk(k-1)-(1-fac(k-1))*sign(ik(k-1));
    sik(k) = sign(ik(k));
    if abs(ik(k))<Q/100, sik(k) = sik(k-1); end
  end
    
  % Compute output equation
  OCV = OCVfromSOCtemp(zk,T,model);
  vk = OCV - irk*RParam' - ik.*R0Param + M*hk + M0*sik;
  
