function [hfinal] = ant(varargin)

% ANT  Simulate Antarctic ice sheet flow using BUILDANT to
% extract data from re-gridded SeaRISE-Antarctic data.  See
% Ant50km.nc for metadata.
% Example: Using Ant50km.nc data for 40 ka run:
%   >> addpath('../')   % so Matlab/Octave can find other codes
%   >> ant
% Using different gridded data and different enhancement factor
% E=5.0, and saving volume time series for reload:
%   >> [t vol] = ant('Ant25km.nc',1,5.0);
%   >> save -v7 vol25km.mat t vol
% (Retrieve the data with "load vol25km.mat".)
% Run times of 9/8/2012 runs with Octave on bueler-lemur:
%   50 km = 373 sec,  25 km = 1829 sec,  20 km = 12065 sec.
% Calls:  BUILDANT, SIAGENERAL, DIFFUSION

addpath('../')

filename = 'Ant50km.nc';
doplot = 1;
E = 3; % default enhancement factor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% added input structure:
% inputs are added in pairs (order doesnt matter):

% 'trend' followed by 'hotter' or 'default' - the increasing temperature
% trend forces the temperature change at the rate of the global average
% rate ove the last century

% 'variable temperature' % followed by 'on' or 'off' implements the
% spatially variable temperature



% added sections: dealing with variable inputs
ind = strcmp(varargin, 'trend'); % followed by 'warmer', 'colder' or 'none'
if any(ind)
    trend = varargin{find(ind)+1}; 
else
    trend = 'default';
end

ind = strcmp(varargin, 'variable temperature'); % followed by 'on' or 'off'
if any(ind)
    varTemp = varargin{find(ind)+1};
else
    varTemp = 'off';
end

ind = strcmp(varargin, 'Change SMB'); % followed by 'yes' or 'no'
if any(ind)
    SMBchange = varargin{find(ind)+1};
else
    SMBchange = 'no';
end

ind = strcmp(varargin, 'temperature'); % followed by 'yes' or 'no'
if any(ind)
    temperature = varargin{find(ind)+1};
else
    temperature = 'default';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for model testing:
ind = strcmp(varargin, 'time blocks');
if any(ind)
    numTimeBlock = varargin{find(ind)+1}; 
else
    numTimeBlock = 'default';
end

ind = strcmp(varargin, 'grid resolution'); % followed by 'warmer', 'colder' or 'none'
if any(ind)
    gridResolution = varargin{find(ind)+1}; 
else
    gridResolution = 'default';
end

% make funbction that looks set the coefficients for the
% Temperature/elevation fit These will be used in the temporal
% iteration and updated accoding to changes in ice elevation
tempFit = makesurftemp();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set up model space

% read input data from NetCDF; no plot
[x,y,~,~,prcp,thk,topg,~] = buildant(0,filename);

% this was added to extrapolate data out of the available mapped ground
% topography domain
topg(topg == min(min(topg))) = nan;
topg = inpaint_nans(topg,4);


% grid info
Lx = (max(x) - min(x)) / 2;    Ly = (max(y) - min(y)) / 2;
J = length(x) - 1;    K = length(y) - 1;
dx = 2 * Lx / J;    dy = 2 * Ly / K;


%% print info in command prompt

fprintf('summary of input data:\n')
fprintf('  J = %d,  K = %d\n',J,K)
fprintf('  dx = %.3f km,  dy = %.3f km\n',dx/1000.0,dy/1000.0)
fprintf('  thickness     [min,max] = [%8.2f,%8.2f] m\n',...
  min(min(thk)),max(max(thk)))
fprintf('  bed elevation [min,max] = [%8.2f,%8.2f] m\n',...
  min(min(topg)),max(max(topg)))
fprintf('  precipitation [min,max] = [%8.5f,%8.5f] m a-1\n',...
  min(min(prcp)),max(max(prcp)))

%% set up time progression (time increments and total length)

% run-time and time-step (in years)
secpera = 31556926;
deltata = 1.0;

% for sensitivity testing
if ~strcmp(numTimeBlock,'default')
    NN = numTimeBlock;
else
    NN = 10;  % number of blocks of length tblock
end

tfa = 20000; % time scale over which the model runs
t = (linspace(0,tfa,NN)*secpera);
tblocka = tfa/NN;
fprintf('doing run of %.3f ka total, in blocks of %.3f a,\n',...
        tfa/1000.0,tblocka)
fprintf('  with max time-steps of %.3f a ...\n  running:\n',deltata)

% fix units and parameters for actual run
deltat = deltata * secpera;  % convert to seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make SMB a function of T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = prcp / secpera;  % alternatively: M = zeros(size(thk));

if any(any(thk<0)), error('negative input thicknesses detected'), end

%% create initial conditions

A = E * 1.0e-16 / secpera;

% solve SIA by doing blocks of tf time, and reporting volume
% robustness test
%[dim1,dim2] = size(thk)
%H = 2*thk;

 H = thk;

vol = zeros(1,NN);
rho = 910.0;    rhow = 1028.0;
hinit = getsurface(H,topg,rho,rhow);

%% run model

% set the the warming rate
if strcmp(varTemp,'on')
    if strcmp(trend,'hotter') % iteratively increase the temp every time increment
        % temperature forcing
        
        if strcmp(temperature,'default')
            %warmingRate = 0.74/100; %k/year based on last 100 years
            warmingRate = 20/20000;
        else
            warmingRate = temperature/20000;
        end
        
        warmingIncrement = warmingRate*tblocka;
    end
end

for k = 1:NN
    
    vol(k) = printvolume(k*tblocka,dx,dy,H);
    if strcmp(varTemp,'on')
        T = 273.15 + polyval(tempFit,topg+H); %K
        if strcmp(trend,'hotter') % iteratively increase the temp every time increment
            % temperature forcing
            T = T + k*warmingIncrement;
            % mean(mean(T)) % for debugging purposes
            if strcmp(SMBchange, 'yes')
                SMBtempChange = 0.06; % 6% change per celcius (RACMO)
                M = M * (1+k*warmingRate*SMBtempChange); % increasing/decreasing precip as a function of temp
            end
            
        elseif strcmp(trend,'colder') % cool down the deperature
            % need to do this
            disp('the colder functionality has not been included yet')
            break
            
        elseif strcmp(trend,'default') % keep going as if nothing was changed
            
        else
            disp('trend must be hotter or colder or default')
            
        end
        
        A = fluidity(T,H,rho);

        
    elseif strcmp(varTemp,'off')
    else
      disp('variable temperature must be on or off')
      break
    end
  
  [H,hfinal,~] = siageneral(Lx,Ly,J,K,H,deltat,tblocka*secpera, ...
                                topg,M,A);
  
  if any(any(H<0)), error('negative output thicknesses detected'), end
  
end

if doplot==0, return; end

%% makeplots

% % location for figures
 targetDir = 'C:\Users\kelian_dc\Documents\school\Masters thesis\Numerical Modelling\project\SIA approx project\';
 runName = 'hotteriteration';

figure(2)
imagesc(x/1000,y/1000,flipud(hinit),[0, 4000]), axis square, colorbar
xlabel('x  (km)','fontsize',14), ylabel('y  (km)','fontsize',14)
%title('initial surface elevation')
print([targetDir,'antinitial-',runName],'-dpng')

figure(3)
imagesc(x/1000,y/1000,flipud(hfinal),[0, 4000]), axis square, colorbar
xlabel('x  (km)','fontsize',14), ylabel('y  (km)','fontsize',14)
%title('final surface elevation')
print([targetDir,'antfinal-',runName],'-dpng')

figure(4)
imagesc(x/1000,y/1000,flipud(H-thk),[floor(min(min(flipud(H-thk)))/1000)*1000, floor(max(max(flipud(H-thk)))/1000)*1000]), ...
    axis square, colorbar
xlabel x, ylabel y, title('thickness change')
print([targetDir, 'antthickchange-', runName],'-dpng')

figure(5)
plot(t/secpera,vol/(1.0e6*1.0e9),'o-','markersize',11,'linewidth',2)
xlabel('t  (a)','fontsize',14)
ylabel('volume  (10^6 km^3)','fontsize',14)
grid on
hold on
print([targetDir,'antvol-',runName], '-dpng')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  function vol = printvolume(timea,dx,dy,thk)
    vol = sum(sum(thk)) * dx * dy;
    fprintf('  ice volume at time %7.3fka  is %.4e km^3\n',timea/1000.0,vol/1.0e9)
  end

  function h = getsurface(H,b,rho,rhow)
    % GETSURFACE  Helper function to build surface elevations carefully,
    % according to grounded or floating.
    f = rho / rhow;                     % fraction of floating ice below surface
    h = H + b;                          % only valid where H > 0 and grounded
    h(H <= 0) = max(b(H <= 0),0.0);     % if no ice
    floating = (H > 0) & (b < - f * H); % points where flotation criterion
    h(floating) = (1-f) * H(floating);  %   ... is applied
  end
  
  function A = fluidity(T,H,rho)
    
    dim             = size(T);
    
    % temperature gradient
    G = 0.01; % k/m
    
    % gravity
    g = 9.807; % kg/m^3
    
    %Pre-exponential constant    
    A0              = zeros(dim(1),dim(2));
    A0(T<=263.15)   = 3.985 * 10^-13; %s^(-1)Pa(-3)
    A0(T>263.15)    = 1.916 * 10^3; %s^(-1)Pa(-3)
    
    % Activation energy
    Q               = zeros(dim(1),dim(2));
    Q(T<=263.15)    = 60 * 10^3; %J/mol
    Q(T>263.15)     = 139 * 10^3; %J/mol
    
    % Universal gas constant
    R           = 8.314;
    
    % clausius clapeyrion constant
    beta        = 9.8*10^-8; 
        
    %find approx temperature at the middie of the glacier
    Tmiddle     = T+G*H/2;
    p           = rho*g*H/2;
    Tprime      = Tmiddle + beta*p;

    % for debugging:
%     figure
%     imagesc(Tprime)
    
    exponent = -Q./(R*Tprime);
    expTerm = exp(exponent);
    A           = A0.*expTerm;
   
    
    
%     figure(6)
%     imagesc(A)
%     colorbar
% 
%  yo  
 
  end
end
