function [trw,varargout] = VSLite_v2_3(syear,eyear,phi,T1,T2,M1,M2,T,P,varargin)
% VSLite_v2_3.m - Simulate tree ring width index given monthly climate inputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic Usage:
%    trw = VSLite_v2_3(syear,eyear,phi,T1,T2,M1,M2,T,P)
%    gives just simulated tree ring as ouput.
%
%   [trw,gT,gM,gE,M] = VSLite_v2_3(syear,eyear,phi,T1,T2,M1,M2,T,P))
%    also includes growth response to temperature, growth response to soil
%    moisture, scaled insolation index, and soil moisture estimate in outputs.
%
% Basic Inputs:
%   syear = start year of simulation.
%   eyear = end year of simulation.
%   phi = latitude of site (in degrees N)
%   T1 = scalar temperature threshold below which temp. growth response is zero (in deg. C)
%   T2 = scalar temperature threshold above which temp. growth response is one (in deg. C)
%   M1 = scalar soil moisture threshold below which moist. growth response is zero (in v/v)
%   M2 = scalar soil moisture threshold above which moist. growth response is one (in v/v)
%     (Note that optimal growth response parameters T1, T2, M1, M2 may be estimated
%      using code estimate_vslite_params_v2_3.m also freely available at
%      the NOAA NCDC Paleoclimatology software library.)
%   T = (12 x Nyrs) matrix of ordered mean monthly temperatures (in degEes C)
%   P = (12 x Nyrs) matrix of ordered accumulated monthly precipitation (in mm)
%
% Advanced Inputs (must be specified as property/value pairs):
%   'lbparams':  Parameters of the Leaky Bucket model of soil moisture.
%                These may be specified in an 8 x 1 vector in the following
%                order (otherwise the default values are read in):
%                   Mmax: scalar maximum soil moisture content (in v/v),
%                     default value is 0.76
%                   Mmin: scalar minimum soil moisture (in v/v), default
%                     value is 0.01
%                   alph: scalar runoff parameter 1 (in inverse months),
%                     default value is 0.093
%                   m_th: scalar runoff parameter 3 (unitless), default
%                     value is 4.886
%                   mu_th: scalar runoff parameter 2 (unitless), default
%                     value is 5.80
%                   rootd: scalar root/"bucket" depth (in mm), default
%                     value is 1000
%                   M0: initial value for previous month's soil moisture at
%                     t = 1 (in v/v), default value is 0.2
%                   substep: logical 1 or 0; perform monthly substepping in
%                     leaky bucket (1) or not (0)? Default value is 0.
%   'intwindow': Integration window. Which months' growth responses should
%                be intregrated to compute the annual ring-width index?
%                Specified as a 2 x 1 vector of integer values. Both
%                elements are given in integer number of months since January
%                (July) 1st of the current year in the Northern (Southern)
%                hemisphere, and specify the beginning and end of the integration
%                window, respectively. Defaults is [1 ; 12] (eg. integrate
%                response to climate over the corresponding calendar year,
%                assuming location is in the northern hemisphere).
%   'hydroclim': Value is a single character either taking value ['P'] or ['M'].
%                If ['M'], then 9th input is interpreted as an estimate of
%                soil moisture content (in v/v) rather than as precipitation.
%                Model default is to read in precipitation and use the CPC's
%                Leaky Bucket model of hydrology to estimate soil moisture,
%                however if soil moisture observations or a more sophisticated
%                estimate of moisture accounting for snow-related processes
%                is available, then using these data directly are recommended
%                (and will also speed up code).
%
% For more detailed documentation, see:
% 1) Tolwinski-Ward et al., An efficient forward model of the climate
% controls on interannual variation in tree-ring width, Climate Dynamics (2011)
% DOI: 10.1007/s00382-010-0945-5
%
% 2) Tolwinski-Ward et al., Erratum to: An efficient forward model of the climate
% controls on interannual variation in tree-ring width, Climate Dynamics (2011)
% DOI: 10.1007/s00382-011-1062-9
%
% 3) Tolwinski-Ward et al., Bayesian parameter estimation and
% interpretation for an intermediate model of tree-ring width, Clim. Past
% (2013), DOI: 10.5194/cp-9-1-2013
%
% 4) Documentation available with the model at http://www.ncdc.noaa.gov/paleo/softlib/softlib.html
%
% Revision History
% v0.1 - Original coding at monthly timestep from full daily timestep model (SETW, 4/09)
% v1.0 - Changed soil moisture module to the CPC Leaky Bucket model (SETW, 5/09)
% v1.1 - No upper parametric bounds for gT, gW as in full model; no density module (SETW, 9/09)
% v1.2 - Added adjustable integration window parameters (SETW, 1/10)
% v2.0 - Minor debugging for Octave compatibility, final version for publication (SETW, 10/10)
% v2.1 - Error in evapotranspiration calculation corrected (SETW, 7/11)
% v2.2 - Add switch to allow for monthly sub-stepping in soil moisture computation (SETW, N.Graham, K.Georgakaos, 9/11)
% v2.3 - Add switch to allow moisture M to be given as input rather than estimated
%        from T and P; add variable input options and improve commenting (SETW, 7/13)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iyear = syear:eyear;
nyrs = length(syear:eyear);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in advanced inputs if user-specified; else read in parameter defaults:
if nargin > 9 
    %%%% First fill parameter values in with defaults: %%%%%
    % Parameters of the Leaky Bucket model:
    Mmax = 0.76;
    Mmin = 0.01;
    alph = 0.093;
    m_th = 4.886;
    mu_th = 5.80;
    rootd = 1000;
    M0 = 0.2;
    substep = 0;
    % Integration window parameters:
    I_0 = 1;
    I_f = 12;
    % Hydroclimate variable:
    hydroclim = 'P';
    for i = 1:Nvararg/2
        namein = varargin{2*(i-1)+1};
        valin = varargin{2*i};
        switch namein
            case 'lbparams'
                Mmax = valin(1);
                Mmin = valin(2);
                alph = valin(3);
                m_th = valin(4);
                mu_th = valin(5);
                rootd = valin(6);
                M0 = valin(7);
                substep = valin(8);
            case 'intwindow'
                I_0 = valin(1);
                I_f = valin(2);
            case 'hydroclim'
                hydroclim = valin;
        end
    end
else % otherwise, read in defaults:
    % Parameters of the Leaky Bucket model:
    Mmax = 0.76;
    Mmin = 0.01;
    alph = 0.093;
    m_th = 4.886;
    mu_th = 5.80;
    rootd = 1000;
    M0 = 0.2;
    substep = 0;
    % Integration window parameters:
    I_0 = 1;
    I_f = 12;
    % Hydroclimate variable:
    hydroclim = 'P';
end
%%% Pre-allocate storage for outputs: %%%%
Gr = NaN(12,nyrs);
gT = NaN(12,nyrs);
gM = NaN(12,nyrs);
M = NaN(12,nyrs);
potEv = NaN(12,nyrs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load in or estimate soil moisture:
if strcmp(hydroclim,'M')
    % Read in soil moisture:
    M = P;
else
    % Compute soil moisture:
    if substep == 1;
        M = leakybucket_submonthly(syear,eyear,phi,T,P,Mmax,Mmin,alph,m_th,mu_th,rootd,M0);
    elseif substep == 0
        M = leakybucket_monthly(syear,eyear,phi,T,P,Mmax,Mmin,alph,m_th,mu_th,rootd,M0);
    elseif substep ~=1 && substep ~= 0
        disp('''substep'' must either be set to 1 or 0.');
        return
    end
end
% Compute gE, the scaled monthly proxy for insolation:
gE = Compute_gE(phi);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now compute growth responses to climate, and simulate proxy:
%%%%%%%%%%%%%%%%
% syear = start (first) year of simulation
% eyear = end (last) year of simulation
% cyear = year the model is currently working on
% iyear = index of simulation year
% Compute monthly growth response to T & M, and overall growth response G:
for cyear=1:length(iyear)      % begin cycling over years
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for t = 1:12  % begin cycling over months in a year
        %%% Calculate Growth Response functions gT(t) and gM(t)
        % First, temperature growth response:
        x = T(t,cyear);
        if (x < T1)
            gT(t,cyear) = 0;
        elseif (x >= T1) && (x <= T2)
            gT(t,cyear) = (x - T1)/(T2 - T1);
        elseif (x >= T2)
            gT(t,cyear) = 1;
        end
        % Next, Soil moisture growth response:
        x = M(t,cyear);
        if (x < M1)
            gM(t,cyear) = 0;
        elseif (x >= M1) && (x <= M2)
            gM(t,cyear) = (x - M1)/(M2 - M1);
        elseif (x >= M2)
            gM(t,cyear) = 1;
        end
    end % end month (t) cycle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute overall growth rate:
    Gr(:,cyear) = gE.*min(gT(:,cyear),gM(:,cyear));
end % end year cycle
%%%%%%%%%%%%%% Compute proxy quantity from growth responses %%%%%%%%%%%%%%%
width = NaN*ones(length(syear:eyear),1);
if phi>0 % if site is in the Northern Hemisphere:
    if I_0<0; % if we include part of the previous year in each year's modeled growth:
        startmo = 13+I_0;
        endmo = I_f;
        % use average of growth data across modeled years to estimate first year's growth due
        % to previous year:
        width(1) = sum(Gr(1:endmo,1)) + sum(mean(Gr(startmo:12,:),2));
        for cyear = 2:length(syear:eyear)
            width(cyear) = sum(Gr(startmo:12,cyear-1)) + sum(Gr(1:endmo,cyear));
        end
    else % no inclusion of last year's growth conditions in estimates of this year's growth:
        startmo = I_0+1;
        endmo = I_f;
        for cyear = 1:length(syear:eyear)
            width(cyear) = sum(Gr(startmo:endmo,cyear));
        end
    end
elseif phi<0 % if site is in the Southern Hemisphere:
    % (Note: in the Southern Hemisphere, ring widths are dated to the year in which growth began!)
    startmo = 7+I_0; % (eg. I_0 = -4 in SH corresponds to starting integration in March of cyear)
    endmo = I_f-6; % (eg. I_f = 12 in SH corresponds to ending integraion in June of next year)
    for cyear = 1:length(syear:eyear)-1
        width(cyear) = sum(Gr(startmo:12,cyear)) + sum(Gr(1:endmo,cyear+1));
    end
    % use average of growth data across modeled years to estimate last year's growth due
    % to the next year:
    width(length(syear:eyear)) = sum(Gr(startmo:12,length(syear:eyear)))+...
        sum(mean(Gr(1:endmo,:),2));
end
%
trw = ((width-mean(width))/std(width))'; % proxy series is standardized width.
%
if nargout >=1
    varargout(1) = {gT};
    varargout(2) = {gM};
    varargout(3) = {gE};
    varargout(4) = {M};
    varargout(5) = {potEv};
    varargout{6} = {width};
    varargout{7} = {mean(width)};
    varargout{8} = {std(width)};
end
%
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SUBROUTINES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% LEAKY BUCKET WITH SUBSTEPPING %%%%%%%%%%%%%%%%%%%%%%%%
function [M,potEv,ndl,cdays] = leakybucket_submonthly(syear,eyear,phi,T,P,...
    Mmax,Mmin,alph,m_th,mu_th,rootd,M0)
% leackybucket_submonthly.m - Simulate soil moisture; substeps within monthly timesteps
% to better capture nonlinearities and improve moisture estimates.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: [M,potEv,ndl,cdays] = leakybucket_submonthly(syear,eyear,phi,T,P,...
%                   Mmax,Mmin,alph,m_th,mu_th,rootd,M0)
%    outputs simulated soil moisture and potential evapotranspiration.
%
% Inputs:
%   syear = start year of simulation.
%   eyear = end year of simulation.
%   phi = latitude of site (in degrees N)
%   T = (12 x Nyrs) matrix of ordered mean monthly temperatures (in degEes C)
%   P = (12 x Nyrs) matrix of ordered accumulated monthly precipitation (in mm)
%   Mmax = scalar maximum soil moisture held by the soil (in v/v)
%   Mmin = scalar minimum soil moisture (for error-catching) (in v/v)
%   alph = scalar runoff parameter 1 (in inverse months)
%   m_th = scalar runoff parameter 3 (unitless)
%   mu_th = scalar runoff parameter 2 (unitless)
%   rootd = scalar root/"bucket" depth (in mm)
%   M0 = initial value for previous month's soil moisture at t = 1 (in v/v)
%
% Outputs:
%   M = soil moisture computed via the CPC Leaky Bucket model (in v/v, 12 x Nyrs)
%   potEv = potential evapotranspiration computed via Thornthwaite's 1947 scheme (in mm)
%
% SETW+ N. Graham and K. Georgakakos 2011

% modified by Nick G. and K. Georgakakos - to sub-step the monthly steps. Also this version has added
% soil moisture initial conditions for restarts, or spin-up.  Hands back monthly soil moisture
% and summer soil moisture as well - see varargout.  Nick G. 2011/06
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iyear = syear:eyear;
nyrs = length(iyear);
% Storage for growth response output variables (size [12 x Nyears]):
M   = NaN(12,nyrs);
potEv = NaN(12,nyrs);

% ADDED BY NICK
if(M0 < 0.)
    M0=200/rootd;
end

% Compute normalized daylength (neglecting small difference in calculation for leap-years)
latr = phi*pi/180;  % change to radians
ndays = [0 31 28 31 30 31 30 31 31 30 31 30 31];
cdays = cumsum(ndays);
sd = asin(sin(pi*23.5/180) * sin(pi * (((1:365) - 80)/180)))';   % solar declination
y = -tan(ones(365,1).* latr) .* tan(sd);
if ~isempty(find(y>=1,1))
    y(y>=1) = 1;
end
if ~isempty(find(y<=-1,1))
    y(y<=-1) = -1;
end
hdl = acos(y);
dtsi = (hdl.* sin(ones(365,1).*latr).*sin(sd))+(cos(ones(365,1).*latr).*cos(sd).*sin(hdl));
ndl=dtsi./max(dtsi); % normalized day length

% calculate mean monthly daylength (used for evapotranspiration in soil moisture calcs)
jday = cdays(1:12) +.5*ndays(2:13);
m_star = 1-tand(phi)*tand(23.439*cos(jday*pi/182.625));
mmm = NaN*ones(1,12);
for mo = 1:12
    if m_star(mo) < 0
        mmm(mo) = 0;
    elseif m_star(mo) >0 && m_star(mo) < 2
        mmm(mo) = m_star(mo);
    elseif m_star(mo) > 2
        mmm(mo) = 2;
    end
end
nhrs = 24*acosd(1-mmm)/180; % the number of hours in the day in the middle of the month
L = (ndays(2:13)/30).*(nhrs/12); % mean daylength in each month.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% -- year cycle -- %%%%
% syear = start (first) year of simulation
% eyear = end (last) year of simulation
% cyear = year the model is currently working on
% iyear = index of simulation year

for cyear=1:length(iyear)      % begin cycling over years
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for t = 1:12  % begin cycling over months in a year
        %%%%% Compute potential evapotranspiration for current month after Thornthwaite:
        if T(t,cyear) < 0
            Ep = 0;
        elseif T(t,cyear)>=0 && T(t,cyear) < 26.5;
            istar = (T(:,cyear)/5); istar(istar<0) = 0;
            I = sum(istar.^1.514);
            a = (6.75e-7)*I^3 - (7.71e-5)*I^2 + (1.79e-2)*I + .49;
            Ep = 16*L(t)*(10*T(t,cyear)/I)^a;
        elseif T(t,cyear) >= 26.5;
            Ep = -415.85 + 32.25*T(t,cyear) - .43* T(t,cyear)^2;
        end
        potEv(t,cyear) = Ep;
        %%%%% Now calculate soil moisture according to the CPC Leaky Bucket model
        %%%%% (see J. Huang et al, 1996). Set n-steps according to 2 mm increments
        %%%%% have to update alpha and Ep as well - 2 mm increments came from
        %%%%% testing by K. Georgakakos, but one could use 5 or more, with less "accurate" results.
        %%%%% Stepping is necessary because the parametization is linearized around init condition.
        %%%%%%%%%%%%%%%%%
        dp = 2.0; % mm of precip per increment
        nstep = floor(P(t,cyear)/dp)+1; % number of sub-monthly substeps
        Pinc = P(t,cyear)/nstep; % precip per substep
        alphinc = alph/nstep; % runoff rate per substep time interval
        Epinc = Ep/nstep; % potential evapotrans per substep.
        %%%%%%%%%%%%%%%%%
        % handling for sm_init
        if (t > 1)
            M0=M(t-1,cyear);
        elseif (t == 1 && cyear > 1)
            M0=M(12,cyear-1);
        end
        sm0=M0;
        
        for istep=1:nstep
            % evapotranspiration:
            Etrans = Epinc*sm0*rootd/(Mmax*rootd);
            % groundwater loss via percolation:
            G = mu_th*alphinc/(1+mu_th)*sm0*rootd;
            % runoff; contributions from surface flow (1st term) and subsurface (2nd term)
            R = Pinc*(sm0*rootd/(Mmax*rootd))^m_th + (alphinc/(1+mu_th))*sm0*rootd;
            dWdt = Pinc - Etrans - R - G;
            sm1 = sm0 + dWdt/rootd;
            %
            sm0=max(sm1,Mmin);
            sm0=min(sm0,Mmax);
        end
        M(t,cyear) = sm0;
        % error-catching:
        if M(t,cyear) <= Mmin; M(t,cyear) = Mmin; end;
        if M(t,cyear) >= Mmax; M(t,cyear) = Mmax; end;
        if isnan(M(t,cyear))==1; M(t,cyear) = Mmin; end;
    end % end month (t) cycle
end % end year cycle

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% LEAKY BUCKET WITHOUT SUBSTEPPING %%%%%%%%%%%%%%%%%%%%%
function [M,potEv,ndl,cdays] =...
    leakybucket_monthly(syear,eyear,phi,T,P,Mmax,Mmin,alph,m_th,mu_th,rootd,M0)
% leackybucket_monthly.m - Simulate soil moisture with coarse monthly time step.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: [M,potEv,ndl,cdays] = leakybucket_monthly(syear,eyear,phi,T,P,Mmax,Mmin,alph,m_th,mu_th,rootd,M0)
%    outputs simulated soil moisture and potential evapotranspiration.
%
% Inputs:
%   syear = start year of simulation.
%   eyear = end year of simulation.
%   phi = latitude of site (in degrees N)
%   T = (12 x Nyrs) matrix of ordered mean monthly temperatures (in degEes C)
%   P = (12 x Nyrs) matrix of ordered accumulated monthly precipitation (in mm)
%   Mmax = scalar maximum soil moisture held by the soil (in v/v)
%   Mmin = scalar minimum soil moisture (for error-catching) (in v/v)
%   alph = scalar runoff parameter 1 (in inverse months)
%   m_th = scalar runoff parameter 3 (unitless)
%   mu_th = scalar runoff parameter 2 (unitless)
%   rootd = scalar root/"bucket" depth (in mm)
%   M0 = initial value for previous month's soil moisture at t = 1 (in v/v)
%
% Outputs:
%   M = soil moisture computed via the CPC Leaky Bucket model (in v/v, 12 x Nyrs)
%   potEv = potential evapotranspiration computed via Thornthwaite's 1947 scheme (in mm)
%
% SETW 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iyear = syear:eyear;
nyrs = length(iyear);
% Storage for output variables (size [12 x Nyears]):
M  = NaN(12,nyrs);
potEv = NaN(12,nyrs);

% Compute normalized daylength (neglecting small difference in calculation for leap-years)
latr = phi*pi/180;  % change to radians
ndays = [0 31 28 31 30 31 30 31 31 30 31 30 31];
cdays = cumsum(ndays);
sd = asin(sin(pi*23.5/180) * sin(pi * (((1:365) - 80)/180)))';   % solar declination
y = -tan(ones(365,1).* latr) .* tan(sd);
if ~isempty(find(y>=1,1))
    y(y>=1) = 1;
end
if ~isempty(find(y<=-1,1))
    y(y<=-1) = -1;
end
hdl = acos(y);
dtsi = (hdl.* sin(ones(365,1).*latr).*sin(sd))+(cos(ones(365,1).*latr).*cos(sd).*sin(hdl));
ndl=dtsi./max(dtsi); % normalized day length

% calculate mean monthly daylength (used for evapotranspiration in soil moisture calcs)
jday = cdays(1:12) +.5*ndays(2:13);
m_star = 1-tand(phi)*tand(23.439*cos(jday*pi/182.625));
mmm = NaN*ones(1,12);
for mo = 1:12
    if m_star(mo) < 0
        mmm(mo) = 0;
    elseif m_star(mo) >0 && m_star(mo) < 2
        mmm(mo) = m_star(mo);
    elseif m_star(mo) > 2
        mmm(mo) = 2;
    end
end
nhrs = 24*acosd(1-mmm)/180; % the number of hours in the day in the middle of the month
L = (ndays(2:13)/30).*(nhrs/12); % mean daylength in each month.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% -- year cycle -- %%%%
% syear = start (first) year of simulation
% eyear = end (last) year of simulation
% cyear = year the model is currently working on
% iyear = index of simulation year

for cyear=1:nyrs     % begin cycling over years
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for t = 1:12  % begin cycling over months in a year
        %%%%% Compute potential evapotranspiration for current month after Thornthwaite:
        if T(t,cyear) < 0
            Ep = 0;
        elseif T(t,cyear)>=0 && T(t,cyear) < 26.5;
            istar = (T(:,cyear)/5); istar(istar<0) = 0;
            I = sum(istar.^1.514);
            a = (6.75e-7)*I^3 - (7.71e-5)*I^2 + (1.79e-2)*I + .49;
            Ep = 16*L(t)*(10*T(t,cyear)/I)^a;
        elseif T(t,cyear) >= 26.5;
            Ep = -415.85 + 32.25*T(t,cyear) - .43* T(t,cyear)^2;
        end
        potEv(t,cyear) = Ep;
        %%%%% Now calculate soil moisture according to the CPC Leaky Bucket model
        %%%%% (see J. Huang et al, 1996).
        if t > 1
            % evapotranspiration:
            Etrans = Ep*M(t-1,cyear)*rootd/(Mmax*rootd);
            % groundwater loss via percolation:
            G = mu_th*alph/(1+mu_th)*M(t-1,cyear)*rootd;
            % runoff; contributions from surface flow (1st term) and subsurface (2nd term)
            R = P(t,cyear)*(M(t-1,cyear)*rootd/(Mmax*rootd))^m_th +...
                (alph/(1+mu_th))*M(t-1,cyear)*rootd;
            dWdt = P(t,cyear) - Etrans - R - G;
            M(t,cyear) = M(t-1,cyear) + dWdt/rootd;
        elseif t == 1 && cyear > 1
            % evapotranspiration:
            Etrans = Ep*M(12,cyear-1)*rootd/(Mmax*rootd);
            % groundwater loss via percolation:
            G = mu_th*alph/(1+mu_th)*M(12,cyear-1)*rootd;
            % runoff; contributions from surface flow (1st term) and subsurface (2nd term)
            R = P(t,cyear)*(M(12,cyear-1)*rootd/(Mmax*rootd))^m_th +...
                (alph/(1+mu_th))*M(12,cyear-1)*rootd;
            dWdt = P(t,cyear) - Etrans - R - G;
            M(t,cyear) = M(12,cyear-1) + dWdt/rootd;
        elseif t == 1 && cyear == 1
            if M0 < 0; M0 = .20; end
            % evapotranspiration (take initial soil moisture value to be 200 mm)
            Etrans = Ep*M0*rootd/(Mmax*rootd);
            % groundwater loss via percolation:
            G = mu_th*alph/(1+mu_th)*(M0*rootd);
            % runoff; contributions from surface flow (1st term) and subsurface (2nd term)
            R = P(t,cyear)*(M0*rootd/(Mmax*rootd))^m_th + (alph/(1+mu_th))*M0*rootd;
            dWdt = P(t,cyear) - Etrans - R - G;
            M(t,cyear) = M0 + dWdt/rootd;
        end
        % error-catching:
        if M(t,cyear) <= Mmin; M(t,cyear) = Mmin; end;
        if M(t,cyear) >= Mmax; M(t,cyear) = Mmax; end;
        if isnan(M(t,cyear))==1; M(t,cyear) = Mmin; end;
    end % end month (t) cycle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end % end year cycle
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SCALED DAYLENGTH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gE] = Compute_gE(phi)
% Just what it sounds like... computes just gE from latitude a la VS-Lite,
% but without all the other stuff.
%
% Usage: gE = Compute_gE(phi)
%
% SETW 3/8/13

%
gE = NaN(12,1);
%
% Compute normalized daylength (neglecting small difference in calculation for leap-years)
latr = phi*pi/180;  % change to radians
ndays = [0 31 28 31 30 31 30 31 31 30 31 30 31];
cdays = cumsum(ndays);
sd = asin(sin(pi*23.5/180) * sin(pi * (((1:365) - 80)/180)))';   % solar declination
y = -tan(ones(365,1).* latr) .* tan(sd);
if ~isempty(find(y>=1,1))
    y(y>=1) = 1;
end
if ~isempty(find(y<=-1,1))
    y(y<=-1) = -1;
end
hdl = acos(y);
dtsi = (hdl.* sin(ones(365,1).*latr).*sin(sd))+(cos(ones(365,1).*latr).*cos(sd).*sin(hdl));
ndl=dtsi./max(dtsi); % normalized day length

% calculate mean monthly daylength (used for evapotranspiration in soil moisture calcs)
jday = cdays(1:12) +.5*ndays(2:13);
m_star = 1-tand(phi)*tand(23.439*cos(jday*pi/182.625));
mmm = NaN*ones(1,12);
for mo = 1:12
    if m_star(mo) < 0
        mmm(mo) = 0;
    elseif m_star(mo) >0 && m_star(mo) < 2
        mmm(mo) = m_star(mo);
    elseif m_star(mo) > 2
        mmm(mo) = 2;
    end
end
%nhrs = 24*acosd(1-mmm)/180; % the number of hours in the day in the middle of the month
%L = (ndays(2:13)/30).*(nhrs/12);
%
for t = 1:12
    gE(t) = mean(ndl(cdays(t)+1:cdays(t+1),1));
end
%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

