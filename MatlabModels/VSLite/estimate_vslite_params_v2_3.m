function [T1,T2,M1,M2,varargout] = estimate_vslite_params_v2_3(T,P,phi,RW,varargin)
% Given calibration-interval temperature, precipitation, and ring-width data,
% and site latitude, estimate_vslite_params_v2_3.m performes a Bayesian parameter
% estimation of the growth response parameters T1, T2, M1 and M2. The
% default prior distributions are based on a survey of current literature
% pertaining to biological growth thresholds in trees; however uniform or
% user-defined four-parameter-beta distributions may also optionally be
% used as priors.  The scheme supports an assumption of either independent, 
% Gaussian model errors, or AR(1) error structure, and in either case the
% parameters of the error model may also be estimated.
% 
% For more details, see Tolwinski-Ward et al., 'Bayesian parameter estimation
% and interpretation for an intermediate model of tree-ring width', Clim. Past, 
% 9, 1-13, 3013, doi: 10.5194/cp-9-1-2013
%
% Basic Usage: [T1,T2,M1,M2] = estimate_vslite_params_v2_3(T,P,phi,RW)
%
% Basic Usage Inputs:
% T: monthly calibration-interval temperature, dimension 12 x number of calibration years.
% P: monthly calibration-interval precipitation, dimension 12 x number of cal. years.
% phi: site latitude in degrees N.
% RW: standardized annual calibration-interval ring-width index.
%
% Basic Usage Ouptuts:
% T1, T2, M1, M2: point estimates given by the median of their respective
%                 posterior distributions.
%
% Advanced Usage: [T1,T2,M1,M2,varargout] = vslite_bayes_param_cal(T,P,phi,RW,varargin)
%
% Advanced Optional Inputs:
%     Must be specified as property/value pairs.  Valid property/value pairs are:
%     'errormod'      Error model. Options are [0], [1], and [2] for white Gaussian
%                     noise, AR(1) model, or AR(2) model.  Default is [0].
%     'gparscalint'   Indices of years to use to estimate the growth response
%                     parameters T1, T2, M1, M2. Default is all years.
%     'eparscalint'   Indices of years to use to estimate the parameters of the
%                     error model if these parameters are to be estimated.
%                     Must be contiguous if using AR(1) error model. Default
%                     is all years. (Note: may underestimate error if not disjoint
%                     from interval used to fit growth response parameters
%                     as specified in 'gparscalint'.)
%     'errorpars'     Vector holding values of error model parameters is user
%                     wishes to fix their values rather than estimate them.
%                     For errormod == 0 (white noise model), values is a scalar
%                     with fixed value for sigma2w; for errormod == 1 (AR(1)),
%                     errorpars = [phi1 tau^2]; for errormod == 2 (AR(2)),
%                     errorpars = [phi1 phi2 tau^2]. No default (since default is
%                     to estimate parameters of a white noise error model).
%     'pt_ests'       Choices are ['mle'] or ['med'] to return either the
%                     posterior draw that maximizes the likelihood of the data
%                     or the marginal posterior medians as the parameter point
%                     estimates. Default is 'mle'.
%     'hydroclim'     Is the hydroclimate input variable (2nd basic input) 
%                     precipitation ['P'] or soil moisture ['M']? Default 
%                     is ['P']; CPC Leaky Bucket model is then used to estimate 
%                     M from input T and input P.
%     'substep'       If hydroclim == 'P', then 'substep' is logical 0/1
%                     depending on whether leaky bucket model without/with
%                     substepping is preferred.  Default is [0].
%     'intwindow'     VS-Lite integration window, specified as vector [I_0 I_f]
%                     Default is [0 12].
%     'nsamp'         200<=integer<=10,000 fixing number of MCMC iterations.
%                     Default is [1000].
%     'nbi'           Number of burn-in samples. Default is [200].
%     'nchain'        Integer number of comp threads for the computation
%                     of Rhat. Default is [3].
%     'gparpriors'    Form of growth parameter priors.  Either can be
%                     ['fourbet'] for a set of informative four-parameter 
%                     beta-distributed priors, or can be ['uniform']. 
%                     Parameterizations of these priors can be specified by
%                     the user, or literature-based choices will be used as
%                     defaults in either case (see following 5 property/value pairs
%                     for more information). Defaults to ['fourbet'].
%     'T1priorsupp'   2x1 vector with elements giving lower and upper bounds
%                     for support of uniform T1 prior. If not included in input
%                     argument list, default used is [0.0 8.5]
%     'T2priorsupp'   " T2 prior. Default is [9.0 20.0]
%     'M1priorsupp'   " M1 prior. Default is [0.01 0.03]
%     'M2priorsupp'   " M2 prior. Default is [0.1 0.5]
%     'fourbetparams' is a 4x4 matrix specifying parameters of the
%                     four-parameter beta distributed priors.  First row
%                     gives parameters of T1, second row gives T2, third
%                     row gives parameters for M1, and 4th row params of
%                     M2. Columns 1 and 2 give the two shape parameters,
%                     while columns 3 and 4 give the lower and upper bounds
%                     on the interval containing the transformed beta
%                     distribution. If not included in input arguent list,
%                     default parameter set based on current literature is
%                               [9   5   0   9
%                                3.5 3.5 10  24
%                                1.5 2.8 0.0 0.1
%                                1.5 2.5 0.1 0.5]
%                     (See Tolwinski-Ward et al., doi: 10.5194/cp-9-1-2013,
%                     for explanation of these choices.)
%     'convthresh'    Scalar value greater than 0.  Threshold for MCMC
%                     convergence; warning is displayed if abs(Rhat-1)>convthresh.
%                     Default value is [0.1].
%     'verbose'       Logical [0] or [1]; print progress to screen? Default 
%                     is [1].
%
% Advanced Optional Ouptuts (must be specified in the following order):
% T1dist, T2dist, M1dist, M2dist: Returns the entire numerical posterior distributions
%                 of the growth response parameters if the user wants to check for
%                 convergence, autocorrelation, multi-modality, etc., or to
%                 use the full distributions of the parameters to quantify
%                 uncertainty.
% Rhats:          Returns the convergence statistics associated with T1, T2, M1, M2,
%                 and sigma2rw if it was estimated.
% convwarning:    Logical [0] or [1] depending on whether any Rhat values were
%                 outside of the threshold distance from 1.
% -- Next ordered outputs for white noise error model (errormod==0): -- %%%
% sig2rw:         Point estimate of model error variance
% sigma2rwdist:   Returns the entire numerical posterior distribution
% Gdist:          Returns the numerical posterior distribution of monthly growth
%                 responses to the input climate for the corresponding posterior
%                 parameter sets; has dimension 12 x Nyrs x Nsamps.
%                 
% -- Next ordered outputs for AR(1) error model (errormod==1): -- %%%%%%%%%
% phi1:           Point estimate of AR(1) coefficient
% phi1dist:       Numerical posterior distribution of AR(1) coefficient 
% tau2:           Point estimate of error model innovation variance
% tau2dist:       Numerical distribution of error model innovation variance
% Gdist:          Returns the numerical posterior distribution of monthly growth
%                 responses to the input climate for the corresponding posterior
%                 parameter sets; has dimension 12 x Nyrs x Nsamps.
%
% SETW 9/20/2011:  version 1.0. Estimates T1, T2, M1, M2 for fixed sigma2w, assuming
%                  normally- and independent and identically distributed model residuals.
% SETW 4/15/2012:  version 2.0. Added simultaneous estimation of sigma2w under assumption
%                  of inverse-gamma prior, literature-based priors, Rhat convergence
%                  metric, and output plots.
% SETW 12/15/2012: version 2.1 for publication; added options for user-defined
%                  prior distributions.
% SETW 5/10/2013:  version 2.2: Revised data-level model and added options for white or
%                  AR(1) error models following reviewer comments in
%                  Climate of the Past Discussions; options for flexible
%                  calibration-intervals for growth parameters and error
%                  model parameters also added; MLE added as option for 
%                  point estimates;
%                  version 2.3: additional commenting added; option to
%                  condition on user-supplied input soil moisture data included 
%                  as opposed to necessarily estimating M from T & P via
%                  Leaky Bucket.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin > 4 % read in advanced options if user-specified:
    % first fill values in with defaults:
    errormod = 0;
    gparscalint = 1:length(RW);
    eparscalint = 1:length(RW);
    pt_ests = 'mle';
    hydroclim = 'P';
    substep = 0;
    intwindow = [0 12];
    nsamp = 1000;
    nbi = 200;
    nchain = 3;
    gparpriors = 'fourbet';
    aT1 = 9; bT1 = (aT1+1)/2;
    slpT1 = 9; intT1 = 0;
    aT2 = 3.5; bT2 = aT2;
    slpT2 = 14; intT2 = 10;
    slpM1 = 0.1; intM1 = 0.0;
    aM1 = 1.5; bM1 = (1-.035/slpM1)*aM1/(.035/slpM1); % .035 is target mean
    slpM2 = 0.4; intM2 = 0.1;
    aM2 = 1.5; bM2 = (1-(.25-intM2)/slpM2)*aM2/((.25-intM2)/slpM2); % .25 is target mean
    convthresh = .1;
    verbose = 1;
    % then over-write defaults if user-specified:
    Nvararg = length(varargin);
    for i = 1:Nvararg/2
        namein = varargin{2*(i-1)+1};
        valin = varargin{2*i};
        switch namein
            case 'errormod'
                errormod = valin;
            case 'gparscalint'
                gparscalint = valin;
            case 'eparscalint'
                eparscalint = valin;
            case 'errorpars'
                errorpars = valin;
            case 'pt_ests'
                pt_ests = valin;
            case 'hydroclim'
                hydroclim = valin;
            case 'substep'
                substep = valin;
            case 'intwindow'
                intwindow = valin;
            case 'nsamp'
                nsamp = valin;
            case 'nbi'
                nbi = valin;
            case 'nchain'
                nchain = valin;
            case 'gparpriors'
                gparpriors = valin;
                if strcmp(gparpriors,'uniform')
                    % read in default supports for uniform priors:
                    aT1 = 0; bT1 = 9;
                    aT2 = 10; bT2 = 24;
                    aM1 = 0; bM1 = .1;
                    aM2 = .1; bM2 = .5;
                end
            case 'T1priorsupp'
                aT1 = valin(1); bT1 = valin(2);
            case 'T2priorsupp'
                aT2 = valin(1); bT2 = valin(2);
            case 'M1priorsupp'
                aM1 = valin(1); bM1 = valin(2);
            case 'M2priorsupp'
                aM2 = valin(1); bM2 = valin(2);
            case 'fourbetparams'
                aT1 = valin(1,1); bT1 = valin(1,2);
                intT1 = valin(1,3); slpT1 = valin(1,4)-valin(1,3);
                aT2 = valin(2,1); bT2 = valin(2,2);
                intT2 = valin(2,3); slpT2 = valin(2,4)-valin(2,3);
                aM1 = valin(3,1); bM1 = valin(3,2);
                intM1 = valin(3,3); slpM1 = valin(3,4)-valin(3,3);
                aM2 = valin(4,1); bM2 = valin(4,2);
                intM2 = valin(4,3); slpM2 = valin(4,4)-valin(4,3);
            case 'convthresh'
                convthresh = valin;
            case 'verbose'
                verbose = valin;
        end
    end
else % otherwise, read in defaults:
    errormod = 0;
    gparscalint = 1:length(RW);
    eparscalint = 1:length(RW);
    pt_ests = 'mle';
    hydroclim = 'P';
    substep = 0;
    intwindow = [0 12];
    nsamp = 1000;
    nbi = 200;
    nchain = 3;
    gparpriors = 'fourbet';
    aT1 = 9; bT1 = (aT1+1)/2;
    slpT1 = 9; intT1 = 0;
    aT2 = 3.5; bT2 = aT2;
    slpT2 = 14; intT2 = 10;
    slpM1 = 0.1; intM1 = 0.0;
    aM1 = 1.5; bM1 = (1-.035/slpM1)*aM1/(.035/slpM1); % .035 is target mean
    slpM2 = 0.4; intM2 = 0.1;
    aM2 = 1.5; bM2 = (1-(.25-intM2)/slpM2)*aM2/((.25-intM2)/slpM2); % .25 is target mean
    convthresh = .1;
    verbose = 1;
end
%
% Take zscore of RW data to fulfill assumptions of model error/noise structure
RW = zscore(RW);
Gterms_dist = NaN(size(T));
%
% Compute soil moisture:
Mmax =.76; % maximum soil moisture; v/v
Mmin =.01; % minimum soil moisture; v/v
muth = 5.8; % mu from thornthwaite's Ep scheme
mth = 4.886; % m from thornthwaite's Ep scheme
alpha = .093;
Minit = 200; % initial value for soil moisture; v/v
dr = 1000; % root depth
%
Nyrs = size(T,2);
%
% Read in or compute estimate of soil moisture M:
if strcmp(hydroclim,'P') % if second input variable is precip,
    % then estimate soil moisture from T and P inputs via Leaky Bucket:
    if substep == 1;
        M = leakybucket_submonthly(1,Nyrs,phi,T,P,Mmax,Mmin,alpha,mth,muth,dr,Minit/dr);
    elseif substep == 0
        M = leakybucket_monthly(1,Nyrs,phi,T,P,Mmax,Mmin,alpha,mth,muth,dr,Minit/dr);
    end
elseif strcmp(hydroclim,'M') % if user supplies soil moisture estimate from elsewhere,
    % read in the soil moisture from the second input variable:
    M = P;
end
% Compute monthly growth response to insolation, gE:
gE = Compute_gE(phi);
%
%%%% Now do the MCMC sampling: %%%%%%%%%%%%%
for chain = 1:nchain
    % Storage space for realizations of parameters
    Tt = NaN(1,nsamp+nbi);
    To = NaN(1,nsamp+nbi);
    Mt = NaN(1,nsamp+nbi);
    Mo = NaN(1,nsamp+nbi);
    logLdata = NaN(1,nsamp+nbi);
    %
    if verbose; disp(['Working on chain ' num2str(chain)...
            ' out of ' num2str(nchain) '...']); end
    %
    % Initialize the MCMC:
    gT = NaN(size(T));
    gM = NaN(size(M));
    sim = 1;
    %
    % Initialize growth response parameters:
    if strcmp(gparpriors,'uniform')
        % Initialize Tt and To with draws from priors:
        Tt(sim) = unifrnd(aT1,bT1);
        To(sim) = unifrnd(aT2,bT2);
        % Initialize Mt and Mo with draws from priors:
        Mt(sim) = unifrnd(aM1,bM1);
        Mo(sim) = unifrnd(aM2,bM2);
    elseif strcmp(gparpriors,'fourbet')
        % Initialize Tt and To with draws from priors:
        Tt(sim) = slpT1*betarnd(aT1,bT1)+intT1;
        To(sim) = slpT2*betarnd(aT2,bT2)+intT2;
        % Initialize Mt and Mo with draws from priors:
        Mt(sim) = slpM1*betarnd(aM1,bM1)+intM1;
        Mo(sim) = slpM2*betarnd(aM2,bM2)+intM2;
    end
    %
    gT(T<Tt(sim)) = 0;
    gT(T>To(sim)) = 1;
    gT(T>Tt(sim)&T<To(sim)) = (T(T>Tt(sim)&T<To(sim))-Tt(sim))/(To(sim)-Tt(sim));
    %
    gM(M<Mt(sim)) = 0;
    gM(M>Mo(sim)) = 1;
    gM(M>Mt(sim)&M<Mo(sim)) = (M(M>Mt(sim)&M<Mo(sim))-Mt(sim))/(Mo(sim)-Mt(sim));
    %
    Gterms = min(gM,gT).*repmat(gE,1,size(T,2));
    %
    sim = sim+1;
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% The main sampling MCMC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create storage space to hold current values of the parameters to pass
    % to auxiliary sampling functions, and storage for MCMC realizations of these
    % parameters. Then initialize values of parameters.
    if errormod == 0;
        sigma2w = NaN(nsamp+nbi,1);
        sigma2w(1) = unifrnd(0,1); % initialize estimate from the prior.
        errorpars = sigma2w(1); % errorpars holds current value of error model parameters.
    elseif errormod == 1;
        tau2 = NaN(nsamp+nbi,1);
        phi1 = NaN(nsamp+nbi,1);
        % initialize estimates from the joint prior:
        phi1(1) = unifrnd(0,1);
        tau2(1) = unifrnd(0,1);
        while tau2(1) > 1-phi1(1)^2
            phi1(1) = unifrnd(0,1);
            tau2(1) = unifrnd(0,1);
        end
        % hold current values of error model parameters:
        errorpars(1) = phi1(1); errorpars(2) = tau2(1);
    end
    %
    while sim < nsamp+nbi+1
        %
        switch gparpriors
            case 'uniform'
                Tt(sim) = Tt_U_aux(Tt(sim-1),T,To(sim-1),gM,RW',errorpars,...
                    gE,Gterms,aT1,bT1,intwindow,gparscalint);
            case 'fourbet'
                Tt(sim) = Tt_lit_aux(Tt(sim-1),T,To(sim-1),gM,RW',errorpars,...
                    gE,Gterms,aT1,bT1,slpT1,intT1,intwindow,gparscalint);
        end
        gT(T<Tt(sim)) = 0;
        gT(T>To(sim-1)) = 1;
        gT(T>Tt(sim)&T<To(sim-1)) = (T(T>Tt(sim)&T<To(sim-1))-Tt(sim))/(To(sim-1)-Tt(sim));
        Gterms = min(gM,gT).*repmat(gE,1,size(T,2));
        %
        switch gparpriors
            case 'uniform'
                To(sim) = To_U_aux(To(sim-1),T,Tt(sim),gM,RW',errorpars,...
                    gE,Gterms,aT2,bT2,intwindow,gparscalint);
            case 'fourbet'
                To(sim) = To_lit_aux(To(sim-1),T,Tt(sim),gM,RW',...
                    errorpars,gE,Gterms,aT2,bT2,slpT2,intT2,intwindow,gparscalint);
        end
        gT(T<Tt(sim)) = 0;
        gT(T>To(sim)) = 1;
        gT(T>Tt(sim)&T<To(sim)) = (T(T>Tt(sim)&T<To(sim))-Tt(sim))/(To(sim)-Tt(sim));
        Gterms = min(gM,gT).*repmat(gE,1,size(T,2));
        %
        switch gparpriors
            case 'uniform'
                Mt(sim) = Mt_U_aux(Mt(sim-1),M,Mo(sim-1),gT,RW',errorpars,...
                    gE,Gterms,aM1,bM1,intwindow,gparscalint);
            case 'fourbet'
                Mt(sim) = Mt_lit_aux(Mt(sim-1),M,Mo(sim-1),gT,RW',...
                    errorpars,gE,Gterms,aM1,bM1,slpM1,intM1,intwindow,gparscalint);
        end
        gM(M<Mt(sim)) = 0;
        gM(M>Mo(sim-1)) = 1;
        gM(M>Mt(sim)&M<Mo(sim-1)) = (M(M>Mt(sim)&M<Mo(sim-1))-Mt(sim))/(Mo(sim-1)-Mt(sim));
        Gterms = min(gM,gT).*repmat(gE,1,size(T,2));
        %
        switch gparpriors
            case 'uniform'
                Mo(sim) = Mo_U_aux(Mo(sim-1),M,Mt(sim),gT,RW',errorpars,...
                    gE,Gterms,aM2,bM2,intwindow,gparscalint);
            case 'fourbet'
                Mo(sim) = Mo_lit_aux(Mo(sim-1),M,Mt(sim),gT,RW',errorpars,...
                    gE,Gterms,aM2,bM2,slpM2,intM2,intwindow,gparscalint);
        end
        gM(M<Mt(sim)) = 0;
        gM(M>Mo(sim)) = 1;
        gM(M>Mt(sim)&M<Mo(sim)) = (M(M>Mt(sim)&M<Mo(sim))-Mt(sim))/(Mo(sim)-Mt(sim));
        Gterms = min(gM,gT).*repmat(gE,1,size(T,2));
        %
        Gterms_dist(:,:,sim) = Gterms; 
        %
        % Now draw from error model parameters:
        if errormod == 0
            [errorpars,logLdata(sim)] = ...
                errormodel0_aux(errorpars,RW,Gterms,intwindow,eparscalint);
            sigma2w(sim) = errorpars;
        elseif errormod == 1
            [errorpars,logLdata(sim)] = ...
                errormodel1_aux(errorpars,RW,Gterms,intwindow,eparscalint);
            phi1(sim) = errorpars(1);
            tau2(sim) = errorpars(2);
        end
        %
        sim = sim+1;
    end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    eval(['Ttchain' num2str(chain) ' = Tt;']);
    eval(['Tochain' num2str(chain) ' = To;']);
    eval(['Mtchain' num2str(chain) ' = Mt;']);
    eval(['Mochain' num2str(chain) ' = Mo;']);
    eval(['Gtermsdist' num2str(chain) '= Gterms_dist;']);
    if errormod == 0
        eval(['sig2rwchain' num2str(chain) ' = sigma2w;']);
    elseif errormod == 1
        eval(['phi1chain' num2str(chain) ' = phi1;']);
        eval(['tau2chain' num2str(chain) ' = tau2;']);
    end
    eval(['logLchain' num2str(chain) ' = logLdata;']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POSTPROCESS SAMPLES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assess convergence:
Ttchains = []; Tochains = []; Mtchains = []; Mochains = [];
Ttensemb = []; Toensemb = []; Mtensemb = []; Moensemb = [];
Gterms_distensemb = [];
if errormod == 0
    sig2rwchains = [];
    sig2rwensemb = [];
elseif errormod == 1
    phi1chains = []; phi1ensemb = [];
    tau2chains = []; tau2ensemb = [];
end
logLchains = [];
logLensemb = [];
for i = 1:nchain
    if i < nchain
        eval(['Ttchains = [Ttchains ''Ttchain' num2str(i) ',''];'])
        eval(['Tochains = [Tochains ''Tochain' num2str(i) ',''];'])
        eval(['Mtchains = [Mtchains ''Mtchain' num2str(i) ',''];'])
        eval(['Mochains = [Mochains ''Mochain' num2str(i) ',''];'])
        if errormod == 0
            eval(['sig2rwchains = [sig2rwchains ''sig2rwchain' num2str(i) ',''];'])
        elseif errormod == 1
            eval(['phi1chains = [phi1chains ''phi1chain' num2str(i) ',''];'])
            eval(['tau2chains = [tau2chains ''tau2chain' num2str(i) ',''];'])
        end
        eval(['logLchains = [logLchains ''logLchain' num2str(i) ',''];'])
    else
        eval(['Ttchains = [Ttchains ''Ttchain' num2str(i) '''];'])
        eval(['Tochains = [Tochains ''Tochain' num2str(i) '''];'])
        eval(['Mtchains = [Mtchains ''Mtchain' num2str(i) '''];'])
        eval(['Mochains = [Mochains ''Mochain' num2str(i) '''];'])
        if errormod == 0
            eval(['sig2rwchains = [sig2rwchains ''sig2rwchain' num2str(i) '''];'])
        elseif errormod == 1
            eval(['phi1chains = [phi1chains ''phi1chain' num2str(i) '''];'])
            eval(['tau2chains = [tau2chains ''tau2chain' num2str(i) '''];'])
        end
        eval(['logLchains = [logLchains ''logLchain' num2str(i) '''];'])
    end
    eval(['Ttensemb = [Ttensemb Ttchain' num2str(i) '(nbi+1:end)];']);
    eval(['Toensemb = [Toensemb Tochain' num2str(i) '(nbi+1:end)];']);
    eval(['Mtensemb = [Mtensemb Mtchain' num2str(i) '(nbi+1:end)];']);
    eval(['Moensemb = [Moensemb Mochain' num2str(i) '(nbi+1:end)];']);
    eval(['Gterms_distensemb = cat(3,Gterms_distensemb,Gtermsdist'...
            num2str(i) '(:,:,nbi+1:end));']);
    if errormod == 0
        eval(['sig2rwensemb = [sig2rwensemb sig2rwchain' num2str(i) '(nbi+1:end)''];']);
    elseif errormod == 1
        eval(['phi1ensemb = [phi1ensemb phi1chain' num2str(i) '(nbi+1:end)''];']);
        eval(['tau2ensemb = [tau2ensemb tau2chain' num2str(i) '(nbi+1:end)''];']);
    end
    eval(['logLensemb = [logLensemb logLchain' num2str(i) '(nbi+1:end)];']);
end
%
if strcmp(pt_ests,'med');
    T1 = median(Ttensemb); T2 = median(Toensemb);
    M1 = median(Mtensemb); M2 = median(Moensemb);
    if errormod == 0
        sig2rw = median(sig2rwensemb);
    elseif errormod == 1
        phi1hat = median(phi1ensemb);
        tau2hat = median(tau2ensemb);
    end
elseif strcmp(pt_ests,'mle')
    mle_ind = find(logLensemb==max(logLensemb));
    if length(mle_ind)>1; mle_ind = mle_ind(1); end
    T1 = Ttensemb(mle_ind); T2 = Toensemb(mle_ind);
    M1 = Mtensemb(mle_ind); M2 = Moensemb(mle_ind);
    if errormod == 0
        sig2rw = sig2rwensemb(mle_ind);
    elseif errormod == 1
        phi1hat = phi1ensemb(mle_ind);
        tau2hat = tau2ensemb(mle_ind);
    end
end
%
eval(['RhatT1 = gelmanrubin92(nsamp,nbi,' Ttchains ');']);
eval(['RhatT2 = gelmanrubin92(nsamp,nbi,' Tochains ');']);
eval(['RhatM1 = gelmanrubin92(nsamp,nbi,' Mtchains ');']);
eval(['RhatM2 = gelmanrubin92(nsamp,nbi,' Mochains ');']);
if errormod == 0
    eval(['Rhatsig2rw = gelmanrubin92(nsamp,nbi,' sig2rwchains ');']);
elseif errormod == 1
    eval(['Rhatphi1 = gelmanrubin92(nsamp,nbi,' phi1chains ');']);
    eval(['Rhattau2 = gelmanrubin92(nsamp,nbi,' tau2chains ');']);
end
%
Rhats = [RhatT1 RhatT2 RhatM1 RhatM2];
if verbose == 1
    if errormod == 0
        Rhats = [Rhats Rhatsig2rw];
        disp('    Rhat for T1, T2, M1, M2, sigma2rw:');
        disp([RhatT1 RhatT2 RhatM1 RhatM2 Rhatsig2rw]);
    elseif errormod == 1
        Rhats = [Rhats Rhatphi1 Rhattau2];
        disp('    Rhat for T1, T2, M1, M2, phi1, tau2:');
        disp([RhatT1 RhatT2 RhatM1 RhatM2 Rhatphi1 Rhattau2]);
    end
end
if any(abs(Rhats-1)>convthresh)
    disp('Gelman and Rubin metric suggests MCMC has not yet converged to within desired threshold;')
    disp('Parameter estimation code should be re-run using a greater number of MCMC iterations.')
    disp('(See ''nsamp'' advanced input option.)')
    convwarning = 1;
else
    convwarning = 0;
end
%
if nargout > 0
    varargout{1} = Ttensemb; varargout{2} = Toensemb;
    varargout{3} = Mtensemb; varargout{4} = Moensemb;
    varargout{5} = Rhats; varargout{6} = convwarning;
    if errormod == 0
        varargout{7} = sig2rw; varargout{8} = sig2rwensemb;
        varargout{9} = Gterms_distensemb;
    elseif errormod == 1
        varargout{7} = phi1hat; varargout{8} = tau2hat;
        varargout{9} = phi1ensemb; varargout{10} = tau2ensemb;
        varargout{11} = Gterms_distensemb;
    end
end
%
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SOIL MOISTURE SUBROUTINES %%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%% CONDITIONAL PARAMETER SAMPLING SUBROUTINES %%%%%%%%%%%
function [Tt] = Tt_U_aux(Ttcurr,T,To,gM,RW,errorpars,gE,Gterms,att,btt,intwindow,cyrs)
% gM/gT is the matrix of gM for all the years and all the months of the previous simulation.
% T is the matrix of temperature for all years and all months of the previous simulation.
% att is the lower bound on the support of the uniform prior distribution for Tt
% btt is the upper bound on the support of the uniform prior distribution for Tt
%
% SETW 6/10/2010
Ny = size(Gterms,2);
I_0 = intwindow(1); I_f = intwindow(2);
%
if 1 % Sample from prior as proposal distribution!
    Ttprop = unifrnd(att,btt);
    gTprop = NaN*ones(12,Ny);
    gTprop(T<Ttprop) = 0;
    gTprop(T>To) = 1;
    gTprop(T<To&T>Ttprop) = (T(T<To&T>Ttprop)-Ttprop)/(To-Ttprop);
    gprop = diag(gE)*min(gM,gTprop);
    gcurr = Gterms;
    %
    %%%%%%%%%% account for variable integration window:
    if I_0<0; % if we include part of the previous year in each year's modeled growth:
        startmo = 13+I_0;
        endmo = I_f;
        prevseas = [mean(gprop(startmo:12,:),2) gprop(startmo:12,1:end-1)];
        gprop = gprop(1:endmo,:);
        gprop = [prevseas; gprop];
        prevseas = [mean(gcurr(startmo:12,:),2) gcurr(startmo:12,1:end-1)];
        gcurr = gcurr(1:endmo,:);
        gcurr = [prevseas; gcurr];
    else % no inclusion of last year's growth conditions in estimates of this year's growth:
        startmo = I_0+1;
        endmo = I_f;
        gprop = gprop(startmo:endmo,:);
        gcurr = gcurr(startmo:endmo,:);
    end
    %%%%%%%%%%%%
    %
    if length(errorpars) == 1 % White noise error model:
        sigma2rw = errorpars;
        expcurr = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr))).^2);
        expprop = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop))).^2);
        HR = exp(-.5*(expprop-expcurr)/sigma2rw);
    elseif length(errorpars) == 2 % AR(1) error model:
        phi1 = errorpars(1); tau2 = errorpars(2);
        sigma2rw = tau2/(1-phi1^2);
        %
        [iSig] = makeAR1covmat(phi1,tau2,length(cyrs));
        %
        Wcurr = ((sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr)))';
        Wprop = ((sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop)))';
        %
        logLprop = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wprop)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wprop);
        logLcurr = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wcurr)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wcurr);
        HR = exp(logLprop-logLcurr);
    end
end
% accept or reject the proposal.
if binornd(1,min(HR,1))==1
    Tt = Ttprop;
else
    Tt = Ttcurr;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Tt] = Tt_lit_aux(Ttcurr,T,To,gM,RW,errorpars,gE,Gterms,att,btt,slp,int,intwindow,cyrs)
% gM/gT is the matrix of gM for all the years and all the months of the previous simulation.
% T is the matrix of temperature for all years and all months of the previous simulation.
% att is the lower bound on the support of the uniform prior distribution for Tt
% btt is the upper bound on the support of the uniform prior distribution for Tt
%
% SETW 6/10/2010
Ny = size(Gterms,2);
I_0 = intwindow(1); I_f = intwindow(2);
%
if 1 % Sample from prior as proposal distribution!
    Ttprop = slp*betarnd(att,btt)+int;
    %     upperlim = normcdf(int+slp,Ttcurr,.5);
    %     lowerlim = normcdf(int,Ttcurr,.5);
    %     U = unifrnd(lowerlim,upperlim);
    %     Ttprop = norminv(U,Ttcurr,.5);
    %
    gTprop = NaN*ones(12,Ny);
    gTprop(T<Ttprop) = 0;
    gTprop(T>To) = 1;
    gTprop(T<To&T>Ttprop) = (T(T<To&T>Ttprop)-Ttprop)/(To-Ttprop);
    gprop = diag(gE)*min(gM,gTprop);
    gcurr = Gterms;
    %
    %%%%%%%%%% account for variable integration window:
    if I_0<0; % if we include part of the previous year in each year's modeled growth:
        startmo = 13+I_0;
        endmo = I_f;
        prevseas = [mean(gprop(startmo:12,:),2) gprop(startmo:12,1:end-1)];
        gprop = gprop(1:endmo,:);
        gprop = [prevseas; gprop];
        prevseas = [mean(gcurr(startmo:12,:),2) gcurr(startmo:12,1:end-1)];
        gcurr = gcurr(1:endmo,:);
        gcurr = [prevseas; gcurr];
    else % no inclusion of last year's growth conditions in estimates of this year's growth:
        startmo = I_0+1;
        endmo = I_f;
        gprop = gprop(startmo:endmo,:);
        gcurr = gcurr(startmo:endmo,:);
    end
    %%%%%%%%%%%%
    if length(errorpars) == 1 % White noise error model:
        sigma2rw = errorpars;
        expcurr = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr))).^2);
        expprop = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop))).^2);
        HR = exp(-.5*(expprop-expcurr)/sigma2rw);
        %
        %         % must also account for relative prior probs if using "geyer default" sampling:
        %         liklicurr = betapdf((Ttcurr-int)/slp,att,btt);
        %         likliprop = betapdf((Ttprop-int)/slp,att,btt);
        %         HR = (likliprop/liklicurr)*exp(-.5*(expprop-expcurr)/sigma2rw);
    elseif length(errorpars) == 2 % AR(1) error model:
        phi1 = errorpars(1); tau2 = errorpars(2);
        sigma2rw = tau2/(1-phi1^2);
        %
        [iSig] = makeAR1covmat(phi1,tau2,length(cyrs));
        %
        Wcurr = ((sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr)))';
        Wprop = ((sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop)))';
        %
        logLprop = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wprop)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wprop);
        logLcurr = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wcurr)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wcurr);
        HR = exp(logLprop-logLcurr);
        %
        %         % must also account for relative prior probs if using "geyer default" sampling:
        %         liklicurr = betapdf((Ttcurr-int)/slp,att,btt);
        %         likliprop = betapdf((Ttprop-int)/slp,att,btt);
        %         HR = (likliprop/liklicurr)*exp(logLprop-logLcurr);
    end
end
% accept or reject the proposal.
if binornd(1,min(HR,1))==1
    Tt = Ttprop;
else
    Tt = Ttcurr;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [To] = To_U_aux(Tocurr,T,Tt,gM,RW,errorpars,gE,Gterms,ato,bto,intwindow,cyrs)
% gM/gT is the matrix of gM for all the years and all the months of the previous simulation.
% T is the matrix of temperature for all years and all months of the previous simulation.
% att is the lower bound on the support of the uniform prior distribution for Tt
% btt is the upper bound on the support of the uniform prior distribution for Tt
%
% SETW 6/10/2010
Ny = size(Gterms,2);
I_0 = intwindow(1); I_f = intwindow(2);
if 1 % Sample from prior as proposal distribution!
    Toprop = unifrnd(ato,bto);
    gTprop = NaN*ones(12,Ny);
    gTprop(T<Tt) = 0;
    gTprop(T>Toprop) = 1;
    gTprop(T<Toprop&T>Tt) = (T(T<Toprop&T>Tt)-Tt)/(Toprop-Tt);
    gprop = diag(gE)*min(gM,gTprop);
    gcurr = Gterms;
    %
    %%%%%%%%%% account for variable integration window:
    if I_0<0; % if we include part of the previous year in each year's modeled growth:
        startmo = 13+I_0;
        endmo = I_f;
        prevseas = [mean(gprop(startmo:12,:),2) gprop(startmo:12,1:end-1)];
        gprop = gprop(1:endmo,:);
        gprop = [prevseas; gprop];
        prevseas = [mean(gcurr(startmo:12,:),2) gcurr(startmo:12,1:end-1)];
        gcurr = gcurr(1:endmo,:);
        gcurr = [prevseas; gcurr];
    else % no inclusion of last year's growth conditions in estimates of this year's growth:
        startmo = I_0+1;
        endmo = I_f;
        gprop = gprop(startmo:endmo,:);
        gcurr = gcurr(startmo:endmo,:);
    end
    %%%%%%%%%%%%
    %
    if length(errorpars) == 1 % White noise error model:
        sigma2rw = errorpars;
        expcurr = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr))).^2);
        expprop = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop))).^2);
        HR = exp(-.5*(expprop-expcurr)/sigma2rw);
    elseif length(errorpars) == 2 % AR(1) error model:
        phi1 = errorpars(1); tau2 = errorpars(2);
        sigma2rw = tau2/(1-phi1^2);
        %
        [iSig] = makeAR1covmat(phi1,tau2,length(cyrs));
        %
        Wcurr = ((sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr)))';
        Wprop = ((sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop)))';
        %
        logLprop = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wprop)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wprop);
        logLcurr = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wcurr)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wcurr);
        HR = exp(logLprop-logLcurr);
    end
end
% accept or reject the proposal.
if binornd(1,min(HR,1))==1
    To = Toprop;
else
    To = Tocurr;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [To] = To_lit_aux(Tocurr,T,Tt,gM,RW,errorpars,gE,Gterms,ato,bto,slp,int,intwindow,cyrs)
% gM/gT is the matrix of gM for all the years and all the months of the previous simulation.
% T is the matrix of temperature for all years and all months of the previous simulation.
% att is the lower bound on the support of the uniform prior distribution for Tt
% btt is the upper bound on the support of the uniform prior distribution for Tt
%
% SETW 6/10/2010
Ny = size(Gterms,2);
I_0 = intwindow(1); I_f = intwindow(2);
if 1 % Sample from prior as proposal distribution!
    Toprop = slp*betarnd(ato,bto)+int;
    %     upperlim = normcdf(int+slp,Tocurr,.5);
    %     lowerlim = normcdf(int,Tocurr,.5);
    %     U = unifrnd(lowerlim,upperlim);
    %     Toprop = norminv(U,Tocurr,.5);
    %
    gTprop = NaN*ones(12,Ny);
    gTprop(T<Tt) = 0;
    gTprop(T>Toprop) = 1;
    gTprop(T<Toprop&T>Tt) = (T(T<Toprop&T>Tt)-Tt)/(Toprop-Tt);
    gprop = diag(gE)*min(gM,gTprop);
    gcurr = Gterms;
    %
    %%%%%%%%%% account for variable integration window:
    if I_0<0; % if we include part of the previous year in each year's modeled growth:
        startmo = 13+I_0;
        endmo = I_f;
        prevseas = [mean(gprop(startmo:12,:),2) gprop(startmo:12,1:end-1)];
        gprop = gprop(1:endmo,:);
        gprop = [prevseas; gprop];
        prevseas = [mean(gcurr(startmo:12,:),2) gcurr(startmo:12,1:end-1)];
        gcurr = gcurr(1:endmo,:);
        gcurr = [prevseas; gcurr];
    else % no inclusion of last year's growth conditions in estimates of this year's growth:
        startmo = I_0+1;
        endmo = I_f;
        gprop = gprop(startmo:endmo,:);
        gcurr = gcurr(startmo:endmo,:);
    end
    %%%%%%%%%%%%
    %
    if length(errorpars) == 1 % White noise error model:
        sigma2rw = errorpars;
        expcurr = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr))).^2);
        expprop = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop))).^2);
        HR = exp(-.5*(expprop-expcurr)/sigma2rw);
        %
        %         % must also account for relative prior probs if using "geyer default" sampling:
        %         liklicurr = betapdf((Tocurr-int)/slp,ato,bto);
        %         likliprop = betapdf((Toprop-int)/slp,ato,bto);
        %         HR = (likliprop/liklicurr)*exp(-.5*(expprop-expcurr)/sigma2rw);
    elseif length(errorpars) == 2 % AR(1) error model:
        phi1 = errorpars(1); tau2 = errorpars(2);
        sigma2rw = tau2/(1-phi1^2);
        %
        [iSig] = makeAR1covmat(phi1,tau2,length(cyrs));
        %
        Wcurr = ((sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr)))';
        Wprop = ((sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop)))';
        %
        logLprop = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wprop)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wprop);
        logLcurr = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wcurr)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wcurr);
        HR = exp(logLprop-logLcurr);
        %
        %         % must also account for relative prior probs if using "geyer default" sampling:
        %         liklicurr = betapdf((Ttcurr-int)/slp,ato,bto);
        %         likliprop = betapdf((Ttprop-int)/slp,ato,bto);
        %         HR = (likliprop/liklicurr)*exp(logLprop-logLcurr);
    end
end
% accept or reject the proposal.
if binornd(1,min(HR,1))==1
    To = Toprop;
else
    To = Tocurr;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Mt] = Mt_U_aux(Mtcurr,M,Mo,gT,RW,errorpars,gE,Gterms,amt,bmt,intwindow,cyrs)
% gM/gT is the matrix of gM for all the years and all the months of the previous simulation.
% M is the matrix of soil moisture for all years and all months of the previous simulation.
% amt is the lower bound on the support of the uniform prior distribution for Mt
% bmt is the upper bound on the support of the uniform prior distribution for Mt
%
% SETW 6/10/2010
Ny = size(Gterms,2);
I_0 = intwindow(1); I_f = intwindow(2);
if 1 % Sample from prior as proposal distribution!
    Mtprop = unifrnd(amt,bmt);
    gWprop = NaN*ones(12,Ny);
    gWprop(M<Mtprop) = 0;
    gWprop(M>Mo) = 1;
    gWprop(M<Mo&M>Mtprop) = (M(M<Mo&M>Mtprop)-Mtprop)/(Mo-Mtprop);
    gprop = diag(gE)*min(gWprop,gT);
    gcurr = Gterms;
    %
    %%%%%%%%%% account for variable integration window:
    if I_0<0; % if we include part of the previous year in each year's modeled growth:
        startmo = 13+I_0;
        endmo = I_f;
        prevseas = [mean(gprop(startmo:12,:),2) gprop(startmo:12,1:end-1)];
        gprop = gprop(1:endmo,:);
        gprop = [prevseas; gprop];
        prevseas = [mean(gcurr(startmo:12,:),2) gcurr(startmo:12,1:end-1)];
        gcurr = gcurr(1:endmo,:);
        gcurr = [prevseas; gcurr];
    else % no inclusion of last year's growth conditions in estimates of this year's growth:
        startmo = I_0+1;
        endmo = I_f;
        gprop = gprop(startmo:endmo,:);
        gcurr = gcurr(startmo:endmo,:);
    end
    %%%%%%%%%%%%
    %
    if length(errorpars) == 1 % White noise error model:
        sigma2rw = errorpars;
        expcurr = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr))).^2);
        expprop = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop))).^2);
        HR = exp(-.5*(expprop-expcurr)/sigma2rw);
    elseif length(errorpars) == 2 % AR(1) error model:
        phi1 = errorpars(1); tau2 = errorpars(2);
        sigma2rw = tau2/(1-phi1^2);
        %
        [iSig] = makeAR1covmat(phi1,tau2,length(cyrs));
        %
        Wcurr = ((sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr)))';
        Wprop = ((sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop)))';
        %
        logLprop = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wprop)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wprop);
        logLcurr = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wcurr)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wcurr);
        HR = exp(logLprop-logLcurr);
    end
end
% accept or reject the proposal.
if binornd(1,min(HR,1))==1
    Mt = Mtprop;
else
    Mt = Mtcurr;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Mt] = Mt_lit_aux(Mtcurr,M,Mo,gT,RW,errorpars,gE,Gterms,amt,bmt,slp,int,intwindow,cyrs)
% gM/gT is the matrix of gM for all the years and all the months of the previous simulation.
% M is the matrix of soil moisture for all years and all months of the previous simulation.
% amt is the lower bound on the support of the uniform prior distribution for Mt
% bmt is the upper bound on the support of the uniform prior distribution for Mt
%
% SETW 6/10/2010
Ny = size(Gterms,2);
I_0 = intwindow(1); I_f = intwindow(2);
if 1 % Sample from prior as proposal distribution!
    Mtprop = slp*betarnd(amt,bmt)+int;
    %     upperlim = normcdf(int+slp,Mtcurr,.01);
    %     lowerlim = normcdf(int,Mtcurr,.01);
    %     U = unifrnd(lowerlim,upperlim);
    %     Mtprop = norminv(U,Mtcurr,.01);
    %
    gWprop = NaN*ones(12,Ny);
    gWprop(M<Mtprop) = 0;
    gWprop(M>Mo) = 1;
    gWprop(M<Mo&M>Mtprop) = (M(M<Mo&M>Mtprop)-Mtprop)/(Mo-Mtprop);
    gprop = diag(gE)*min(gWprop,gT);
    gcurr = Gterms;
    %
    %%%%%%%%%% account for variable integration window:
    if I_0<0; % if we include part of the previous year in each year's modeled growth:
        startmo = 13+I_0;
        endmo = I_f;
        prevseas = [mean(gprop(startmo:12,:),2) gprop(startmo:12,1:end-1)];
        gprop = gprop(1:endmo,:);
        gprop = [prevseas; gprop];
        prevseas = [mean(gcurr(startmo:12,:),2) gcurr(startmo:12,1:end-1)];
        gcurr = gcurr(1:endmo,:);
        gcurr = [prevseas; gcurr];
    else % no inclusion of last year's growth conditions in estimates of this year's growth:
        startmo = I_0+1;
        endmo = I_f;
        gprop = gprop(startmo:endmo,:);
        gcurr = gcurr(startmo:endmo,:);
    end
    %%%%%%%%%%%%
    %
    if length(errorpars) == 1 % White noise error model:
        sigma2rw = errorpars;
        expcurr = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr))).^2);
        expprop = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop))).^2);
        HR = exp(-.5*(expprop-expcurr)/sigma2rw);
        %
        %         % must also account for relative prior probs if using "geyer default" sampling:
        %         liklicurr = betapdf((Mtcurr-int)/slp,amt,bmt);
        %         likliprop = betapdf((Mtprop-int)/slp,amt,bmt);
        %         HR = (likliprop/liklicurr)*exp(-.5*(expprop-expcurr)/sigma2rw);
    elseif length(errorpars) == 2 % AR(1) error model:
        phi1 = errorpars(1); tau2 = errorpars(2);
        sigma2rw = tau2/(1-phi1^2);
        %
        [iSig] = makeAR1covmat(phi1,tau2,length(cyrs));
        %
        Wcurr = ((sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr)))';
        Wprop = ((sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop)))';
        %
        logLprop = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wprop)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wprop);
        logLcurr = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wcurr)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wcurr);
        HR = exp(logLprop-logLcurr);
        %
        %         % must also account for relative prior probs if using "geyer default" sampling:
        %         liklicurr = betapdf((Ttcurr-int)/slp,amt,bmt);
        %         likliprop = betapdf((Ttprop-int)/slp,amt,bmt);
        %         HR = (likliprop/liklicurr)*exp(logLprop-logLcurr);
    end
end
% accept or reject the proposal.
if binornd(1,min(HR,1))==1
    Mt = Mtprop;
else
    Mt = Mtcurr;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Mo] = Mo_U_aux(Mocurr,M,Mt,gT,RW,errorpars,gE,Gterms,amo,bmo,intwindow,cyrs)
% gM/gT is the matrix of gM for all the years and all the months of the previous simulation.
% M is the matrix of soil moisture for all years and all months of the previous simulation.
% amt is the lower bound on the support of the uniform prior distribution for Mt
% bmt is the upper bound on the support of the uniform prior distribution for Mt
%
% SETW 6/10/2010
Ny = size(Gterms,2);
I_0 = intwindow(1); I_f = intwindow(2);
if 1 % Sample from prior as proposal distribution!
    Moprop = unifrnd(amo,bmo);
    gWprop = NaN*ones(12,Ny);
    gWprop(M<Mt) = 0;
    gWprop(M>Moprop) = 1;
    gWprop(M<Moprop&M>Mt) = (M(M<Moprop&M>Mt)-Mt)/(Moprop-Mt);
    gprop = diag(gE)*min(gWprop,gT);
    gcurr = Gterms;
    %
    %%%%%%%%%% account for variable integration window:
    if I_0<0; % if we include part of the previous year in each year's modeled growth:
        startmo = 13+I_0;
        endmo = I_f;
        prevseas = [mean(gprop(startmo:12,:),2) gprop(startmo:12,1:end-1)];
        gprop = gprop(1:endmo,:);
        gprop = [prevseas; gprop];
        prevseas = [mean(gcurr(startmo:12,:),2) gcurr(startmo:12,1:end-1)];
        gcurr = gcurr(1:endmo,:);
        gcurr = [prevseas; gcurr];
    else % no inclusion of last year's growth conditions in estimates of this year's growth:
        startmo = I_0+1;
        endmo = I_f;
        gprop = gprop(startmo:endmo,:);
        gcurr = gcurr(startmo:endmo,:);
    end
    %%%%%%%%%%%%
    %
    if length(errorpars) == 1 % White noise error model:
        sigma2rw = errorpars;
        expcurr = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr))).^2);
        expprop = sum((RW(cyrs)'-sqrt(1-sigma2rw)*(sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop))).^2);
        HR = exp(-.5*(expprop-expcurr)/sigma2rw);
    elseif length(errorpars) == 2 % AR(1) error model:
        phi1 = errorpars(1); tau2 = errorpars(2);
        sigma2rw = tau2/(1-phi1^2);
        %
        [iSig] = makeAR1covmat(phi1,tau2,length(cyrs));
        %
        Wcurr = ((sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr)))';
        Wprop = ((sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop)))';
        %
        logLprop = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wprop)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wprop);
        logLcurr = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wcurr)'*iSig*(RW(cyrs) - sqrt(1-sigma2rw)*Wcurr);
        HR = exp(logLprop-logLcurr);
    end
end
% accept or reject the proposal.
if binornd(1,min(HR,1))==1
    Mo = Moprop;
else
    Mo = Mocurr;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Mo] = Mo_lit_aux(Mocurr,M,Mt,gT,RW,errorpars,gE,Gterms,amo,bmo,slp,int,intwindow,cyrs)
% gM/gT is the matrix of gM for all the years and all the months of the previous simulation.
% M is the matrix of soil moisture for all years and all months of the previous simulation.
% amt is the lower bound on the support of the uniform prior distribution for Mt
% bmt is the upper bound on the support of the uniform prior distribution for Mt
%
% SETW 6/10/2010
Ny = size(Gterms,2);
I_0 = intwindow(1); I_f = intwindow(2);
if 1 % Sample from prior as proposal distribution!
    Moprop = slp*betarnd(amo,bmo)+int;
    %     upperlim = normcdf(int+slp,Mocurr,.03);
    %     lowerlim = normcdf(int,Mocurr,.03);
    %     U = unifrnd(lowerlim,upperlim);
    %     Moprop = norminv(U,Mocurr,.03);
    %
    gWprop = NaN*ones(12,Ny);
    gWprop(M<Mt) = 0;
    gWprop(M>Moprop) = 1;
    gWprop(M<Moprop&M>Mt) = (M(M<Moprop&M>Mt)-Mt)/(Moprop-Mt);
    gprop = diag(gE)*min(gWprop,gT);
    gcurr = Gterms;
    %
    %%%%%%%%%% account for variable integration window:
    if I_0<0; % if we include part of the previous year in each year's modeled growth:
        startmo = 13+I_0;
        endmo = I_f;
        prevseas = [mean(gprop(startmo:12,:),2) gprop(startmo:12,1:end-1)];
        gprop = gprop(1:endmo,:);
        gprop = [prevseas; gprop];
        prevseas = [mean(gcurr(startmo:12,:),2) gcurr(startmo:12,1:end-1)];
        gcurr = gcurr(1:endmo,:);
        gcurr = [prevseas; gcurr];
    else % no inclusion of last year's growth conditions in estimates of this year's growth:
        startmo = I_0+1;
        endmo = I_f;
        gprop = gprop(startmo:endmo,:);
        gcurr = gcurr(startmo:endmo,:);
    end
    %%%%%%%%%%%%
    %
    if length(errorpars) == 1 % White noise error model:
        sigma2rw = errorpars;
        expcurr = sum((RW(cyrs)'-sqrt(1-sigma2rw)*...
            (sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr))).^2);
        expprop = sum((RW(cyrs)'-sqrt(1-sigma2rw)*...
            (sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop))).^2);
        HR = exp(-.5*(expprop-expcurr)/sigma2rw);
        %
        %         % must also account for relative prior probs if using "geyer default" sampling:
        %         liklicurr = betapdf((Mocurr-int)/slp,amo,bmo);
        %         likliprop = betapdf((Moprop-int)/slp,amo,bmo);
        %         HR = (likliprop/liklicurr)*exp(-.5*(expprop-expcurr)/sigma2rw);
    elseif length(errorpars) == 2 % AR(1) error model:
        phi1 = errorpars(1); tau2 = errorpars(2);
        sigma2rw = tau2/(1-phi1^2);
        %
        [iSig] = makeAR1covmat(phi1,tau2,length(cyrs));
        %
        Wcurr = ((sum(gcurr(:,cyrs))-mean(sum(gcurr)))/std(sum(gcurr)))';
        Wprop = ((sum(gprop(:,cyrs))-mean(sum(gprop)))/std(sum(gprop)))';
        %
        logLprop = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wprop)'*iSig*(RW(cyrs)-sqrt(1-sigma2rw)*Wprop);
        logLcurr = -.5*(RW(cyrs)-sqrt(1-sigma2rw)*Wcurr)'*iSig*(RW(cyrs)-sqrt(1-sigma2rw)*Wcurr);
        HR = exp(logLprop-logLcurr);
        %
        %         % must also account for relative prior probs if using "geyer default" sampling:
        %         liklicurr = betapdf((Mocurr-int)/slp,amo,bmo);
        %         likliprop = betapdf((Moprop-int)/slp,amo,bmo);
        %         HR = (likliprop/liklicurr)*exp(logLprop-logLcurr);
    end
end
% accept or reject the proposal.
if binornd(1,min(HR,1))==1
    Mo = Moprop;
else
    Mo = Mocurr;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigma2rw,logLdata] = errormodel0_aux(sigma2rwcurr,RW,Gterms,intwindow,cyrs)
% RW = vector of observed annual ring widths
% Gterms = vector of terms that sum together to give the simulated raw ring with index for all
% months (rows) and all years (columns)
% SETW 4/5/2013
%
%%%%%%%%%% account for variable integration window:
I_0 = intwindow(1); I_f = intwindow(2);
if I_0<0; % if we include part of the previous year in each year's modeled growth:
    startmo = 13+I_0;
    endmo = I_f;
    prevseas = [mean(Gterms(startmo:12,:),2) Gterms(startmo:12,1:end-1)];
    Gterms = Gterms(1:endmo,:);
    Gterms = [prevseas; Gterms];
else % no inclusion of last year's growth conditions in estimates of this year's growth:
    startmo = I_0+1;
    endmo = I_f;
    Gterms = Gterms(startmo:endmo,:);
end
%%%%%%%%%%%%
% % sample proposal from the prior:
% sigma2rwprop = unifrnd(0,1); % restricted to the unit interval since sigma^2_rw = (1)/(1+SNR^2) in this model

% try in terms of a uniform prior on the std dev. of the model error,
% rather than a uniform prior on the variance....
sigma2rwprop = (unifrnd(0,1))^2;

%
% accept or reject?
Nyrs = length(cyrs);
Gamma = squeeze(sum(Gterms));
Gammabar = mean(Gamma);
siggamma = std(Gamma);
%
logprop = -.5*sum((RW(cyrs)-sqrt(1-sigma2rwprop)*(Gamma(cyrs)-Gammabar)/siggamma).^2)/sigma2rwprop;
logcurr = -.5*sum((RW(cyrs)-sqrt(1-sigma2rwcurr)*(Gamma(cyrs)-Gammabar)/siggamma).^2)/sigma2rwcurr;
HR = ((sigma2rwcurr/sigma2rwprop)^(Nyrs/2))*exp(logprop-logcurr);
if binornd(1,min(HR,1))==1
    sigma2rw = sigma2rwprop;
    logLdata = logprop-Nyrs/2*log(sigma2rwprop);
else
    sigma2rw = sigma2rwcurr;
    logLdata = logcurr-Nyrs/2*log(sigma2rwcurr);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pars,logLdata] = errormodel1_aux(currpars,RW,Gterms,intwindow,cyrs)
% RW = vector of observed annual ring widths
% Gterms = vector of terms that sum together to give the simulated raw ring with index for all
% months (rows) and all years (columns)
% SETW 4/5/2013
%
%%%%%%%%%% account for variable integration window:
I_0 = intwindow(1); I_f = intwindow(2);
if I_0<0; % if we include part of the previous year in each year's modeled growth:
    startmo = 13+I_0;
    endmo = I_f;
    prevseas = [mean(Gterms(startmo:12,:),2) Gterms(startmo:12,1:end-1)];
    Gterms = Gterms(1:endmo,:);
    Gterms = [prevseas; Gterms];
else % no inclusion of last year's growth conditions in estimates of this year's growth:
    startmo = I_0+1;
    endmo = I_f;
    Gterms = Gterms(startmo:endmo,:);
end
%%%%%%%%%%%%
% read current values of parameters:
phi1curr = currpars(1);
tau2curr = currpars(2);
% if 0 % sample proposal from the prior:
phi1prop = unifrnd(0,1);
tau2prop = unifrnd(0,1);
while tau2prop > 1-phi1prop^2
    % satisfy conditions for stationarity, causality, and also
    % sigma2_w = tau2/(1-phi1^2) <= 1 since sigma2_w = 1/(1+SNR^2) in this model
    phi1prop = unifrnd(0,1);
    tau2prop = unifrnd(0,1);
end

% try in terms of a uniform prior on the std dev. of the model error,
% rather than a uniform prior on the variance....
% sig_w = unifrnd(0,1);
% phi1prop = unifrnd(0,1);
% tau2prop = (1-phi1prop^2)*sig_w^2;



% else % sample proposal using random walk step from current location:
%     phi1prop = phi1curr + .25*randn;
%     tau2prop = tau2curr + .1*randn;
%     while tau2prop > 1-phi1prop^2 || tau2prop <0 || phi1prop < 0
%         % satisfy conditions for stationarity, causality, and also
%         % sigma2_w = tau2/(1-phi1^2) <= 1 since sigma2_w = 1/(1+SNR^2) in this model
%         phi1prop = phi1curr + .25*randn;
%         tau2prop = tau2curr + .1*randn;
%     end
% end
%
% accept or reject?
Ny = length(cyrs);
Gamma = squeeze(sum(Gterms));
Gammabar = mean(Gamma);
siggamma = std(Gamma);
% VS-lite estimate of TRW at current parameter values:
What = ((Gamma(cyrs)-Gammabar)/siggamma)';
%
[iSigprop,detSigprop] = makeAR1covmat(phi1prop,tau2prop,Ny);
[iSigcurr,detSigcurr] = makeAR1covmat(phi1curr,tau2curr,Ny);
alphaprop = sqrt(1-tau2prop/(1-phi1prop^2));
alphacurr = sqrt(1-tau2curr/(1-phi1curr^2));
%
logLprop = -.5*(RW(cyrs)'-alphaprop*What)'*iSigprop*(RW(cyrs)'-alphaprop*What);
logLcurr = -.5*(RW(cyrs)'-alphacurr*What)'*iSigcurr*(RW(cyrs)'-alphacurr*What);
HR = sqrt(detSigcurr/detSigprop)*exp(logLprop-logLcurr);

if binornd(1,min(HR,1))==1
    pars(1) = phi1prop;
    pars(2) = tau2prop;
    logLdata = logLprop-log(detSigprop);
else
    pars(1) = phi1curr;
    pars(2) = tau2curr;
    logLdata = logLcurr-log(detSigcurr);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [invSigma,detSigma] = makeAR1covmat(phi1,tau2,N)
%%% [Sigma,invSigma,detSigma] = makeAR1covmat(phi1,tau2,N)
% Make approximate joint covariance matrix, inverse covariace matrix, and covariance matrix
% determinant for N sequential observations that follow the AR(1) model
% X_t = phi1 X_t-1 + eps_t, where eps_t ~ N(0,tau^2)

A = -phi1*eye(N);
superdiag = sub2ind([N N],(1:(N-1))',(2:N)');
A(superdiag) = 1;
% Now Var(A* e) \approx tau2*I, so
% Sigma \approx tau2* inv(A)*inv(A')
% and
% invSigma \approx (1/tau2)* A'*A

%Sigma = tau2*(A\eye(N))*(A'\eye(N));
invSigma = (1/tau2)*(A')*A;
detSigma = (tau2/phi1^2)^N;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Rhat] = gelmanrubin92(Nsamp,Nbi,varargin)
% Usage: Rhat = gelmanrubin92(Nsamp,Nbi,chain1,chain2,...,chainN)
% Nsamp = number of iterations in each sample
% Nbi = number to consider "burn-in".
% chain1, ..., chainN must have dimension Nsamp x 1.
% SETW 1/26/2011

% Number of chains:
m = nargin-2;
for i = 1:m
    eval(['chain' num2str(i) ' = varargin{i};']);
end
% number of non-burn-in iterations:
n = Nsamp-Nbi;
%
Xbar = NaN(1,m);
Xs2 = NaN(1,m);
allX = NaN(1,n*m);
for i = 1:m
    eval(['X = chain' num2str(i) ';']);
    % within-chain means of X:
    Xbar(i) = nanmean(X(Nbi+1:length(X)));
    % within-chain variances of X:
    Xs2(i) = nanvar(X(Nbi+1:length(X)));
    allX((i-1)*n+1:i*n) = X(Nbi+1:Nsamp);
end
% mean across chains of mean X in each month:
Xbarbar = mean(Xbar);
%
BX = n*(sum(((Xbar-repmat(Xbarbar,1,m)).^2),2))/(m-1);
%
WX = nanmean(Xs2,2);
%
muhatX = nanmean(allX,2);
%
sig2hatX = (n-1)*WX/n + BX/n; % G&R92 eqn. 3
%
VhatX = sig2hatX + BX/(m*n);
varhatXs2 = var(Xs2,0,2);
%
covhatXs2Xbar2 = sum((Xs2-nanmean(Xs2)).*(Xbar.^2-nanmean(Xs2.^2)))/m; 
covhatXs2Xbar = sum((Xs2-nanmean(Xs2)).*(Xbar-nanmean(Xs2)))/m;
%
covhatXs2Xbar2 = covhatXs2Xbar2';covhatXs2Xbar = covhatXs2Xbar';
%
varhatVhatX = (((n-1)/n)^2)*varhatXs2/m + (((m+1)/(m*n))^2)*2*BX.^2/(m-1) +...
    2*((m+1)*(n-1)/(m*n^2))*n*(covhatXs2Xbar2-2*Xbarbar.*covhatXs2Xbar)/m;
dfX = 2*(VhatX.^2)./varhatVhatX;
%
Rhat = (VhatX./WX).*(dfX./(dfX-2));
end

