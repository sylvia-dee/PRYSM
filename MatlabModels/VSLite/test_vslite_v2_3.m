function varargout = test_vslite_v2_3(site)
% Basic usage:  test_vslite(site)
% Input: 'site' is either 1, 2, 3, or 4, runs VS-Lite at one of four test sites.
%   Site 1: 36.45N, -118.22E, 'ca530'
%   Site 2: 34.17N, -117.12E, 'ca544'
%   Site 3: 39.02N, -122.82E, 'ca615'
%   Site 4: 39.32N, -106.08E, 'co523'
% 
% Or can specify outputs: [trw,gT,gM,gE,M] = test_vslite_v2_3
%
% Monthly input climate data is from the OSU PRISM data product for the years 1895 through 1984. 
% Observed chronologies are from the International Tree Ring Data Bank, but have been standardized
% for the period 1895-1984.
% Parameters used are the modes of those found to be 'optimal' over 100 calibration/verification
% tests in Tolwinski-Ward et al., 2010
% SETW 11/2010
% SETW: Updated 11/2011 to run with updated VSLite version 2.2 
% SETW: Updated 7/2013 to run with VS-Lite v2.3, and to additionally 
% estimate parameters of the model at each site using Bayesian scheme 
% detailed in Tolwinski-Ward et al 2013 Climate of the Past paper.

% If user doesn't specify site, just run the first one.
if nargin == 0; site = 1; end

% load the test data:
load vslite_testdata.mat

% select climate and location data for the chosen test site:
T = m08clim(site).T;
P = m08clim(site).P;
phi = sitecoords(site,1);

% estimate the climate response parameters:
disp('Performing Bayesian estimation of VS-Lite parameters for chosen site.')
[T1,T2,M1,M2] = estimate_vslite_params_v2_3(T,P,phi,trw_obs(:,site)','nsamp',2000);

% Run VS-Lite.
[trw,gT,gM,gE,M] = VSLite_v2_3(syear,eyear,phi,T1,T2,M1,M2,T,P);
% Draw some output.
figure;
set(gcf,'units','normalized','position',[.25 .25 .5 .4])

if exist('worldmap')
    subplot(2,2,2);
    plot(mean(gM,2),'b--'); xlim([1 12]); hold on;
    plot(mean(gT,2),'r');
    title('Mean gT (red) and gM (blue)')
    subplot(2,2,1)
    worldmap([27 50],[-127 -65])
    load coast
    plotm(lat,long,'k')
    hold on; plotm(sitecoords(site,1),sitecoords(site,2),'*r','markersize',8)
    title('Site location')
    subplot(2,1,2);
    plot(syear:eyear,trw,'rx-'); hold on
    plot(syear:eyear,trw_obs(:,site),'kd-'); xlim([syear eyear]);
    eval(['title(''Simulated RW (red) and observed (black)'')']);
    ylim([min(min(trw),min(trw_obs(:,site))) max(max(trw),max(trw_obs(:,site)))]);
elseif exist('m_proj')
    subplot(2,2,2);
    plot(mean(gM,2),'b--'); xlim([1 12]); hold on;
    plot(mean(gT,2),'r');
    title('Mean gT (red) and gM (blue)')
    subplot(2,2,1)
    m_proj('Equidistant Cylindrical','longitude',[-127 -65],'latitude',[27 50])
    m_coast('color','k')
    m_grid
    m_line(sitecoords(site,2),sitecoords(site,1),'marker','*','color','r','markersize',8)
    title('Site location')
    subplot(2,1,2);
    plot(syear:eyear,trw,'rx-'); hold on
    plot(syear:eyear,trw_obs(:,site),'kd-'); xlim([syear eyear]);
    eval(['title(''Simulated RW (red) and observed (black)'')']);
    ylim([min(min(trw),min(trw_obs(:,site))) max(max(trw),max(trw_obs(:,site)))]);
else
    subplot(2,1,1);
    plot(mean(gM,2),'b--'); xlim([1 12]); hold on;
    plot(mean(gT,2),'r');
    title('Mean gT (red) and gM (blue)')
    subplot(2,1,2);
    plot(syear:eyear,trw,'rx-'); hold on
    plot(syear:eyear,trw_obs(:,site),'kd-'); xlim([syear eyear]);
    eval(['title(''Simulated RW (red) and observed (black)'')']);
end

if nargout > 0 
   varargout{1} = trw;
   varargout{2} = gT;
   varargout{3} = gM;
   varargout{4} = gE;
   varargout{5} = M;
end
    