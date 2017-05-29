function [oriabu,bioabu,oriiso,bioiso] = turbo2(abu,iso,mxl,numb)
% The MATLAB program TURBO2 can be used to simulate the effects of
% bioturbation on single sediments particles.
%
% ABU = series of abundances of carrier type 1 down core
% ISO = series of isotopes measured on carrier
% MXL = series of mixed layer thicknesses down core
% NUMB = number of carriers to be measured
% ORIABU = original abundances of both carrier types 1 and 2
% BIOABU = bioturbated abundances of both carriers types 1 and 2
% ORIISO = original isotope signature of both carrier  types 1 and 2
% BIOISO = bioturbated isotope signature of both carrier types 1 and 2
% Martin Trauth 18 July 2012

% CITATION:
% "TURBO2: A MATLAB simulation to study the effects of bioturbation on paleoceanographic time series" 
%  by Martin H. Trauth, University of Potsdam, published in Computers & Geosciences, 2013

% FULL GITHUB EXAMPLES: https://github.com/cageo/Trauth-2013

nrows = max(mxl)+0;
ncols = max(abu)+50;
sedabu = NaN(nrows,ncols);
sediso = NaN(nrows,ncols);
%h = waitbar(0,'Mixing Process ...');
for i = 1:length(abu)
  %waitbar(i/length(abu))
  rncols = randperm(ncols);
  sedabu(size(sedabu,1)+1,1:ncols) = 2*ones(1,ncols);
  sedabu(size(sedabu,1),1:abu(i)) = ones(1,abu(i));
  sedabu(size(sedabu,1),:) = sedabu(size(sedabu,1),rncols);
  sediso(size(sediso,1)+1,1:ncols) = iso(i)*ones(1,ncols);
  for j = 1 : ncols
      z = randperm(mxl(i));
      ns(1:mxl(i),j) = sedabu(size(sedabu,1)-mxl(i)+1:size(sedabu,1),j);
      sedabu(size(sedabu,1)-mxl(i)+1:size(sedabu,1),j) = ns(z,j);
      ni(1:mxl(i),j) = sediso(size(sediso,1)-mxl(i)+1:size(sediso,1),j);
      sediso(size(sediso,1)-mxl(i)+1:size(sediso,1),j) = ni(z,j);
  end
end
%close(h)
%
%flipud(sedabu(:,1:5))
%flipud(sediso(:,1:5))
%
clear ns ni i j mxl z
sedabu = sedabu(nrows+1:end,:);
sediso = sediso(nrows+1:end,:);
oriabu(:,1) = abu;
oriabu(:,2) = ncols-abu;
bioabu(:,1) = sum(sedabu==1,2);
bioabu(:,2) = sum(sedabu==2,2);
oriiso(:,1) = iso;
oriiso(:,2) = iso;
bioiso1 = sediso;
bioiso2 = sediso;
bioiso1(sedabu~=1) = NaN;
bioiso2(sedabu~=2) = NaN;
biopart1 = NaN(size(bioiso1));
biopart2 = NaN(size(bioiso2));

for i = 1:2
    biopart1(i,1:bioabu(i,1)) = bioiso1(i,isnan(bioiso1(i,:))==0);
    biopart2(i,1:bioabu(i,2)) = bioiso2(i,isnan(bioiso2(i,:))==0);
end


biopart1 = biopart1(:,1:numb);
biopart2 = biopart2(:,1:numb);

for i = 1:length(abu)
    bioiso(i,1) = nanmean(biopart1(i,:));
    bioiso(i,2) = nanmean(biopart2(i,:));
end
% oriabu = flipud(oriabu);
% oriiso = flipud(oriiso);
% bioabu = flipud(bioabu);
% bioiso = flipud(bioiso);
% %% example MATLAB script to run TURBO2
% clear
% data = load('turbo2input_termination.txt');
% %%
% age   = data(:,1);
% mxl   = data(:,2);
% abu   = data(:,3);
% iso   = data(:,4);
% lngth = length(data(:,1));
% numb = 50;
% %%
% [oriabu,bioabu,oriiso,bioiso] = turbo2(abu,iso,mxl,numb);
% 
% %%
% mxltext = num2str(mean(mxl));
% numbtxt = num2str(numb,2);
% %%
% subplot(2,2,1)
% plot(1:lngth,oriabu(:,1),'k'), hold on
% plot(1:lngth,bioabu(:,1),'b')
% plot(1:lngth,numb*ones(lngth),'g')
% set(gca,'XGrid','On','YGrid','On')
% title('Abundance of Carrier 1')
% subplot(2,2,2)
% plot(1:lngth,oriabu(:,2),'k'), hold on
% plot(1:lngth,bioabu(:,2),'r')
% plot(1:lngth,numb*ones(lngth),'g')
% set(gca,'XGrid','On','YGrid','On')
% title('Abundance of Carrier 2')
% subplot(2,2,3)
% plot(1:lngth,oriiso(:,1),'k'), hold on
% plot(1:lngth,bioiso(:,1),'b')
% set(gca,'YDir','Reverse','XGrid','On','YGrid','On')
% title('Isotopes of Carrier 1')
% subplot(2,2,4)
% plot(1:lngth,oriiso(:,2),'k'), hold on
% plot(1:lngth,bioiso(:,2),'r')
% set(gca,'YDir','Reverse','XGrid','On','YGrid','On')
% title('Isotopes of Carrier 2')
% %%
% printfilename = ['turbo2_fig1_',mxltext,'cm_',numbtxt,'carriers.tiff'];
% print('-dtiff',printfilename)
% %%
% figure, hold on
% plot(1:lngth,oriiso(:,1),'k')
% plot(1:lngth,bioiso(:,1),'b')
% plot(1:lngth,bioiso(:,2),'r')
% set(gca,'YDir','Reverse','XGrid','On','YGrid','On','Box','On')
% titletxt = ['Isotopes of Carriers 1+2, ',mxltext,...
%     ' cm Mixed Layer, ',numbtxt,' Carriers'];
% title(titletxt)
% legend('Original Isotopes','carriers 1','carriers 2')
% %%
% printfilename = ['turbo2_fig2_',mxltext,'cm_',numbtxt,'carriers.tiff'];
% print('-dtiff',printfilename)
