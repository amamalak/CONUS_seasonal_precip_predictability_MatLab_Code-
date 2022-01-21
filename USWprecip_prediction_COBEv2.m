clear 
clc
close all 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This is the code script for the study published in Mamalakis et al. (2022). 
% The exact data sources used for the study are free online and urls are provided.
% However, if you want to use different data (or updated versions of it), 
% you may need to modify certain parts of the script accordingly. 
% Please make sure you cite this study in case of using this code or the 
% methodology that is based on. 
% Questions might be addressed to: amamalak@rams.colostate.edu
%
% Editor: Antonios Mamalakis, PhD.
%
% Citation:
% Mamalakis, A., A. AghaKouchak, J.T. Randerson, E. Foufoula-Georgiou
% (2022) Hotspots of predictability: Identifying regions of high
% precipitation predictability at seasonal timescales from limited time
% series observations, Water Resources Research.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loading precipitation data

% precipitation data are downloaded from CPC US at https://psl.noaa.gov/data/gridded/data.unified.daily.conus.html
 
% precip units: mm/d 
% period coverage: 01/1948 - 06/2018 
ncid0= netcdf.open('precipitation/precip.V1.0.mon.mean.nc','NC_NOWRITE'); %change this file directory if necessary 
varid0= netcdf.inqVarID(ncid0,'precip');
monthlyprecip=netcdf.getVar(ncid0,varid0,'double');
% latitude/longitude
varid_lat= netcdf.inqVarID(ncid0,'lat');
lat_prcp=netcdf.getVar(ncid0,varid_lat,'double');
varid_lon= netcdf.inqVarID(ncid0,'lon');
lon_prcp=netcdf.getVar(ncid0,varid_lon,'double');
netcdf.close(ncid0);

monthlyprecip=permute(monthlyprecip,[2,1,3]); %array rearrangement
for i=1:length(lat_prcp)
    for j=1:length(lon_prcp)
        
monthlyprecip(i,j,monthlyprecip(i,j,:)<0)=NaN; %everything negative is a nan (missing) value

    end
end

load coastlines
coastlon(coastlon<0)=coastlon(coastlon<0)+360;

% Plot precipitation map as an example to double-check you read the data well from the website!!
% compare with the plots in the website -->  https://psl.noaa.gov/data/gridded/data.unified.daily.conus.html
figure()
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'}); %do not plot Alaska and Hawaii
contourfm(lat_prcp,lon_prcp,monthlyprecip(:,:,end-6),0:1:12,'LineStyle','none')  % example plot is 12/2017
geoshow(ax, states,'FaceColor','none')
caxis([0 12])
contourcbar

% also compare the min and max values with the ones shown in the website
max(max(monthlyprecip(:,:,end-6)))
min(min(monthlyprecip(:,:,end-6)))

%% Performing spatial upscaling (so that analysis is computationally cheaper)

% consider to disregard this code cell if the upscaling is not warranted in
% your analysis.

lat_prcp_raw=lat_prcp; % original grid at 0.25 by 0.25
lon_prcp_raw=lon_prcp;
clear lat_prcp lon_prcp
lat_prcp=[20:50]'; % new grid; new resolution is 1 by 1 
lon_prcp=[230:305]';
dlat=lat_prcp(2)-lat_prcp(1); % 
dlon=lon_prcp(2)-lon_prcp(1);
monthlyprecip_raw=monthlyprecip;
clear monthlyprecip

for i=1:length(lat_prcp)
    for j=1:length(lon_prcp)
        
        [~, lat1]=min(abs(lat_prcp(i)-dlat/2-lat_prcp_raw));
        [~, lat2]=min(abs(lat_prcp(i)+dlat/2-lat_prcp_raw));
        [~, lon1]=min(abs(lon_prcp(j)-dlon/2-lon_prcp_raw));
        [~, lon2]=min(abs(lon_prcp(j)+dlon/2-lon_prcp_raw));
        % spatial averaging
        monthlyprecip(i,j,:)= nanmean(nanmean(monthlyprecip_raw(lat1:lat2,lon1:lon2,:),1),2);

        clear lat1 lat2 lon1 lon2
    end
end

% Plot again the same example map that was plotted in the previous code
% cell to verify the upscale looks right
figure()
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,monthlyprecip(:,:,end-6),0:1:12,'LineStyle','none')  % 12/2017
geoshow(ax, states,'FaceColor','none')
caxis([0 12])
contourcbar

max(max(monthlyprecip(:,:,end-6)))
min(min(monthlyprecip(:,:,end-6)))

%% Consider only the last 50 years of data 
%(to avoid issues of unreliable earlier obs and changes in teleconnections)

% keep only 01/1967 - 12/2017
% Here you may need to modify. The version of the data I am using is during 01/1948 - 06/2018.
% If your data cover a different period, modify accordingly.
monthlyprecip(:,:,end-5:end)=[];
temp=monthlyprecip;
clear monthlyprecip
monthlyprecip=temp(:,:,19*12+1:end);
clear temp

%% Study's focus: Nov-Mar season (modify the next lines accordingly if focusing on a different season)

for t=1:size(monthlyprecip,3)/12-1
    NMprecip(:,:,t)=sum(monthlyprecip(:,:,(t-1)*12+11:(t-1)*12+15),3)*30; %Nov-Mar (units: mm)
    MSprecip(:,:,t)=sum(monthlyprecip(:,:,(t-1)*12+17:(t-1)*12+21),3)*30; %May-Sep (units: mm) this is not used in the study; consider to disregard.
    
    ANNprecip(:,:,t)=sum(monthlyprecip(:,:,(t-1)*12+10:(t-1)*12+21),3)*30; %ANN (units: mm)
end

%%%%%%%%%%%%% Figure 2a %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,nanmean(NMprecip./ANNprecip,3)*100,0:10:100,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
caxis([0 100])
contourcbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure() % not used in the study %
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,nanmean(MSprecip./ANNprecip,3)*100,0:10:100,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
caxis([0 100])
contourcbar

%% Other precipitation statistics 

acfUS=nan*NMprecip(:,:,1); 

for i=1:length(lat_prcp)
    for j=1:length(lon_prcp)
        
        temp=squeeze(NMprecip(i,j,:));
        NMprecip2(i,j,:)= detrend(temp)+nanmean(temp); %linear detrending
       
        
        if sum(isnan(temp))==0
        acf=autocorr(squeeze(NMprecip2(i,j,:)));
        acfUS(i,j)=acf(2); % lag-1 autocorrelation in the precip series 
        clear acf
        end
        clear temp
        
    end
end

%%%%%%%%%%%%% Figure 2b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,nanmean(NMprecip2,3),0:100:1200,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
caxis([0 1200])
contourcbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Figure 2c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,nanstd(NMprecip2,0,3)./nanmean(NMprecip2,3),0:0.1:1,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
caxis([0 1])
contourcbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Figure 2d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[LON_prcp,LAT_prcp]=meshgrid(lon_prcp,lat_prcp);
figure()
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,acfUS,-1:0.2:1,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
hold on 
plotm(reshape(LAT_prcp(acfUS>2/sqrt(size(NMprecip2,3))),[],1),reshape(LON_prcp(acfUS>2/sqrt(size(NMprecip2,3))),[],1),'k.')    
caxis([-1 1])
contourcbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Loading sst data

% SST data is downloaded from COBESST2 https://www.esrl.noaa.gov/psd/data/gridded/data.cobe2.html

close all

% SST units: Celsius
% period coverage: 01/1850 - 12/2018
ncid0= netcdf.open('sst_data/sst.mon.mean_cobe_v2.nc','NC_NOWRITE'); %change this file directory if necessary 
varid0= netcdf.inqVarID(ncid0,'sst');
monthlySST=netcdf.getVar(ncid0,varid0,'double');

varid_lat= netcdf.inqVarID(ncid0,'lat');
lat_sst=netcdf.getVar(ncid0,varid_lat,'double');
varid_lon= netcdf.inqVarID(ncid0,'lon');
lon_sst=netcdf.getVar(ncid0,varid_lon,'double');
[LON_sst,LAT_sst]=meshgrid(lon_sst,lat_sst);
netcdf.close(ncid0);

monthlySST=permute(monthlySST,[2,1,3]); %array rearrangement
for i=1:length(lat_sst)
    for j=1:length(lon_sst)
        
monthlySST(i,j,abs(monthlySST(i,j,:))>500)=NaN; %everything above/below +/-500 Celsius is a nan (missing) value

    end
end

%
load coastlines
coastlon(coastlon<0)=coastlon(coastlon<0)+360;

% Plot SST map as an example to double-check you read the data well from the website!!
% compare with the plots in the website -->  COBESST2 https://www.esrl.noaa.gov/psd/data/gridded/data.cobe2.html

figure()
contourf(LON_sst,LAT_sst,monthlySST(:,:,end))% example is 12/2018
hold on 
plot(coastlon,coastlat,'.')

% also compare the min and max values with the ones shown in the website
max(max(monthlySST(:,:,end)))
min(min(monthlySST(:,:,end)))

% keep only 01/1967- 12/2016 (study's focus)
% you may need to modify these lines depending on what is your period of
% focus or what the original coverage of your data is!
monthlySST(:,:,end-23:end)=[];
temp=monthlySST;
clear monthlySST
monthlySST=temp(:,:,117*12+1:end);
clear temp
%% definition of SST indices

% NZI 
% first determine index boundaries
[~, lat1]=min(abs(lat_sst+25)) ;
[~, lat2]=min(abs(lat_sst+40)) ;
[~, lon1]=min(abs(lon_sst-170)) ;
[~, lon2]=min(abs(lon_sst-200)) ;
NZI=squeeze(nanmean(nanmean(monthlySST(lat1:lat2,lon1:lon2,:),1),2));% average without considering area differences
for i=1:size(monthlySST,3)
% in the latitude direction the grids become smaller as we move poleward
COSINE=cosd(repmat(lat_sst(lat1:lat2),1,lon2-lon1+1));
sstt1=monthlySST(lat1:lat2,lon1:lon2,i).*COSINE;
% area-weighted average based on the cosine function
NZI_w(i,1)=sum(sum(sstt1(isnan(sstt1)==0)))/sum(sum(COSINE(isnan(sstt1)==0)));
end
clear COSINE sstt1
% plot the index series based on the two ways and compare
figure()
plot(NZI)
hold on
plot(NZI_w,'r')
title(['r = ',num2str(corr(NZI,NZI_w))])
max(abs(NZI-NZI_w))
% for all indices, we should be using the index that was calculated considering the area
% differences (although most of the times there is no big difference
% between the two)


% EPI 
[~, lat1]=min(abs(lat_sst-20)) ;
[~, lat2]=min(abs(lat_sst-5)) ;
[~, lon1]=min(abs(lon_sst-130)) ;
[~, lon2]=min(abs(lon_sst-160)) ;
EPI=squeeze(nanmean(nanmean(monthlySST(lat1:lat2,lon1:lon2,:),1),2));% without considering area differences
for i=1:size(monthlySST,3)
%area-weighted
COSINE=cosd(repmat(lat_sst(lat1:lat2),1,lon2-lon1+1));
sstt1=monthlySST(lat1:lat2,lon1:lon2,i).*COSINE;
EPI_w(i,1)=sum(sum(sstt1(isnan(sstt1)==0)))/sum(sum(COSINE(isnan(sstt1)==0)));
end
clear COSINE sstt1
figure()
plot(EPI)
hold on
plot(EPI_w,'r')
title(['r = ',num2str(corr(EPI,EPI_w))])
max(abs(EPI-EPI_w))


% Nino3.4
[~, lat1]=min(abs(lat_sst-5)) ;
[~, lat2]=min(abs(lat_sst+5)) ;
[~, lon1]=min(abs(lon_sst-190)) ;
[~, lon2]=min(abs(lon_sst-240)) ;
Nino34=squeeze(nanmean(nanmean(monthlySST(lat1:lat2,lon1:lon2,:),1),2));% without considering area differences
for i=1:size(monthlySST,3)
%area-weighted
COSINE=cosd(repmat(lat_sst(lat1:lat2),1,lon2-lon1+1));
sstt1=monthlySST(lat1:lat2,lon1:lon2,i).*COSINE;
Nino34_w(i,1)=sum(sum(sstt1(isnan(sstt1)==0)))/sum(sum(COSINE(isnan(sstt1)==0)));
end
clear COSINE sstt1
figure()
plot(Nino34)
hold on
plot(Nino34_w,'r')
title(['r = ',num2str(corr(Nino34,Nino34_w))])
max(abs(Nino34-Nino34_w))

% Nino4
[~, lat1]=min(abs(lat_sst-5)) ;
[~, lat2]=min(abs(lat_sst+5)) ;
[~, lon1]=min(abs(lon_sst-160)) ;
[~, lon2]=min(abs(lon_sst-210)) ;
Nino4=squeeze(nanmean(nanmean(monthlySST(lat1:lat2,lon1:lon2,:),1),2));% without considering area differences
for i=1:size(monthlySST,3)
%area-weighted
COSINE=cosd(repmat(lat_sst(lat1:lat2),1,lon2-lon1+1));
sstt1=monthlySST(lat1:lat2,lon1:lon2,i).*COSINE;
Nino4_w(i,1)=sum(sum(sstt1(isnan(sstt1)==0)))/sum(sum(COSINE(isnan(sstt1)==0)));
end
clear COSINE sstt1
figure()
plot(Nino4)
hold on
plot(Nino4_w,'r')
title(['r = ',num2str(corr(Nino4,Nino4_w))])
max(abs(Nino4-Nino4_w))

% Nino3
[~, lat1]=min(abs(lat_sst-5)) ;
[~, lat2]=min(abs(lat_sst+5)) ;
[~, lon1]=min(abs(lon_sst-210)) ;
[~, lon2]=min(abs(lon_sst-270)) ;
Nino3=squeeze(nanmean(nanmean(monthlySST(lat1:lat2,lon1:lon2,:),1),2));% without considering area differences
for i=1:size(monthlySST,3)
%area-weighted
COSINE=cosd(repmat(lat_sst(lat1:lat2),1,lon2-lon1+1));
sstt1=monthlySST(lat1:lat2,lon1:lon2,i).*COSINE;
Nino3_w(i,1)=sum(sum(sstt1(isnan(sstt1)==0)))/sum(sum(COSINE(isnan(sstt1)==0)));
end
clear COSINE sstt1
figure()
plot(Nino3)
hold on
plot(Nino3_w,'r')
title(['r = ',num2str(corr(Nino3,Nino3_w))])
max(abs(Nino3-Nino3_w))

% Nino12
[~, lat1]=min(abs(lat_sst-0)) ;
[~, lat2]=min(abs(lat_sst+10)) ;
[~, lon1]=min(abs(lon_sst-270)) ;
[~, lon2]=min(abs(lon_sst-280)) ;
Nino12=squeeze(nanmean(nanmean(monthlySST(lat1:lat2,lon1:lon2,:),1),2));% without considering area differences
for i=1:size(monthlySST,3)
%area-weighted
COSINE=cosd(repmat(lat_sst(lat1:lat2),1,lon2-lon1+1));
sstt1=monthlySST(lat1:lat2,lon1:lon2,i).*COSINE;
Nino12_w(i,1)=sum(sum(sstt1(isnan(sstt1)==0)))/sum(sum(COSINE(isnan(sstt1)==0)));
end
clear COSINE sstt1
figure()
plot(Nino12)
hold on
plot(Nino12_w,'r')
title(['r = ',num2str(corr(Nino12,Nino12_w))])
max(abs(Nino12-Nino12_w))

% TNI
TNI_w=Nino12_w-Nino4_w;

% TNA
[~, lat1]=min(abs(lat_sst-25)) ;
[~, lat2]=min(abs(lat_sst-5)) ;
[~, lon1]=min(abs(lon_sst-305)) ;
[~, lon2]=min(abs(lon_sst-345)) ;
TNA=squeeze(nanmean(nanmean(monthlySST(lat1:lat2,lon1:lon2,:),1),2));% without considering area differences
for i=1:size(monthlySST,3)
%area-weighted
COSINE=cosd(repmat(lat_sst(lat1:lat2),1,lon2-lon1+1));
sstt1=monthlySST(lat1:lat2,lon1:lon2,i).*COSINE;
TNA_w(i,1)=sum(sum(sstt1(isnan(sstt1)==0)))/sum(sum(COSINE(isnan(sstt1)==0)));
end
clear COSINE sstt1
figure()
plot(TNA)
hold on
plot(TNA_w,'r')
title(['r = ',num2str(corr(TNA,TNA_w))])
max(abs(TNA-TNA_w))

% TSA
[~, lat1]=min(abs(lat_sst-0)) ;
[~, lat2]=min(abs(lat_sst+20)) ;
[~, lon1]=min(abs(lon_sst-10)) ;
[~, lon2]=min(abs(lon_sst-330)) ;
TSA=squeeze(nanmean(nanmean(monthlySST(lat1:lat2,[1:lon1,lon2:length(lon_sst)],:),1),2));% without considering area differences
for i=1:size(monthlySST,3)
%area-weighted
COSINE=cosd(repmat(lat_sst(lat1:lat2),1,lon1-1+length(lon_sst)-lon2+2));
sstt1=monthlySST(lat1:lat2,[1:lon1,lon2:end],i).*COSINE;
TSA_w(i,1)=sum(sum(sstt1(isnan(sstt1)==0)))/sum(sum(COSINE(isnan(sstt1)==0)));
end
clear COSINE sstt1
figure()
plot(TSA)
hold on
plot(TSA_w,'r')
title(['r = ',num2str(corr(TSA,TSA_w))])
max(abs(TSA-TSA_w))

% NAT
[~, lat1]=min(abs(lat_sst-20)) ;
[~, lat2]=min(abs(lat_sst-5)) ;
[~, lon1]=min(abs(lon_sst-320)) ;
[~, lon2]=min(abs(lon_sst-340)) ;
NAT=squeeze(nanmean(nanmean(monthlySST(lat1:lat2,lon1:lon2,:),1),2));% without considering area differences
for i=1:size(monthlySST,3)
%area-weighted
COSINE=cosd(repmat(lat_sst(lat1:lat2),1,lon2-lon1+1));
sstt1=monthlySST(lat1:lat2,lon1:lon2,i).*COSINE;
NAT_w(i,1)=sum(sum(sstt1(isnan(sstt1)==0)))/sum(sum(COSINE(isnan(sstt1)==0)));
end
clear COSINE sstt1
figure()
plot(NAT)
hold on
plot(NAT_w,'r')
title(['r = ',num2str(corr(NAT,NAT_w))])
max(abs(NAT-NAT_w))

% SAT
[~, lat1]=min(abs(lat_sst+5)) ;
[~, lat2]=min(abs(lat_sst+20)) ;
[~, lon1]=min(abs(lon_sst-5)) ;
[~, lon2]=min(abs(lon_sst-345)) ;
SAT=squeeze(nanmean(nanmean(monthlySST(lat1:lat2,[1:lon1,lon2:length(lon_sst)],:),1),2));% without considering area differences
for i=1:size(monthlySST,3)
%area-weighted
COSINE=cosd(repmat(lat_sst(lat1:lat2),1,lon1-1+length(lon_sst)-lon2+2));
sstt1=monthlySST(lat1:lat2,[1:lon1,lon2:end],i).*COSINE;
SAT_w(i,1)=sum(sum(sstt1(isnan(sstt1)==0)))/sum(sum(COSINE(isnan(sstt1)==0)));
end
clear COSINE sstt1
figure()
plot(SAT)
hold on
plot(SAT_w,'r')
title(['r = ',num2str(corr(SAT,SAT_w))])
max(abs(SAT-SAT_w))

% TASI
TASI_w=NAT_w-SAT_w;

% AMO
[~, lat1]=min(abs(lat_sst-70)) ;
[~, lat2]=min(abs(lat_sst-0)) ;
[~, lon1]=min(abs(lon_sst-280)) ;
[~, lon2]=min(abs(lon_sst-360)) ;
AMO=squeeze(nanmean(nanmean(monthlySST(lat1:lat2,lon1:lon2,:),1),2));% without considering area differences
for i=1:size(monthlySST,3)
%area-weighted
COSINE=cosd(repmat(lat_sst(lat1:lat2),1,lon2-lon1+1));
sstt1=monthlySST(lat1:lat2,lon1:lon2,i).*COSINE;
AMO_w(i,1)=sum(sum(sstt1(isnan(sstt1)==0)))/sum(sum(COSINE(isnan(sstt1)==0)));
end
clear COSINE sstt1
figure()
plot(AMO)
hold on
plot(AMO_w,'r')
title(['r = ',num2str(corr(AMO,AMO_w))])
max(abs(AMO-AMO_w))


% WTIO
[~, lat1]=min(abs(lat_sst-10)) ;
[~, lat2]=min(abs(lat_sst+10)) ;
[~, lon1]=min(abs(lon_sst-50)) ;
[~, lon2]=min(abs(lon_sst-70)) ;
WTIO=squeeze(nanmean(nanmean(monthlySST(lat1:lat2,lon1:lon2,:),1),2));% without considering area differences
for i=1:size(monthlySST,3)
%area-weighted
COSINE=cosd(repmat(lat_sst(lat1:lat2),1,lon2-lon1+1));
sstt1=monthlySST(lat1:lat2,lon1:lon2,i).*COSINE;
WTIO_w(i,1)=sum(sum(sstt1(isnan(sstt1)==0)))/sum(sum(COSINE(isnan(sstt1)==0)));
end
clear COSINE sstt1
figure()
plot(WTIO)
hold on
plot(WTIO_w,'r')
title(['r = ',num2str(corr(WTIO,WTIO_w))])
max(abs(WTIO-WTIO_w))

% SETIO
[~, lat1]=min(abs(lat_sst-0)) ;
[~, lat2]=min(abs(lat_sst+10)) ;
[~, lon1]=min(abs(lon_sst-90)) ;
[~, lon2]=min(abs(lon_sst-110)) ;
SETIO=squeeze(nanmean(nanmean(monthlySST(lat1:lat2,lon1:lon2,:),1),2));% without considering area differences
for i=1:size(monthlySST,3)
%area-weighted
COSINE=cosd(repmat(lat_sst(lat1:lat2),1,lon2-lon1+1));
sstt1=monthlySST(lat1:lat2,lon1:lon2,i).*COSINE;
SETIO_w(i,1)=sum(sum(sstt1(isnan(sstt1)==0)))/sum(sum(COSINE(isnan(sstt1)==0)));
end
clear COSINE sstt1
figure()
plot(SETIO)
hold on
plot(SETIO_w,'r')
title(['r = ',num2str(corr(SETIO,SETIO_w))])
max(abs(SETIO-SETIO_w))

% DMI
DMI_w=WTIO_w-SETIO_w;

% SWIO
[~, lat1]=min(abs(lat_sst+25)) ;
[~, lat2]=min(abs(lat_sst+32)) ;
[~, lon1]=min(abs(lon_sst-31)) ;
[~, lon2]=min(abs(lon_sst-45)) ;
SWIO=squeeze(nanmean(nanmean(monthlySST(lat1:lat2,lon1:lon2,:),1),2));% without considering area differences
for i=1:size(monthlySST,3)
%area-weighted
COSINE=cosd(repmat(lat_sst(lat1:lat2),1,lon2-lon1+1));
sstt1=monthlySST(lat1:lat2,lon1:lon2,i).*COSINE;
SWIO_w(i,1)=sum(sum(sstt1(isnan(sstt1)==0)))/sum(sum(COSINE(isnan(sstt1)==0)));
end
clear COSINE sstt1
figure()
plot(SWIO)
hold on
plot(SWIO_w,'r')
title(['r = ',num2str(corr(SWIO,SWIO_w))])
max(abs(SWIO-SWIO_w))


%% definition of predictors
% Note 1: we focus on the months Aug-Sep (WE WANT NO OVERLAP WITH THE SEASON WE ARE PREDICTING PRECIP FOR) 
% Modify accordingly if necessary.
% Note 2: we are using the indices that were calculated considering the
% area differences of the grids

for t=1:length(NZI_w)/12

Nino34pr(t,1)=nanmean(Nino34_w((t-1)*12+8:(t-1)*12+10)); %Aug-Oct
Nino4pr(t,1)=nanmean(Nino4_w((t-1)*12+8:(t-1)*12+10));
Nino3pr(t,1)=nanmean(Nino3_w((t-1)*12+8:(t-1)*12+10));
Nino12pr(t,1)=nanmean(Nino12_w((t-1)*12+8:(t-1)*12+10));
NZIpr(t,1)=nanmean(NZI_w((t-1)*12+8:(t-1)*12+10));
EPIpr(t,1)= nanmean(EPI_w((t-1)*12+8:(t-1)*12+10));
%TNIpr(t,1)=nanmean(TNI_w((t-1)*12+8:(t-1)*12+10));


TNApr(t,1)=nanmean(TNA_w((t-1)*12+8:(t-1)*12+10));
TSApr(t,1)=nanmean(TSA_w((t-1)*12+8:(t-1)*12+10));
TASIpr(t,1)=nanmean(TASI_w((t-1)*12+8:(t-1)*12+10));
AMOpr(t,1)=nanmean(AMO_w((t-1)*12+8:(t-1)*12+10));

WTIOpr(t,1)=nanmean(WTIO_w((t-1)*12+8:(t-1)*12+10));
SETIOpr(t,1)=nanmean(SETIO_w((t-1)*12+8:(t-1)*12+10));
SWIOpr(t,1)=nanmean(SWIO_w((t-1)*12+8:(t-1)*12+10));
%DMIpr(t,1)=nanmean(DMI_w((t-1)*12+8:(t-1)*12+10));
end

%% Loading other SST indices from NOAA PSL https://www.esrl.noaa.gov/psd/data/climateindices/list/

% change the directory if necessary
PDO=load(['Indices/PDO_1948_2018.txt']); % PDO is not calculated as a simple spatial average, so we will be using the one from NOAA PSL --> https://www.esrl.noaa.gov/psd/data/climateindices/list/

% keep only 01/1967- 12/2016 (study's focus)
% modify accordingly if necessary
PDO(1:19,:)=[];
PDO(end-1:end,:)=[];
% we focus on the months Aug-Sep 
PDOpr=nanmean(PDO(:,9:11),2); %Aug-Oct % 1st column is the year

%% Building predictor matrix
% build the matrix that carries all predictors: PACIFIC, ATLANTIC, INDIAN 
X=[Nino4pr,Nino34pr,Nino3pr,Nino12pr,NZIpr,EPIpr,PDOpr,TNApr,TSApr,AMOpr,TASIpr,WTIOpr,SETIOpr,SWIOpr];

%% Detrending all indices: The linearly detrended series will be used for the analysis
close all
for i=1:size(X,2)
    
X2(:,i)=detrend(X(:,i))+mean(X(:,i)); 

figure()
autocorr(X2(:,i)) % as shown the majortiy of the indices show no significant autocorrelation   
end

% example of detrending
figure()
plot([1967:2016]',X(:,5),'b')
hold on 
plot([1967:2016]',X2(:,5),'k')
hold on 
plot([1967:2016]',X(:,5)-detrend(X(:,5)),'--m')

%% Basic Statistics and Distributions

close all

% Check the Gamma assumption for precipitation series
pprecip=nan*NMprecip2(:,:,1);
negatives=nan*NMprecip2(:,:,1);

for i=1:length(lat_prcp)
    for j=1:length(lon_prcp)
        
        temp=squeeze(NMprecip2(i,j,:));
        
        if sum(isnan(temp))==0 && sum(temp<0)==0
        negatives(i,j)=0;    
        
        phat=gamfit(temp); %fit a gamma distribution 
        [h,pprecip(i,j)]=kstest(norminv(gamcdf(temp,phat(1),phat(2)),0,1)); %test kolmogorov-Smirnov for a gamma
        
        % un-comment to use chi square statistic
        %[h,pprecip(i,j)] = chi2gof(norminv(gamcdf(temp,phat(1),phat(2)),0,1));
        
        end
        
        if sum(isnan(temp))==0 && sum(temp<0)>0 % there are some negative values in the series due to the detrending
        negatives(i,j)=1;
        NMprecip2(i,j,:)=NMprecip2(i,j,:)+abs(min(temp))+1; % we add a constant to the whole series so there is no negative values and we can fit the gamma distribution (this does not affect the analysis of predictability)
        clear temp
        %now perform the test to the shifted sample
        temp=squeeze(NMprecip2(i,j,:));
        phat=gamfit(temp); %fit a gamma distribution 
        [h,pprecip(i,j)]=kstest(norminv(gamcdf(temp,phat(1),phat(2)),0,1)); %test kolmogorov-Smirnov for a gamma
        
        end
        
        clear temp phat h 
    end
end

%check to see where the Gamma distribution is an acceptable model 
% plot the p-value
figure()
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,pprecip,0:0.1:1,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
caxis([0 1])
contourcbar

%check to see in what fraction of the points the p-value for a Gamma distribution is smaller than 0.05
sum(sum(pprecip<0.05))/sum(sum(pprecip>0))
% the result is 0.0011 which is much smaller 0.05 so the Gamma assumption
% is a good model

%hist of the p-values
figure()
hist(reshape(pprecip,[],1))

%check to see in which points we needed to add a constant value 
figure()
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
geoshow(ax, states)
hold on
contourfm(lat_prcp,lon_prcp,negatives,0:0.5:1,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
caxis([0 1])
contourcbar

%check to see in how manay points we needed to add a constant value
sum(sum(negatives==1))
% only in 15 points we needed to add a constant value to the precip series
% that is a very small number. Also, recall that adding a constant does not
% affect covariability and our predictability assessmnet; it is just a
% data shift
%% Check the Normal assumption for the predictors' series
for i=1:size(X2,2)
   
[hh(i),pp(i)]=kstest((X2(:,i)-mean(X2(:,i)))/std(X2(:,i)));

end

% in all cases, the normal assumption is not rejected at 0.05 level (all hh
% values are zeros)

%% Principal Component Analysis of the indices

X2_temp=(X2-repmat(mean(X2,1),size(X2,1),1))./repmat(std(X2),size(X2,1),1);
[wcoeff,~,~,~,explained] = pca(X2_temp)
% the first 2 PCs account for about 60% of the variance 
% and the first 5 PCs acount for about 85% of the variance
X2_pcs=X2_temp*wcoeff;

figure() % not used in the study %
plot(X2_pcs)

%% covariance matrix of the PCs and the indices

%%%%%%%%%%%%%%%%%% Figure 3a of the study %%%%%%%%%%%%%%%%%%
figure() 
CorrMat=corrcoef([X2_pcs,X2]);
image(CorrMat(1:14,15:end),'CDataMapping','scaled')
colorbar
[xxx,yyy]=meshgrid([1:size(X2,2)],[1:size(X2,2)]);
hold on
plot(xxx(abs(CorrMat(1:14,15:end))>=0.278),yyy(abs(CorrMat(1:14,15:end))>=0.278),'x')
clear xxx yyy
caxis([-1 1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot the PCA patterns of the indices 

% focus on Aug-Oct
for t=1:size(monthlySST,3)/12
SST(:,:,t)=nanmean(monthlySST(:,:,(t-1)*12+8:(t-1)*12+10),3); 
end


for i=1:length(lat_sst)
    for j=1:length(lon_sst)
        temp=squeeze(SST(i,j,:));
        temp=detrend(temp); %linearly detrend
        patternPC(i,j,1)=corr(temp,X2_pcs(:,1));
        patternPC(i,j,2)=corr(temp,X2_pcs(:,2));
        patternPC(i,j,3)=corr(temp,X2_pcs(:,3));
        patternPC(i,j,4)=corr(temp,X2_pcs(:,4));
        patternPC(i,j,5)=corr(temp,X2_pcs(:,5));
        clear temp 
    end
end

%%%%%%%%%%%%%%%% Figure 3 b-f %%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:5
figure()
temp=patternPC(:,:,i);
contourf(LON_sst,LAT_sst,temp)
hold on 
plot(coastlon,coastlat,'.')
hold on 
plot(LON_sst(abs(temp)>0.278),LAT_sst(abs(temp)>0.278),'.')
axis([0 360 -70 70])
caxis([-0.8 0.8])
colorbar
clear temp
end

clear SST
clear patternPC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Prediction of precipitation (using simple regression): NO cross-validation

close all
clc

% we use the label "fake" in these variables because this is a within
% sample prediction and does not represent the real predictive skill of the
% model; the effect of a possible overfitting cannot be detected in a within-sample prediction. 

Fake_r_Nino34=nan*NMprecip2(:,:,1);
FakeR2_Nino34=nan*NMprecip2(:,:,1);
FakeR2_Nino34_NZI=nan*NMprecip2(:,:,1);
FakeR2_Pac=nan*NMprecip2(:,:,1);
FakeR2_Pac_Atl=nan*NMprecip2(:,:,1);
FakeR2_all=nan*NMprecip2(:,:,1);

for ii=1:length(lat_prcp)
    ii/length(lat_prcp) % tracking progress
    for jj=1:length(lon_prcp)
        mean_depth2=squeeze(NMprecip2(ii,jj,:));
        if sum(isnan(mean_depth2))==0 && sum(mean_depth2<0)==0    
            
            Fake_r_Nino34(ii,jj)=corr(mean_depth2,X2(:,2));
             
            % linear regression; least squares solution 
            b = regress(mean_depth2,[ones(length(mean_depth2),1),X2(:,2)]);
            FakeR2_Nino34(ii,jj)=1-sum(([ones(length(mean_depth2),1),X2(:,2)]*b-mean_depth2).^2)/sum((mean_depth2-mean(mean_depth2)).^2); 
            clear b 
            
            b = regress(mean_depth2,[ones(length(mean_depth2),1),X2(:,[2,5])]);
            FakeR2_Nino34_NZI(ii,jj)=1-sum(([ones(length(mean_depth2),1),X2(:,[2,5])]*b-mean_depth2).^2)/sum((mean_depth2-mean(mean_depth2)).^2); 
            clear b 
            
            b = regress(mean_depth2,[ones(length(mean_depth2),1),X2(:,[1:7])]);
            FakeR2_Pac(ii,jj)=1-sum(([ones(length(mean_depth2),1),X2(:,[1:7])]*b-mean_depth2).^2)/sum((mean_depth2-mean(mean_depth2)).^2); 
            clear b
            
            b = regress(mean_depth2,[ones(length(mean_depth2),1),X2(:,[1:11])]);
            FakeR2_Pac_Atl(ii,jj)=1-sum(([ones(length(mean_depth2),1),X2(:,[1:11])]*b-mean_depth2).^2)/sum((mean_depth2-mean(mean_depth2)).^2); 
            clear b
            
            b = regress(mean_depth2,[ones(length(mean_depth2),1),X2]);
            FakeR2_all(ii,jj)=1-sum(([ones(length(mean_depth2),1),X2]*b-mean_depth2).^2)/sum((mean_depth2-mean(mean_depth2)).^2); 
            clear b
            
        end
        clear mean_depth2
    end
end

%% plotting results

figure() % not used in the study
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,Fake_r_Nino34,-1:0.1:1,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
caxis([-1 1])
colormap('jet')
contourcbar

%%%%%%%%%%%%% Figure 5b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
k=1;
R2_critical=1/(1+1/finv(0.95,k,50-k-1)*((50-k-1)/k));
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,FakeR2_Nino34,0:0.01:0.5,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
plotm(reshape(LAT_prcp(FakeR2_Nino34>R2_critical),[],1),reshape(LON_prcp(FakeR2_Nino34>R2_critical),[],1),'k.')    
caxis([0 0.5])
colormap('jet')
contourcbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure() % not used in the study
k=2;
R2_critical=1/(1+1/finv(0.95,k,50-k-1)*((50-k-1)/k));
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,FakeR2_Nino34_NZI,0:0.01:0.5,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
plotm(reshape(LAT_prcp(FakeR2_Nino34_NZI>R2_critical),[],1),reshape(LON_prcp(FakeR2_Nino34_NZI>R2_critical),[],1),'k.')    
caxis([0 0.5])
colormap('jet')
contourcbar

%%%%%%%%% Figure 5c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
k=14;
R2_critical=1/(1+1/finv(0.95,k,50-k-1)*((50-k-1)/k));
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,FakeR2_all,0:0.01:0.5,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
plotm(reshape(LAT_prcp(FakeR2_all>R2_critical),[],1),reshape(LON_prcp(FakeR2_all>R2_critical),[],1),'k.')    
caxis([0 0.5])
colormap('jet')
contourcbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
hist(reshape(FakeR2_Nino34_NZI-FakeR2_Nino34,[],1))
figure()
hist(reshape(FakeR2_all-FakeR2_Nino34_NZI,[],1))

%% Prediction of precipitation (using simple regression): 5 fold cross-validation

close all
clc
n_periods=5; % 5-fold cross validation
X2_pcs_std=(X2_pcs-repmat(mean(X2_pcs,1),size(X2_pcs,1),1))./repmat(std(X2_pcs),size(X2_pcs,1),1); % standardized PC series


R2=repmat(nan*NMprecip2(:,:,1),1,1,10);
b2PCs=repmat(nan*NMprecip2(:,:,1),1,1,n_periods,2);

R2best2=nan*NMprecip2(:,:,1);
idbest2=repmat(nan*NMprecip2(:,:,1),1,1,2);
bbest2=repmat(nan*NMprecip2(:,:,1),1,1,n_periods,2);

for ii=1:length(lat_prcp)
    ii/length(lat_prcp) % tracking progress
    for jj=1:length(lon_prcp)
        mean_depth2=squeeze(NMprecip2(ii,jj,:));
        if sum(isnan(mean_depth2))==0 && sum(mean_depth2<0)==0               
                Pest=nan*ones(4,5,size(NMprecip2,3));
                best=nan*ones(4,5,n_periods,2);
                
                l_period=ceil(length(mean_depth2)/n_periods);

                for ff=1:n_periods 
                    test_id=[(ff-1)*l_period+1:(ff-1)*l_period+l_period];   test_id=unique(min(test_id, length(mean_depth2)));
                    train_id=[1:length(mean_depth2)];   train_id(test_id)=[];

                    
                    % linear regression least squares 
                    % using indices
                    b = regress(mean_depth2(train_id),[ones(length(train_id),1),X2(train_id,2)]); % fit model
                    Pest_Nino34(test_id,1)=[ones(length(test_id),1),X2(test_id,2)]*b; % predict out-of-sample
                    clear b  
                    b = regress(mean_depth2(train_id),[ones(length(train_id),1),X2(train_id,[2,5])]);
                    Pest_Nino34_NZI(test_id,1)=[ones(length(test_id),1),X2(test_id,[2,5])]*b;
                    clear b
                    b = regress(mean_depth2(train_id),[ones(length(train_id),1),X2(train_id,[1:7])]);
                    Pest_Pac(test_id,1)=[ones(length(test_id),1),X2(test_id,[1:7])]*b;
                    clear b
                    b = regress(mean_depth2(train_id),[ones(length(train_id),1),X2(train_id,[1:11])]);
                    Pest_Pac_Atl(test_id,1)=[ones(length(test_id),1),X2(test_id,[1:11])]*b;
                    clear b
                    b = regress(mean_depth2(train_id),[ones(length(train_id),1),X2(train_id,:)]);
                    Pest_all(test_id,1)=[ones(length(test_id),1),X2(test_id,:)]*b;
                    clear b 
                    % using PCs
                    b = regress(mean_depth2(train_id),[ones(length(train_id),1),X2_pcs_std(train_id,1)]);
                    Pest_pca1(test_id,1)=[ones(length(test_id),1),X2_pcs_std(test_id,1)]*b;
                    clear b 
                    b = regress(mean_depth2(train_id),[ones(length(train_id),1),X2_pcs_std(train_id,1:2)]);
                    Pest_pca2(test_id,1)=[ones(length(test_id),1),X2_pcs_std(test_id,1:2)]*b;
                    b2PCs(ii,jj,ff,:)=b(2:end);
                    clear b 
                    b = regress(mean_depth2(train_id),[ones(length(train_id),1),X2_pcs_std(train_id,1:3)]);
                    Pest_pca3(test_id,1)=[ones(length(test_id),1),X2_pcs_std(test_id,1:3)]*b;
                    clear b
                    b = regress(mean_depth2(train_id),[ones(length(train_id),1),X2_pcs_std(train_id,1:4)]);
                    Pest_pca4(test_id,1)=[ones(length(test_id),1),X2_pcs_std(test_id,1:4)]*b;
                    clear b
                    b = regress(mean_depth2(train_id),[ones(length(train_id),1),X2_pcs_std(train_id,1:5)]);
                    Pest_pca5(test_id,1)=[ones(length(test_id),1),X2_pcs_std(test_id,1:5)]*b;
                    clear b
                    
                    % searching for the best 2-predictor model using the
                    % first 5 PCS
                    % we build all possible combinations: 5*4/2=10 couples
                
                    for kk=1:4
                        for ll=kk+1:5
                           b = regress(mean_depth2(train_id),[ones(length(train_id),1),X2_pcs_std(train_id,[kk,ll])]);
                           Pest(kk,ll,test_id)=[ones(length(test_id),1),X2_pcs_std(test_id,[kk,ll])]*b;
                           best(kk,ll,ff,:)=b(2:end);
                           clear b 
                        end
                    end
                    
                clear test_id train_id 
                
                end
                
                % here we detect the best performing combination for the
                % specific grid point we are considering in each loop
                R22=1-sum((Pest-repmat(permute(mean_depth2,[3,2,1]),4,5,1)).^2,3)/sum((mean_depth2-mean(mean_depth2)).^2);
                [mR22,idy]=max(R22);
                [R2best2(ii,jj),idx]=max(mR22);
                idbest2(ii,jj,1)=idy(idx); idbest2(ii,jj,2)=idx;
                bbest2(ii,jj,:,:)=best(idy(idx),idx,:,:);
                
                R2(ii,jj,1)=1-sum((Pest_Nino34-mean_depth2).^2)/sum((mean_depth2-mean(mean_depth2)).^2);
                R2(ii,jj,2)=1-sum((Pest_Nino34_NZI-mean_depth2).^2)/sum((mean_depth2-mean(mean_depth2)).^2);              
                R2(ii,jj,3)=1-sum((Pest_Pac-mean_depth2).^2)/sum((mean_depth2-mean(mean_depth2)).^2);
                R2(ii,jj,4)=1-sum((Pest_Pac_Atl-mean_depth2).^2)/sum((mean_depth2-mean(mean_depth2)).^2);
                R2(ii,jj,5)=1-sum((Pest_all-mean_depth2).^2)/sum((mean_depth2-mean(mean_depth2)).^2);
                R2(ii,jj,6)=1-sum((Pest_pca1-mean_depth2).^2)/sum((mean_depth2-mean(mean_depth2)).^2);
                R2(ii,jj,7)=1-sum((Pest_pca2-mean_depth2).^2)/sum((mean_depth2-mean(mean_depth2)).^2);
                R2(ii,jj,8)=1-sum((Pest_pca3-mean_depth2).^2)/sum((mean_depth2-mean(mean_depth2)).^2);
                R2(ii,jj,9)=1-sum((Pest_pca4-mean_depth2).^2)/sum((mean_depth2-mean(mean_depth2)).^2);
                R2(ii,jj,10)=1-sum((Pest_pca5-mean_depth2).^2)/sum((mean_depth2-mean(mean_depth2)).^2);
                
                clear Pest_Nino34 Pest_Nino34_NZI Pest_Pac Pest_Pac_Atl Pest_all Pest_pca1 Pest_pca2 Pest_pca3 Pest_pca4 Pest_pca5
                clear Pest best R22 mR22 idy idx
        end
        clear mean_depth2
    end
end

%% plotting results 

%%%%%%%%%%%% Figure 5d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
k=1;
R2_critical=1/(1+1/finv(0.95,k,50-k-1)*((50-k-1)/k));
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,R2(:,:,1),0:0.01:0.5,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
plotm(reshape(LAT_prcp(R2(:,:,1)>R2_critical),[],1),reshape(LON_prcp(R2(:,:,1)>R2_critical),[],1),'k.')    
caxis([0 0.5])
colormap('jet')
contourcbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure() % not used in the study
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,R2(:,:,2),0:0.01:0.5,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
caxis([0 0.5])
colormap('jet')
contourcbar

figure() % not used in the study
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,R2(:,:,3),0:0.01:0.5,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
caxis([0 0.5])
colormap('jet')
contourcbar

figure() %not used in the study
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,R2(:,:,4),0:0.01:0.5,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
caxis([0 0.5])
colormap('jet')
contourcbar

%%%%%%%%%%%%% Figure 5e %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
k=14;
R2_critical=1/(1+1/finv(0.95,k,50-k-1)*((50-k-1)/k));
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,R2(:,:,5),0:0.01:0.5,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
%plotm(reshape(LAT_prcp(R2(:,:,5)>R2_critical),[],1),reshape(LON_prcp(R2(:,:,5)>R2_critical),[],1),'k.')    
caxis([0 0.5])
colormap('jet')
contourcbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Figure 6b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
k=1;
R2_critical=1/(1+1/finv(0.95,k,50-k-1)*((50-k-1)/k));
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,R2(:,:,6),0:0.01:0.5,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
plotm(reshape(LAT_prcp(R2(:,:,6)>R2_critical),[],1),reshape(LON_prcp(R2(:,:,6)>R2_critical),[],1),'k.') 
caxis([0 0.5])
colormap('jet')
contourcbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Figure 6c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
k=2;
R2_critical=1/(1+1/finv(0.95,k,50-k-1)*((50-k-1)/k));
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,R2(:,:,7),0:0.01:0.5,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
plotm(reshape(LAT_prcp(R2(:,:,7)>R2_critical),[],1),reshape(LON_prcp(R2(:,:,7)>R2_critical),[],1),'k.') 
caxis([0 0.5])
colormap('jet')
contourcbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Figure 6d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
k=3;
R2_critical=1/(1+1/finv(0.95,k,50-k-1)*((50-k-1)/k));
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,R2(:,:,8),0:0.01:0.5,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
plotm(reshape(LAT_prcp(R2(:,:,8)>R2_critical),[],1),reshape(LON_prcp(R2(:,:,8)>R2_critical),[],1),'k.') 
caxis([0 0.5])
colormap('jet')
contourcbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% Figure 6e %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure()
k=2;
R2_critical=1/(1+1/finv(0.95,k,50-k-1)*((50-k-1)/k));
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,R2best2,0:0.01:0.5,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
plotm(reshape(LAT_prcp(R2best2>R2_critical),[],1),reshape(LON_prcp(R2best2>R2_critical),[],1),'k.') 
caxis([0 0.5])
colormap('jet')
contourcbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% plots for overfitting

%%%%%%%% Figure 5a %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
plot([1:4]',[nanmean(reshape(FakeR2_Nino34,[],1)),nanmean(reshape(FakeR2_Pac,[],1)),nanmean(reshape(FakeR2_Pac_Atl,[],1)),nanmean(reshape(FakeR2_all,[],1))],'-o')
hold on
%errorbar([1:4]',[nanmean(reshape(FakeR2_Nino34,[],1)),nanmean(reshape(FakeR2_Pac,[],1)),nanmean(reshape(FakeR2_Pac_Atl,[],1)),nanmean(reshape(FakeR2_all,[],1))],...
%    [nanstd(reshape(FakeR2_Nino34,[],1)),nanstd(reshape(FakeR2_Pac,[],1)),nanstd(reshape(FakeR2_Pac_Atl,[],1)),nanstd(reshape(FakeR2_all,[],1))])
hold on
plot([1:4]',[nanmean(reshape(R2(:,:,1),[],1)),nanmean(reshape(R2(:,:,3),[],1)),nanmean(reshape(R2(:,:,4),[],1)),nanmean(reshape(R2(:,:,5),[],1))],'-o')
hold on
%errorbar([1:4]',[nanmean(reshape(R2(:,:,1),[],1)),nanmean(reshape(R2(:,:,3),[],1)),nanmean(reshape(R2(:,:,4),[],1)),nanmean(reshape(R2(:,:,5),[],1))],...
%    [nanstd(reshape(R2(:,:,1),[],1)),nanstd(reshape(R2(:,:,3),[],1)),nanstd(reshape(R2(:,:,4),[],1)),nanstd(reshape(R2(:,:,5),[],1))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%% Figure 6a %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure()
plot([1:5]',[nanmean(reshape(R2(:,:,6),[],1)),nanmean(reshape(R2(:,:,7),[],1)),nanmean(reshape(R2(:,:,8),[],1)),nanmean(reshape(R2(:,:,9),[],1)),nanmean(reshape(R2(:,:,10),[],1))],'-o')
hold on
%errorbar([1:5]',[nanmean(reshape(R2(:,:,6),[],1)),nanmean(reshape(R2(:,:,7),[],1)),nanmean(reshape(R2(:,:,8),[],1)),nanmean(reshape(R2(:,:,9),[],1)),nanmean(reshape(R2(:,:,10),[],1))],...
%    [nanstd(reshape(R2(:,:,6),[],1)),nanstd(reshape(R2(:,:,7),[],1)),nanstd(reshape(R2(:,:,8),[],1)),nanstd(reshape(R2(:,:,9),[],1)),nanstd(reshape(R2(:,:,10),[],1))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Which PC component was the most important for prediction?

import_max_best=nan*NMprecip2(:,:,1);
id_max_best=nan*NMprecip2(:,:,1:2);

for ii=1:length(lat_prcp)
    for jj=1:length(lon_prcp)
        mean_depth2=squeeze(NMprecip2(ii,jj,:));
        if sum(isnan(mean_depth2))==0 && sum(mean_depth2<0)==0  
        tempb=permute(bbest2(ii,jj,:,:),[3,4,1,2]);
        
        temp1=nanmean(abs(tempb)./repmat(sum(abs(tempb),2),1,2),1);
        clear tempb
        tempb=temp1; clear temp1
        [import_max_best(ii,jj),temp1]=max(tempb);
        id_max_best(ii,jj,1)=idbest2(ii,jj,temp1); % the most important component
        id_max_best(ii,jj,2)=idbest2(ii,jj,(temp1==1)+1); % the second most important component
        
        clear temp1
        
        clear tempb    
        end
    end
end

% Plotting results
close all
k=2;
R2_critical=1/(1+1/finv(0.95,k,50-k-1)*((50-k-1)/k));

%%%%%%%%%%%% Figure 9a %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
%temp=id_max_best(:,:,1);
temp=id_max_best(:,:,1).*(R2best2>R2_critical);
temp(temp==0)=nan;
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,temp,0.5:5.5,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
caxis([0.5 5.5])
contourcbar
clear temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Figure 9b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure()
%temp=id_max_best(:,:,2);
temp=id_max_best(:,:,2).*(R2best2>R2_critical);
temp(temp==0)=nan;
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,temp,0.5:5.5,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
caxis([0.5 5.5])
contourcbar
clear temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure() % not used in the study
%temp=import_max_best*100;
temp=import_max_best*100.*(R2best2>R2_critical);
temp(temp==0)=nan;
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,temp,0:10:100,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
caxis([0 100])
contourcbar
clear temp

%% Prediction of precipitation using COPULAS: using a 5 fold cross-validation

% this code cell takes a while to finish running
% be aware!
% the results of this cel are already saved in the mat file "USresults.mat"

close all
clc
n_periods=5; % 5-fold cross validation


R2=repmat(nan*NMprecip2(:,:,1),1,1,3);
CSI_dry=repmat(nan*NMprecip2(:,:,1),1,1,3,3);
CSI_wet=repmat(nan*NMprecip2(:,:,1),1,1,3,3);
CSI_nor=repmat(nan*NMprecip2(:,:,1),1,1,3,3);
LogL=repmat(nan*NMprecip2(:,:,1),1,1,3,3);


for ii=1:length(lat_prcp)
    ii/length(lat_prcp)
    for jj=1:length(lon_prcp)
        mean_depth2=squeeze(NMprecip2(ii,jj,:));
        if sum(isnan(mean_depth2))==0 && sum(mean_depth2<0)==0               

                x=[0:2*max(mean_depth2)/100:2*max(mean_depth2)];
                l_period=ceil(length(mean_depth2)/n_periods);

                for ff=1:n_periods 
                    test_id=[(ff-1)*l_period+1:(ff-1)*l_period+l_period];   test_id=unique(min(test_id, length(mean_depth2)));
                    train_id=[1:length(mean_depth2)];   train_id(test_id)=[];

                    %fitting distributions
                    phat=gamfit(mean_depth2(train_id)); % precip
                    for i=1:size(X2,2)
                    [muHat(i),sigmaHat(i)] = normfit(X2(train_id,i)); % indices
                    end
                    for i=1:size(X2_pcs,2)
                    [muHatpc(i),sigmaHatpc(i)] = normfit(X2_pcs(train_id,i)); %principal components
                    end
                    
                    for t=1:length(test_id)
                           for i=1:size(X2,2)
                                ecdf(i) = normcdf(X2(test_id(t),i),muHat(i),sigmaHat(i));
                           end
                           for i=1:size(X2_pcs,2)
                                ecdfpc(i) = normcdf(X2_pcs(test_id(t),i),muHatpc(i),sigmaHatpc(i));
                           end

                           for i=1:length(x)

                                PRepdf=gampdf(x(i),phat(1),phat(2));
                                PRecdf=gamcdf(x(i),phat(1),phat(2));    


                                %predicting
                                
                                %for j=1:size(X2,2)
                                %PRpdf(i,test_id(t),j) = copulapdf('Gaussian',[PRecdf,ecdf(j)],corr(mean_depth2(train_id),X2(train_id,j)))*PRepdf;   
                                %end
                                
                                
                                % COPULA PREDICTION MODEL
                                PRpdf_Nino34(i,test_id(t)) = copulapdf('Gaussian',[PRecdf,ecdf(2)],corr(mean_depth2(train_id),X2(train_id,2)))*PRepdf;
                                PRpdf_pca2(i,test_id(t))= copulapdf('Gaussian',[PRecdf,ecdfpc(1:2)],nearcorr(corrcoef([mean_depth2(train_id),X2_pcs(train_id,[1,2])])))/copulapdf('Gaussian',ecdfpc(1:2),nearcorr(corrcoef([X2_pcs(train_id,[1,2])])))*PRepdf;
                                PRpdf_pca_best2(i,test_id(t))= copulapdf('Gaussian',[PRecdf,ecdfpc(idbest2(ii,jj,:))],nearcorr(corrcoef([mean_depth2(train_id),X2_pcs(train_id,idbest2(ii,jj,:))])))/copulapdf('Gaussian',ecdfpc(idbest2(ii,jj,:)),nearcorr(corrcoef([X2_pcs(train_id,idbest2(ii,jj,:))])))*PRepdf;

                                clear PRepdf PRecdf
                           end
                           
                           clear ecdf ecdfpc
                           
                    end


                clear test_id train_id phat muHat sigmaHat muHatpc sigmaHatpc
                
                end
                
                % coef of determination
                R2(ii,jj,1)=1-sum((sum(PRpdf_Nino34.*repmat(x',1,length(mean_depth2))*(x(2)-x(1)),1)'-mean_depth2).^2)/sum((mean_depth2-mean(mean_depth2)).^2);
                R2(ii,jj,2)=1-sum((sum(PRpdf_pca2.*repmat(x',1,length(mean_depth2))*(x(2)-x(1)),1)'-mean_depth2).^2)/sum((mean_depth2-mean(mean_depth2)).^2);
                R2(ii,jj,3)=1-sum((sum(PRpdf_pca_best2.*repmat(x',1,length(mean_depth2))*(x(2)-x(1)),1)'-mean_depth2).^2)/sum((mean_depth2-mean(mean_depth2)).^2);
                
                
                PDFs=cat(3,PRpdf_Nino34,PRpdf_pca2,PRpdf_pca_best2);
                
                % calculating categorical statistics
                
                for i=1:size(PDFs,3)

                    pdf=PDFs(:,:,i);
                    pdf=pdf+10^(-10);

                    cdf=cumsum(pdf,1)*(x(2)-x(1));

                    % Randomly sample 1000 realizations from the predictive
                    % distributions.
                    % Modify if you want to sample more realizations.
                    % Results will be more accurate but more
                    % time-consuming too.
                    for j=1:length(mean_depth2)
                    real1(:,j)=interp1(cdf(:,j),x',rand(1000,1));
                    end 

                    %Wet and Dry Hits
                    j=0;
                    for p=[1/3,0.25,0.15]
                    j=j+1;
                    Q1=quantile(mean_depth2,p);
                    Q2=quantile(mean_depth2,1-p);

                    Idry=find(mean_depth2<=Q1);
                    Iwet=find(mean_depth2>Q2);
                    Inorm=find((mean_depth2>Q1).*(mean_depth2<=Q2));

                    temp=sum(sum(real1(:,Idry)<=Q1))/sum(sum(real1(:,Idry)>-1)); %dry detection
                    temp1=temp*p/(sum(sum(real1<=Q1))/sum(sum(real1>-1))); %dry success
                    CSI_dry(ii,jj,i,j)=1/((1/temp)+(1/temp1)-1); %critical success index dry
                    clear temp temp1

                    temp=sum(sum(real1(:,Iwet)>Q2))/sum(sum(real1(:,Iwet)>-1)); %wet detection
                    temp1=temp*p/(sum(sum(real1>Q2))/sum(sum(real1>-1))); %wet success
                    CSI_wet(ii,jj,i,j)=1/((1/temp)+(1/temp1)-1); %critical success index wet
                    clear temp temp1

                    temp=sum(sum((real1(:,Inorm)>Q1).*(real1(:,Inorm)<=Q2)))/sum(sum(real1(:,Inorm)>-1)); %normal detection
                    temp1=temp*(1-2*p)/(sum(sum((real1>Q1).*(real1<=Q2)))/sum(sum(real1>-1))); %normal success
                    CSI_nor(ii,jj,i,j)=1/((1/temp)+(1/temp1)-1); %critical success index normal
                    clear temp temp1

                    % Categorical loglikelihood
                    ID=(mean_depth2<=Q1)*(-1)+(mean_depth2>Q2)*1;
                    temp=prod([sum(real1(:,ID==-1)<=Q1,1),sum((real1(:,ID==0)>Q1).*(real1(:,ID==0)<=Q2),1),sum(real1(:,ID==1)>Q2,1)]/size(real1,1));
                    LogL(ii,jj,i,j)= exp(log(temp)/length(mean_depth2));
                    clear temp

                    clear Q1 Q2 Idry Iwet Inorm ID
                    end

                    clear pdf cdf real1  
                    
                end
                
                clear x PDFs PRpdf_Nino34 PRpdf_pca2 PRpdf_pca_best2
        end
        clear mean_depth2
    end
end



%% saving the results
save('Wprecip prediction results_US_COBEv2/USresults.mat','R2','CSI_dry','CSI_wet','CSI_nor','LogL','lat_prcp','lon_prcp')

%% loading the results
clear
close all
clc
load('Wprecip prediction results_US_COBEv2/stats_limits_MC.mat') % loading results from Monte Carlo
load('Wprecip prediction results_US_COBEv2/USresults.mat') % loading results from copula prediction across CONUS
[LON_prcp,LAT_prcp]=meshgrid(lon_prcp,lat_prcp);

%% plotting the R2 results: 
% These are exactly the same results with the ones that we produced when we used least-squares prediction in lines 879-1113
% Any slight difference might be attributed to imperfect model fit, but in
% theory these are the same things. 

figure()
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,R2(:,:,1),0:0.01:0.5,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
caxis([0 0.5])
contourcbar

figure()
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,R2(:,:,2),0:0.01:0.5,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
caxis([0 0.5])
contourcbar

figure()
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
contourfm(lat_prcp,lon_prcp,R2(:,:,3),0:0.01:0.5,'LineStyle','none') 
geoshow(ax, states,'FaceColor','none')
caxis([0 0.5])
contourcbar

%% Plotting Results for Categorical metrics 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                           %
% These plots are presented in Figures 7 and 8 in the study %
%                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


p=[1/3,0.25,0.15];

i=1;
figure()
    temp=CSI_dry(:,:,3,i);
    ax = usamap('conus');
    states = shaperead('usastatelo', 'UseGeoCoords', true,...
       'Selector',...
       {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
    contourfm(lat_prcp,lon_prcp,temp,p(i)/(2-p(i)):0.01:0.45,'LineStyle','none') 
    geoshow(ax, states,'FaceColor','none')
    hold on 
    plotm(reshape(LAT_prcp(temp>quantile(table1(1,:,i),0.95)),[],1),reshape(LON_prcp(temp>quantile(table1(1,:,i),0.95)),[],1),'k.')
    caxis([p(i)/(2-p(i)) 0.45])
    colormap('jet')
    contourcbar
    clear temp

figure()
    temp=CSI_wet(:,:,3,i);
    ax = usamap('conus');
    states = shaperead('usastatelo', 'UseGeoCoords', true,...
      'Selector',...
      {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
    contourfm(lat_prcp,lon_prcp,temp,p(i)/(2-p(i)):0.01:0.45,'LineStyle','none') 
    geoshow(ax, states,'FaceColor','none')
    hold on 
    plotm(reshape(LAT_prcp(temp>quantile(table1(1,:,i),0.95)),[],1),reshape(LON_prcp(temp>quantile(table1(1,:,i),0.95)),[],1),'k.')
    caxis([p(i)/(2-p(i)) 0.45])
    colormap('jet')
    contourcbar
    clear temp

figure()
    temp=CSI_nor(:,:,3,i);
    ax = usamap('conus');
    states = shaperead('usastatelo', 'UseGeoCoords', true,...
      'Selector',...
      {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
    contourfm(lat_prcp,lon_prcp,temp,(1-2*p(i))/(1+2*p(i)):0.01:0.4,'LineStyle','none') 
    geoshow(ax, states,'FaceColor','none')
    hold on 
    plotm(reshape(LAT_prcp(temp>quantile(table1(2,:,i),0.95)),[],1),reshape(LON_prcp(temp>quantile(table1(2,:,i),0.95)),[],1),'k.')
    caxis([(1-2*p(i))/(1+2*p(i)) 0.4])
    colormap('jet')
    contourcbar
    clear temp

figure()
    temp=LogL(:,:,3,i);
    ax = usamap('conus');
    states = shaperead('usastatelo', 'UseGeoCoords', true,...
      'Selector',...
      {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
    contourfm(lat_prcp,lon_prcp,temp,exp(2*p(i)*log(p(i))+(1-2*p(i))*log(1-2*p(i))):0.01:0.4,'LineStyle','none') 
    geoshow(ax, states,'FaceColor','none')
    hold on 
    plotm(reshape(LAT_prcp(temp>quantile(table1(3,:,i),0.95)),[],1),reshape(LON_prcp(temp>quantile(table1(3,:,i),0.95)),[],1),'k.')
    caxis([exp(2*p(i)*log(p(i))+(1-2*p(i))*log(1-2*p(i))) 0.4])
    colormap('jet')
    contourcbar
    clear temp
    
    
i=2;
figure()
    temp=CSI_dry(:,:,3,i);
    ax = usamap('conus');
    states = shaperead('usastatelo', 'UseGeoCoords', true,...
       'Selector',...
       {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
    contourfm(lat_prcp,lon_prcp,temp,p(i)/(2-p(i)):0.01:0.4,'LineStyle','none') 
    geoshow(ax, states,'FaceColor','none')
    hold on 
    plotm(reshape(LAT_prcp(temp>quantile(table1(1,:,i),0.95)),[],1),reshape(LON_prcp(temp>quantile(table1(1,:,i),0.95)),[],1),'k.')
    caxis([p(i)/(2-p(i)) 0.4])
    colormap('jet')
    contourcbar
    clear temp

figure()
    temp=CSI_wet(:,:,3,i);
    ax = usamap('conus');
    states = shaperead('usastatelo', 'UseGeoCoords', true,...
      'Selector',...
      {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
    contourfm(lat_prcp,lon_prcp,temp,p(i)/(2-p(i)):0.01:0.4,'LineStyle','none') 
    geoshow(ax, states,'FaceColor','none')
    hold on 
    plotm(reshape(LAT_prcp(temp>quantile(table1(1,:,i),0.95)),[],1),reshape(LON_prcp(temp>quantile(table1(1,:,i),0.95)),[],1),'k.')
    caxis([p(i)/(2-p(i)) 0.4])
    colormap('jet')
    contourcbar
    clear temp

figure()
    temp=CSI_nor(:,:,3,i);
    ax = usamap('conus');
    states = shaperead('usastatelo', 'UseGeoCoords', true,...
      'Selector',...
      {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
    contourfm(lat_prcp,lon_prcp,temp,(1-2*p(i))/(1+2*p(i)):0.01:0.45,'LineStyle','none') 
    geoshow(ax, states,'FaceColor','none')
    hold on 
    plotm(reshape(LAT_prcp(temp>quantile(table1(2,:,i),0.95)),[],1),reshape(LON_prcp(temp>quantile(table1(2,:,i),0.95)),[],1),'k.')
    caxis([(1-2*p(i))/(1+2*p(i)) 0.45])
    colormap('jet')
    contourcbar
    clear temp

figure()
    temp=LogL(:,:,3,i);
    ax = usamap('conus');
    states = shaperead('usastatelo', 'UseGeoCoords', true,...
      'Selector',...
      {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
    contourfm(lat_prcp,lon_prcp,temp,exp(2*p(i)*log(p(i))+(1-2*p(i))*log(1-2*p(i))):0.01:0.45,'LineStyle','none') 
    geoshow(ax, states,'FaceColor','none')
    hold on 
    plotm(reshape(LAT_prcp(temp>quantile(table1(3,:,i),0.95)),[],1),reshape(LON_prcp(temp>quantile(table1(3,:,i),0.95)),[],1),'k.')
    caxis([exp(2*p(i)*log(p(i))+(1-2*p(i))*log(1-2*p(i))) 0.45])
    colormap('jet')
    contourcbar
    clear temp

















