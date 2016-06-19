function [T_ff,ff,T_lu,lu,T_pCo2,pCo2]=get_emissions_pCO2atm(emissionScenario,emissionDataPath);

%
% Create arrays for past, present and future emissions and pCO2 in the
% atmosphere
% 
% The arrays are constructed as follows:
%
% PAST TO PRESENT
%
% fossil fuel 1751-2005 [PgC/y] from Marland et al, 2007 - type help read_fossil_fuel_data for more info
% land usage  1850-2005 [PgC/y] from Hougthon,  2008     - type help read_land_use_data for more info
% Global atmospheric pCo2 1765-2008 [ppmv] from OCMIP - updated yearly with global mean surface pCO2 from NOAA 
%
% FUTURE
%
% fossil fuel 1990-2100 [PgC/y] from IPCC SRES A2 or B1 (Hougthon) - type help read_scenarios for more info
% land usage  1990-2100 [PgC/y] from IPCC SRES A2 or B1 (Hougthon) - type help read_scenarios for more info
%
% - How past to present real data are treated:
%
% 1) Take as reference preindustrial year the ones for pCo2atm [=1765]- Because the
%    real fossil fuel emissions stops in 2005, take that as final year and disregard 
%    real data from years after 2005 and before 1765
%
% 2) Interpolate land use emissions to zero for years before 1850
% 
% fossil fuel 1765-2005
% land usage  1765-2005
% atmospheric pCo2 1765-2005 
%
% - How future data are treated:
%
% 3) Land-use and Fossil fuel emissions are given every 10 years starting
% from 1990. So first interpolate those data for every year. These data
% are NOT going to match 2005 values (because they are projections made in
% the 90-ies which do not coincide with reality!). So rescale 2005 values
% in order to match *real* data and take date from that year onwwards to
% 2100. NOTE: future emissions depend on what model within the A2 scenarios
% is chosen. This is set with a key. Type read_A2_scenarios for more info
% 
% 4) Future guess for pCo2atm is found assuming a linear increase for the 
%    perturbed pCO2atm, starting from real value at 2005. 
%
%
% - How past to present and future data matched:
%
% Create single arrays matching past to present and future emissions and pCo2 
% 




%key=1;%key for future emissions (type "help read_scenarios" for legend)

Ti=1765;%Initial Tpco2a
Tf=2005;%Final fossil fuel emissions data

%Get atmospheric pCo2 time history up to present 
%load('CO2_atm.mat','year','pCo2');%
load(fullfile(emissionDataPath,'CO2_atm.mat'),'year','pCo2');%
%[data]=load('noces_co2_mlo.txt');
%year=data(:,1);
%pCo2=data(:,2);

index=find(Ti<=year & year<=Tf);
year=year(index);
pCo2=pCo2(index);



%Read fossil fuel and industrial emission scenario past to present
[year1,ff1]=read_fossil_fuel_data(emissionDataPath);
index=find(Ti<=year1 & year1<=Tf);
year1=year1(index);
ff1=ff1(index);

%Read fossil fuel and industrial emission scenario past to present
[year2,lu2]=read_land_use_data(emissionDataPath);
index=find(Ti<=year2 & year2<=Tf);
year2=year2(index);
lu2=lu2(index);
lu2=[0; lu2];
year2=[Ti;year2];
%lu = spline(year2,lu2,year1);

lu2=interp1(year2,lu2,year1,'linear','extrap');%needs to be addressed
year2=year1;

getFuture=0;

if nargin>0
  if ~isempty(emissionScenario)
    getFuture=1;
  end
end

if getFuture


%Read fossil fuel+industrial emission & land use scenario for the future
%[year3,ff3,lu3]=read_A2_scenarios(key);
[year3,ff3,lu3]=read_scenarios(emissionScenario,emissionDataPath);

% Interpolate future emission scenarios to get data each year instead of
% every 10 years
year3i=[year3(1):1:year3(end)];
ff3i=interp1(year3,ff3,year3i);
lu3i=interp1(year3,lu3,year3i);

year3=year3i;
ff3=ff3i;
lu3=lu3i;

%Match data for Tf and pick up only future emissions from Tf onwards
index=find(year==Tf);
index3=find(year3==Tf);
ff3=ff1(index).*ff3./ff3(index3);
index=find(year2==Tf);
lu3=lu2(index).*lu3./lu3(index3);
year3=year3(index3:end)';
ff3=ff3(index3:end)';
lu3=lu3(index3:end)';


%FIRST GUESS for pCo2a in the future

pCo2p=pCo2-pCo2(1);
index=find(year==Tf);
pCo2pobs=pCo2p(index);%perturbation in final year of data array (=Tf)
pCo2pf=1.*(year3-year3(1))+pCo2pobs;%assume a rate of increas of 1ppm/year
pCo2f=pCo2pf+pCo2(1);


%Merge all these data
T_pCo2=[year;year3(2:end)];
pCo2=[pCo2;pCo2f(2:end)];
T_ff=[year1;year3(2:end)];
ff=[ff1;ff3(2:end)];
T_lu=[year2;year3(2:end)];
lu=[lu2;lu3(2:end)];

else

T_pCo2=[year];
pCo2=[pCo2];
T_ff=[year1];
ff=[ff1];
T_lu=[year2];
lu=[lu2];

end
