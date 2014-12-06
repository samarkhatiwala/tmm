% Annual Net Flux of Carbon to the Atmosphere from Land-Use Change:  1850-2005
%    are taken from Houghton @
% 
% http://cdiac.esd.ornl.gov/trends/landuse/houghton/1850-2005.txt
% Units = Tg C (1 teragram = 10^12 g)
% CITE AS: Houghton, R.A. 2008. Carbon Flux to the Atmosphere from Land-Use Changes: 1850-2005. In TRENDS: A Compendium of Data on Global Change. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, U.S. Department of Energy, Oak Ridge, Tenn., U.S.A.

function [year,lu,units,author]=read_land_use_data(data_folder);

filein=fullfile(data_folder,'global_land_use_1850_2005.dat');

[data]=textread(filein);
year=data(:,1);
lu=data(:,2);%TgC/year=1e3 GtC/year
lu=lu*1e-3;%Converted land use Emissions into [GtC/year]
units='given in TgC/year and transformed in GtC/year';
author='Houghton';


