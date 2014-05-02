%    Global CO2 Emissions from Fossil-Fuel Burning,      
%    Cement Manufacture, and Gas Flaring: 1751-2004      
%    are taken from Marland et al @:
%    http://cdiac.ornl.gov/trends/emis/tre_glob.htm
%    All emission estimates are expressed in million metric tons of carbon.
%    CITE AS: Marland, G., T.A. Boden, and R. J. Andres. 2007. Global, Regional, and National CO2 Emissions. In Trends: A Compendium of Data on Global Change. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, U.S. Department of Energy, Oak Ridge, Tenn., U.S.A.
%  
%    File=global.1751_2004.ems 
%    in '/Users/fterenzi/work_giss/data/emissions'
  
function [year,ff,units,author]=read_fossil_fuel_data(data_folder);

filein=fullfile(data_folder,'global_emissions_1751_2005.dat');

[data]=textread(filein);
year=data(:,1);
ff=data(:,2);% Million of tons of carbon/year=TgC/year=1e3 GtC/year
ff=ff*1e-3;%Converted fossil fuel and industry Emissions into [GtC/year]
units='Given in Millions on tons/year and transformed in GtC/year';
author='Marland et al 2007';
