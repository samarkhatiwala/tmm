function [year,ff,lu,tot,units,author]=read_scenarios(emissionScenario,emissionDataPath);

% Read different IPCC SRES A2 and B1 scenarios from Hougthon
% Data in GtC/year
%
%LEGEND: A2
%
%1=> ASF - what we have used so far
%2=>AIM
%3=> IMAGE - what i think Friedlingstein et al 2006 used
%4=> MESSAGE
%5=> MINICAM
%
%LEGEND: B1
%
%6=> ASF - what we have used so far
%7=>AIM
%8=> IMAGE - what i think Friedlingstein et al 2006 used
%9=> MESSAGE
%10=> MINICAM
%

filein=['future_' emissionScenario '_emissions'];
  
% if key==1
%   filein='future_A2_ASF_emissions.mat';
% end
% if key==2
%   filein='future_A2_AIM_emissions.mat';
% end
% if key==3
%   filein='future_A2_IMAGE_emissions.mat';
% end
% if key==4
%   filein='future_A2_MESSAGE_emissions.mat';
% end
% if key==5
%   filein='future_A2_MINICAM_emissions.mat';
% end
% if key==6
%   filein='future_B1_ASF_emissions.mat';
% end
% if key==7
%   filein='future_B1_AIM_emissions.mat';
% end
% if key==8
%   filein='future_B1_IMAGE_emissions.mat';
% end
% if key==9
%   filein='future_B1_MESSAGE_emissions.mat';
% end
% if key==10
%   filein='future_B1_MINICAM_emissions.mat';
% end

load(fullfile(emissionDataPath,filein))
year=time;
ff=ff;
lu=lu;
tot=tot;
units=units;
author=author;

