clear;
clc;

addpath('..');

%************INPUT*************
sites = 50;
tA = 0.5*exp(0.2i);
tB = 1;
stag_mu = 1;
open = true;
site1 = 10;
site2 = 39;
%******************************

sys = latticeSC_SSH.create_SSH(tA,tB,stag_mu,sites,open);

sys.spectrum

sys.string_order(site1,site2,zeros(1,sites))