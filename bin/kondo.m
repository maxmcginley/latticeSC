clear;
clc;

addpath('..');

%************INPUT*************
sites = 100;
t = 1;
kondo_site = 50;
U_imp = 5;
mu = 0.5;
open = true;
%******************************

mus = ones(1,sites).*mu;
mus(kondo_site) = 0;

AB_ups = ones(1,sites).*t;
AB_ups(kondo_site-1) = 0;

BA_ups = ones(1,sites).*t;
BA_ups(kondo_site) = 0;

imp_terms = zeros(1,sites);
imp_terms(kondo_site) = U_imp;

bonds = cell(0);

bonds = latticeSC.add_to_bonds(bonds,[1,1],mus);
bonds = latticeSC.add_to_bonds(bonds,[2,2],mus);
bonds = latticeSC.add_to_bonds(bonds,[1,3],AB_ups);
bonds = latticeSC.add_to_bonds(bonds,[2,4],BA_ups);
bonds = latticeSC.add_to_bonds(bonds,[3,1],BA_ups);
bonds = latticeSC.add_to_bonds(bonds,[4,2],AB_ups);
bonds = latticeSC.add_to_bonds(bonds,[1,2],imp_terms);

sys = latticeSC(sites,2,bonds,open);

sys.spectrum