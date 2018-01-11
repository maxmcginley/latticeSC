classdef latticeSC_SSH < latticeSC
    %latticeSC_SSH A subclass of latticeSC with added functionality
    %specific to the SSH model 
    
    properties
    end
    
    methods
        function obj = latticeSC_SSH(sites, species, bonds, open)
            obj@latticeSC(sites,species,bonds,open);
        end
        
        %FUNCTION STRING_ORDER
        %Function that evaluates the string order parameter from (site1) to
        %(site2) for an SSH model.
        function val = string_order(obj,site1,site2,occs)
            assert(mod(site1,2) == 0,'String order: site 1 must be even');
            assert(mod(site2,2) == 1,'String order: site 1 must be odd');
            assert(site2 > site1,'String order: site 2 must be after site 1');
            A_inds = (site1+1):2:site2;
            B_inds = site1:2:(site2-1);
            coeff = (-1i)^((site2+1-site1)/2);
            val = coeff*obj.evaluate_operator_expectation(occs,A_inds,B_inds);
        end
    end
    
    methods(Static)
        %FUNCTION CREATE_SSH
        %Creates an instance of the SSH model given hopping amplitudes (tA)
        %and (tB), (sites), and whether the chain is open.
        function obj = create_SSH(tA,tB,stag_mu,sites,open)
            assert(mod(sites,2) == 0,'SSH model: must have even sites');
            hopping_amplitudes = repmat([tA,tB],1,sites/2);
            chem_potential_amplitudes = repmat([stag_mu,-stag_mu],1,sites/2);

            bonds = cell(0);
            bonds = latticeSC.add_to_bonds(bonds,[1,2],hopping_amplitudes);
            bonds = latticeSC.add_to_bonds(bonds,[2,1],hopping_amplitudes);
            bonds = latticeSC.add_to_bonds(bonds,[1,1],chem_potential_amplitudes);

            obj = latticeSC_SSH(sites,1,bonds,open);
        end
    end
end

