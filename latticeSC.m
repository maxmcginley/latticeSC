classdef latticeSC
    %latticeSC - Class for diagonalising and extracting data from 1D
    %non-interacting superconducting lattice models
    %   The class allows for Hamiltonian terms in the form of Majorana
    %   fermion bilinears between flavours A and B.
    
    %     ****************** PROPERTIES **********************
    
    properties
        sites
            %PROPERTY SITES
            %Number of sites in the lattice model
        species
            %PROPERTY SPECIES
            %Number of fermions per site, e.g. 2 for spin-half fermions
        hamiltonian
            %PROPERTY HAMILTONIAN (species * sites) by (species * sites)
            %matrix representing the Hamiltonian in the form H = \gamma^T_A
            %. hamiltonian . \gamma_B mixing only A and B flavour Majorana
            %operators. Indexing is done by site first, then by species.
        spectrum
            %PROPERTY SPECTRUM
            %Energy eigenvalues of the hamiltonian
        U
            %PROPERTY U
            %Left SVD matrix where h = U' * spectrum * V
        V
            %PROPERTY V
            %Right SVD matrix where h = U' * spectrum * V
    end
    
    methods
        function obj = latticeSC(sites, species, bonds, open)
            %latticeSC Construct an instance of this class
            %   Sites and species is provided. Bonding data is provided as
            %   the [] by 2 cell (bonds). Within (bonds), the first column
            %   contains 1 by 2 matrices which define the location of the
            %   bonds up to translation by any number of sites. The second
            %   column contains 1 by (sites) matrices (species) by (sites)
            %   matrices. The first index of this data encodes what species
            %   the bond is between. The order is 1-1, 1-2,  ...,
            %   1-(species), followed by 2-1, 2-2, ... and so on. In the
            %   case of periodic (open = false) BCs, the last column of the
            %   nearest-neighbour data completes the loop from (sites) to
            %   1.
            
            obj.sites = sites;
            obj.species = species;
            obj.hamiltonian = latticeSC.hamiltonian_from_data(...
                sites,species,bonds,open);
            
            [u_tmp, s, v_tmp] = svd(obj.hamiltonian);
            obj.U  = u_tmp';
            obj.V = v_tmp';
            obj.spectrum = diag(s);
        end
        
        %FUNCTION correlation_matrix
        %Gives a matrix of the expectation value of \gamma_a \gamma_b
        %in the eigenstate specified by the occupation numbers occs. 
        function corrmat = correlation_matrix(obj,occs)
            corrmat = 1i * obj.U.' * diag(arrayfun(@(x) 1 - 2*x,occs)) * conj(obj.V);
        end
        
        %FUNCTION EVALUATE_OPERATOR_EXPECTATION
        %Returns eigenstate expectation of an operator, see
        %EVALUATE_OPERATOR_FROM_CORRELATION_MATRIX for input details.
        function val = evaluate_operator_expectation(obj,occs,A_inds,B_inds)
            val = latticeSC.evaluate_operator_from_correlarion_matrix(...
                obj.correlation_matrix(occs),A_inds,B_inds);
        end
    end
    
    methods(Static)
        %FUNCTION EVALUATE_OPERATOR_FROM_CORRELATION_MATRIX
        %Given an operator described by (A_inds and B_inds), and a
        %two-particle correlation matrix corrmat, returns the expectation
        %of the operator as calculated by Wick's theorem. NOTE: the implied
        %order of operators is that of interleaved A and B operators i.e.
        %\gamma_(A_inds 1)^A \gamma_(B_inds 1)^B \cdots. 
        function val = evaluate_operator_from_correlarion_matrix(corrmat,A_inds,B_inds)
            assert(numel(A_inds) == numel(B_inds),...
                'Unequal number of A and B operators in operator expectation');
            val = - det(corrmat(A_inds,B_inds));
        end
        
        %FUNCTION HAMILTONIAN_FROM_DATA
        %Returns a hamiltonian ready for use in a latticeSC instance given
        %the data in bonds (see constructor for data format). 
        function ham = hamiltonian_from_data(sites,species,bonds,open)
            bonds = latticeSC.verify_bonds(bonds,sites,species);
            
            ham = zeros(sites*species);
            
            for b_index = 1:size(bonds,1)
                bond_indices = bonds{b_index,1};
                diag_index = bond_indices(2) - bond_indices(1);
                if diag_index >= 0
                    diag_position = bond_indices(1);
                else
                    diag_position = bond_indices(2);
                end
                vars = zeros(1,sites*species);
                vars(diag_position:species:(sites*species)) = bonds{b_index,2};
                ham = ham + latticeSC.off_diagonal_matrix(diag_index,vars,1,open);
            end
            
        end
        
        %FUNCTION VERIFY_BONDS
        %Checks that the variable bonds is of the correct form to specify a
        %latticeSC Hamiltonian (see constructor for data format). Also
        %expands any data given as a single number per bond type
        function bonds = verify_bonds(bonds,sites,species)
            assert(iscell(bonds),'Bonds must be provided as cell');
            assert(size(bonds,2) == 2,'Bonds cell must have 2 columns');
            assert(all(cellfun(@(x) all(size(x) == [1,2]),bonds(:,1))),...
                'First column of bonds must contain 1x2 arrays');
            for b_index = 1:size(bonds,1)
                assert(any(bonds{b_index,1} <= species),...
                    'Bonds between two off-site Majoranas should be rewritten');
                if numel(bonds{b_index,2}) == 1
                    bonds{b_index,2} = ones(1,sites).*bonds{b_index,2};
                else
                    assert(numel(bonds{b_index,2}) >= sites,...
                        'Number of bonds of given type less than sites');
                end
            end
        end
        
        %FUNCTION OFF_DIAGONAL_MATRIX
        %Produces a matrix with entries on the Nth diagonal, where N =
        %(off_index). The values of these entries are (vals), which are
        %repeated (cells) times. If open is false, the matrix wraps back
        %around.
        function mat = off_diagonal_matrix(off_index,vals,cells,open)
            if size(vals,1) ~= 1
                error('Values should be given as row matrix');
            end
            if off_index == 0
                mat = diag(repmat(vals,1,cells));
                return;
            end
            cell_size = size(vals,2);
            mat_size = cell_size*cells;
            corner_size = abs(off_index);
            principal_cells = idivide(mat_size - corner_size,uint32(cell_size),'floor');
            resid = mat_size - corner_size - principal_cells*cell_size;
            principal_vec = [repmat(vals,1,principal_cells),vals(1:(resid))];
            
            mat = diag(principal_vec,off_index);
            if ~open
                rem_vec = [vals((resid+1):cell_size),repmat(vals,1,cells - principal_cells - 1)];
                if off_index > 0
                    mat((mat_size+1-corner_size):mat_size,1:corner_size) = diag(rem_vec);
                else
                    mat(1:corner_size,(mat_size+1-corner_size):mat_size) = diag(rem_vec);
                end
            end
        end
        
        %FUNCTION ADD_TO_BONDS
        %Given a cell that specifies the Hamiltonian bonds (bonds), and new
        %data describing the new bonds (bond_index and vals), returns an
        %updated (bonds) cell which should be reassigned to bonds at
        %output.
        function bonds = add_to_bonds(bonds,bond_index,vals)
            if isempty(bonds)
                bonds = {bond_index,vals};
            else
                bonds(end+1,:) = {bond_index,vals};
            end
        end
    end
end

