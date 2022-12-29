% Synthetic datasets SBM
% Apply efficiency_plot_ART_2 for 
% Parameters: m=4 communities 125 nodes each
%             K=1 number of layers (single layer)
%             p_in = 0.2 probability arc intercommunity
%             p = 2 in p-Laplacian regularizer
%             rho = [3.5, 3, 2.5, 2] where rho =p_in/p_out with p_out probability arc intracommunity 
%             NL = [3,6,9,12] percentage know labels 
% Methods: Cyclic Coordinate Descent CCD, Random Coordinate Descent RCD, Gauss-Southwell Coordinate Descent GCD, Gradient Descent GD
% Output: efficency plots objective function and accuracy in terms of number of flops

function experiments_ART_2
    
    m = 4; %number of communities (if all same size)
    L = 125*ones(1,m); %array dimentions of communities
    p_in = 0.2; 
    K = 1; %number of layers
    p = 2;
    for rho = [3.5, 3, 2.5, 2] %p_in/p_out
        for NL = [3,6,9,12]%percentage of known labels
                fprintf('rho:%f NL:%d\n',rho,NL)
                efficiency_plot_ART_2(L,NL,K,p_in,rho,p);
        end
    end

end

