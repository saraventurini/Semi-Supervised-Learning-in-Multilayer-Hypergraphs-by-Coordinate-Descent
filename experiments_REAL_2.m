% Real-world datasets (multilayer and hypergraphs)
% Apply efficiency_plot_REAL_2 for 
% Parameters: A cell, A{k} adjacency matrix layer k
%             labels ground truth communities 
%             NL = [3,6,9,12] percentage know labels (NL = [15,18,21,24] Wikipedia)         
%             dataset_name for input and output files' names 
%             p = 2 in p-Laplacian regularizer    
% Methods: Cyclic Coordinate Descent CCD, Random Coordinate Descent RCD, Gauss-Southwell Coordinate Descent GCD, Gradient Descent GD
% Output: efficency plots objective function and accuracy in terms of number of flops

function experiments_REAL_2

p = 2;

%multi-layer graph %Tudisco 
dataset_name_list = ["3sources","BBCSport2view_544","UCI_mfeat","cora"];
%dataset_name_list = ["WikipediaArticles"];
%single layer hypergraphs %Austin Benson
%dataset_name_list = ["contact-high-school","contact-primary-school"];
for dataset_name = dataset_name_list

   %multi-layer graph %Tudisco 
   load('Real_datasets/'+dataset_name,'W_cell','labels') 
   %single layer hypergraphs %Austin Benson
   %load('Data/Single_layer/'+dataset_name+'/A_'+dataset_name+'.mat','A') %adjacency matrix hypergraph 
   %load('Data/Single_layer/'+dataset_name+'/labels_'+dataset_name+'.mat','labels') %labels 

   for NL = [3,6,9,12]%percentage of known labels
   %for NL = [15,18,21,24]%Wikipedia
        fprintf('dataset:%s NL:%d\n',dataset_name,NL)
        %multi-layer graph %Tudisco 
        efficiency_plot_REAL_2(W_cell,labels,NL,dataset_name,p);
        %single layer hypergraphs %Austin Benson
        %efficiency_plot_REAL_2(A,labels,NL,dataset_name,p);
   end
end
