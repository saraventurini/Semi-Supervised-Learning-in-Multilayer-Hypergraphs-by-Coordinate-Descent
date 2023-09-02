% From files downloaded from Austin Benson page https://www.cs.cornell.edu/~arb/data/
% dataset: trivago-clicks,walmart-trips,amazon-reviews,contact-high-school,contact-primary-school,house-bills,house-committees,senate-bills,senate-committees
% no dataset: stackoverflow-answers and mathoverflow-answers because overlapping communities 
% INPUT: dataset name
% OUTPUT: print number communities, number nodes, number hyperedges, mean and median and max size hyperedges
%         save labels, incident matrix, adjacency matrix

%AB_hypergraphs("trivago-clicks")
%AB_hypergraphs("walmart-trips")
%AB_hypergraphs("amazon-reviews")
%AB_hypergraphs("contact-high-school")
%AB_hypergraphs("contact-primary-school")
%AB_hypergraphs("house-bills")
%AB_hypergraphs("house-committees")
%AB_hypergraphs("senate-bills")
%AB_hypergraphs("senate-committees")
function AB_hypergraphs(dataset_name)

%from file to matrix
labels = readmatrix(dataset_name+"/node-labels-"+dataset_name+".txt")'; 
number_communties = length(unique(labels)); %number of communities
N = size(labels,2); %number of nodes
fprintf('%s %d \n', 'Number of communities: ',number_communties);
fprintf('%s %d \n', 'Number of nodes: ',N);

fid = fopen(dataset_name+"/hyperedges-"+dataset_name+".txt");

%number of line/hyperedges
g = textscan(fid,'%s','delimiter','\n');
number_hyperedges = length(g{1});
fprintf('%s %d \n', 'Number of hyperedges: ',number_hyperedges);

%create incident matrix 
I = sparse(N,number_hyperedges); %incident matrix hypergraph 
fid = fopen(dataset_name+"/hyperedges-"+dataset_name+".txt");
tline = fgetl(fid);  %each row: nodes in one hyperedge
i = 0;
while ischar(tline) %unless last row 
    %disp(tline)
    %newStr = split(tline,','); %split nodes in the hyperedge
    nodes_hyperedge = str2num(tline);
    I(nodes_hyperedge,i+1) = 1;
    tline = fgetl(fid); %newt line
    i = i+1;
end
fclose(fid);

%mean hyperedge size
size_hyperedges = sum(I,1); 
mean_hyperedge_size = sum(size_hyperedges)/number_hyperedges;
fprintf('%s %.4f \n', 'Mean hyperedge size: ',full(mean_hyperedge_size));
%median hyperedge size
size_hyperedges_order = sort(size_hyperedges);
if  rem(number_hyperedges,2)==0 %even number of hyperedges 
    median_hyperedge_size = (size_hyperedges_order(number_hyperedges/2)+size_hyperedges_order(number_hyperedges/2+1))/2;
else %odd number of edges 
    median_hyperedge_size = size_hyperedges_order(fix(number_hyperedges/2)+1);
end
fprintf('%s %.4f \n', 'Median hyperedge size: ',full(median_hyperedge_size));
%rank of hypergraph (maximum hyperedge size)
max_hyperedge_size = size_hyperedges_order(end); %max(size_hyperedges);
fprintf('%s %d \n', 'Maximum hyperedge size: ',full(max_hyperedge_size));

%adjacency matrix A = II' - D;
A = I*I'; 
A = A - diag(diag(A));
fprintf('%s %d \n', 'Density adjacency matrix: ',nnz(A)/(N*N-N)); %density A matrix


%save matrices in dataset folder
%current_folder = cd;
dataset_folder = fullfile('\Data\', dataset_name);
save(fullfile(dataset_folder,'labels_'+dataset_name+'.mat'),'labels')
save(fullfile(dataset_folder,'I_'+dataset_name+'.mat'),'I')
save(fullfile(dataset_folder,'A_'+dataset_name+'.mat'),'A')
end