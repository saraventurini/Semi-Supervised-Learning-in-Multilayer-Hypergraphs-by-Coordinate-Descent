% Real-world datasets (multilayer and hypergraphs)
% Methods: Cyclic Coordinate Descent CCD, Random Coordinate Descent RCD, Gauss-Southwell Coordinate Descent GCD, Gradient Descent GD
% Input: for input and output files' names  dataset_name
%                                           N number of nodes 
%                                           NL percentage know labels
%                                           runs2 %number of internal runs methods
%                                           qq %number of iterations calculus accuracy and function for efficiency plot
% Output: efficency plots objective function and accuracy in terms of number of flops

function plot_eff_REAL(dataset_name,N,NL,runs2,qq)

%iterations x accuracy
acc = readtable('Results/REAL/efficiency_plot_REAL_acc_' + dataset_name + '_' + num2str(NL) + '.csv');
%plot(1:10:num_it,acc{1:N,1},1:10:num_it,acc{1:N,2},1:N:num_it,acc{1:10,3},'-o',1:N:num_it,acc{1:10,4},'-o')
plot(cat(2,0,(qq:qq:runs2)),acc{:,1},cat(2,0,(qq:qq:runs2)),acc{:,2},cat(2,0,(qq:qq:runs2)),acc{:,3},cat(2,0,(N:N:runs2),runs2),acc{1:length(cat(2,0,(N:N:runs2),runs2)),4},'-o')
title("Efficency plot real - "+dataset_name +" - NL "+num2str(NL))
xlabel('Iterations') 
ylabel('Accuracy')
legend('CCD','RCD','GSCD','GD')
plot_name = 'Results/REAL/efficiency_plot_REAL_acc_' + dataset_name + '_' + num2str(NL) +'.fig'; 
saveas(gca,plot_name)

%iterations x function
fun = readtable('Results/REAL/efficiency_plot_REAL_fun_' + dataset_name + '_' + num2str(NL) + '.csv');
plot(cat(2,0,(qq:qq:runs2)),fun{:,1},cat(2,0,(qq:qq:runs2)),fun{:,2},cat(2,0,(qq:qq:runs2)),fun{:,3},cat(2,0,(N:N:runs2),runs2),fun{1:length(cat(2,0,(N:N:runs2),runs2)),4},'-o')
title("Efficency plot real - "+dataset_name +" - NL "+num2str(NL))
xlabel('Iterations') 
ylabel('Objective function')
legend('CCD','RCD','GSCD','GD')
plot_name = 'Results/REAL/efficiency_plot_REAL_fun_' + dataset_name + '_' + num2str(NL) +'.fig'; 
saveas(gca,plot_name)

end