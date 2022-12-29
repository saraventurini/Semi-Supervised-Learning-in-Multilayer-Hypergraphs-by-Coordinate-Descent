% Synthetic datasets SBM %efficiency plots
% Methods: Cyclic Coordinate Descent CCD, Random Coordinate Descent RCD, Gauss-Southwell Coordinate Descent GCD, Gradient Descent GD
% Input : for input and output files' names     dataset_name 
%                                               N number nodes 
%                                               m number communities   
%                                               p in p-Laplacian regularizer
%                                               NL percentage know labels
%                                               runs2 %number of internal runs methods
%                                               qq %number of iterations calculus accuracy and function for efficiency plot
%                                               K number of layers 
%                                               pp = p_in probability arc intercommunity
%                                               rho = p_in/p_out where p_out probability arc intracommunity 
% Output: efficency plots objective function and accuracy in terms of number of flops

function plot_eff_ART(dataset_name,N,m,p,NL,runs2,qq,K,pp,rho)

%iterations x accuracy
acc = readtable('Results/ART/efficiency_plot_ART_acc_'+dataset_name+'.csv');
plot(cat(2,0,(qq:qq:runs2)),acc{:,1},cat(2,0,(qq:qq:runs2)),acc{:,2},cat(2,0,(qq:qq:runs2)),acc{:,3},cat(2,0,(N:N:runs2),runs2),acc{1:length(cat(2,0,(N:N:runs2),runs2)),4},'-o')
title("Efficency plot artificial - N "+num2str(N)+" - m "+num2str(m)+" - p "+num2str(pp)+" - p_{in} "+num2str(p)+" - p_{in}/p_{out} "+num2str(rho)+" - NL "+num2str(NL)+" - K "+num2str(K))
xlabel('Iterations') 
ylabel('Accuracy')
legend('CCD','RCD','GSCD','GD')
plot_name = 'Results/ART/efficiency_plot_ART_acc_'+dataset_name+'.fig'; 
saveas(gca,plot_name)

%iterations x function
fun = readtable('Results/ART/efficiency_plot_ART_fun_'+dataset_name+'.csv');
plot(cat(2,0,(qq:qq:runs2)),fun{:,1},cat(2,0,(qq:qq:runs2)),fun{:,2},cat(2,0,(qq:qq:runs2)),fun{:,3},cat(2,0,(N:N:runs2),runs2),fun{1:length(cat(2,0,(N:N:runs2),runs2)),4},'-o')
title("Efficency plot artificial - N "+num2str(N)+" - m "+num2str(m)+" - p "+num2str(pp)+" - p_{in} "+num2str(p)+" - p_{in}/p_{out} "+num2str(rho)+" - NL "+num2str(NL)+" - K "+num2str(K))
xlabel('Iterations') 
ylabel('Objective function')
legend('CCD','RCD','GSCD','GD')
plot_name = 'Results/ART/efficiency_plot_ART_fun_'+dataset_name+'.fig'; 
saveas(gca,plot_name)

end