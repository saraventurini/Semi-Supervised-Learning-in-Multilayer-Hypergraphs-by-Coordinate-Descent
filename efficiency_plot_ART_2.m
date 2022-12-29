% Synthetic datasets SBM
% Methods: Cyclic Coordinate Descent CCD, Random Coordinate Descent RCD, Gauss-Southwell Coordinate Descent GCD, Gradient Descent GD
% Input: L array where L(c) is size of community c 
%        (m*sum(L) matrix size)
%        NL percentage know labels 
%        K number of layers 
%        p_in probability arc intercommunity
%        rho = p_in/p_out where p_out probability arc intracommunity 
%        p in p-Laplacian regularizer
% Output: efficency plots objective function and accuracy in terms of number of flops (csv files and plots)

function efficiency_plot_ART_2(L,NL,K,p_in,rho,p) 

q = p_in/rho;

%fix random sead
rng('default');

N = sum(L); %number of nodes
m = length(L); %number of communities
runs1 = 5; %5 random matrices  
molt = 10; 
runs2 = molt*N-1;%10*N-1;%20*N-1; %number of internal runs methods %30*N-1 
IT = (runs2+1)/N - 1; %IT = 10-1;%20-1; %iterations GD  
qq_fun = 1; %fix(N/6); %number of iterations calculus accuracy and function for efficiency plot
runs2 = fix(runs2/qq_fun)*qq_fun;

labels = []; %labels
for com=1:m
   labels = cat(2,labels,com*ones(1,L(com)));
end

%M_COM cell communities indeces
%SC size communities 
M_COM=cell(1,m);
SC=zeros(1,m);
for com=1:m
    M_COM{com}=find(labels==com);
    %SC(com)=sum(labels'==com);
    SC(com)=length(M_COM{com});
end 

c = fix((SC/100)*NL); %number of known labels in each community
c(c==0 | c==1)=2; %at least 2 node labeled each community
c(SC==1)=1; %unless community just one node 

%M_COM_NL cell indices know labels each community 
M_COM_NL=cell(1,m);
for com=1:m
    com_R = M_COM{com}(randperm(SC(com))); 
    M_COM_NL{com}=com_R(1:c(com));
end

%Y matrix
Y = zeros(m,N);
for com=1:m
    Y(com,M_COM_NL{com})=1;
end
Y = sparse(Y);
%normalized
Y = Y./vecnorm(Y,2,2); %vecnorm(A,p,dim)

%expected communities (without known labels) 
RP = cell2mat(M_COM_NL); 
c_exp=labels;
c_exp(RP)=[];
%indices of labeled nodes
I=RP;

%final matrix: matrix_acc(i,j,k) = accuracy of matrix j at iteration i with method k 
matrix_acc=zeros(runs2+1,runs1,4); %accuracy
matrix_fun=zeros(runs2+1,runs1,4); %function

for r=1:runs1 
    fprintf('samples labels:%d\n',r)

    %connected graph
    M = cell(1,K);
    for k=1:K
        M{k}=0;
        while ~isempty(find(sum(M{k})==0, 1))
            M{k} = sparse(adjacent_matrix_generator(L,p_in,q));
        end
    end

    %%delta_M(i,l) = root of the degree of node i on layer l respect to graph M 
    delta_M=zeros(N,K);
    %matrices to calculate the function 
    M_bar = cell(1,K);
    for k=1:K  
        delta_M(:,k) = sqrt(sum(M{k},2));
        delta_M_2 = 1./delta_M(:,k); 
        M_bar{k} = delta_M_2'.*M{k}.*delta_M_2;
    end

    %memorize matrix B*D^(-1/2) for iterative calculation function, gradient, hessian
    size_M = size(M{1},1); 
    B = cell(1,K);
    D = cell(1,K);
    BD = cell(1,K);
    DB = cell(1,K);
    for k=1:K
        %degree
        D_ = sparse(1:size_M,1:size_M,1./delta_M(:,k),size_M,size_M);
        D{k}=D_;
        G = graph(M{k},'omitselfloops'); 
        %incidente matrix
        B_ = -incidence(G)';
        B{k}=B_;
        BD{k} = B_*D_; 
        DB{k} = BD{k}';
    end

    %parameters
    %fixed
    val = 1/(K*2); %sum all alpha has to be 1 
    Alpha_CCD = [val val];
    La_CCD=repmat(Alpha_CCD(1)/(K*Alpha_CCD(2)),K,1);
    Alpha_RCD = [val val];
    La_RCD=repmat(Alpha_RCD(1)/(K*Alpha_RCD(2)),K,1);
    Alpha_GSCD = [val val];
    La_GSCD=repmat(Alpha_GSCD(1)/(K*Alpha_GSCD(2)),K,1);
    Alpha_GD = [val val];
    La_GD=repmat(Alpha_GD(1)/(K*Alpha_GD(2)),K,1);

    %step for p=2 %use Q 
    Q_CCD = sparse(N,N); %I need just the diagonal 
    for k=1:K
        Q_CCD = Q_CCD + 2*(La_CCD(k))*(speye(N) - M_bar{k}) ;
    end
    Q_CCD = 2*speye(N) + Q_CCD; 

    Q_RCD = sparse(N,N);
    for k=1:K
        Q_RCD = Q_RCD + 2*(La_RCD(k))*(speye(N) - M_bar{k}) ;
    end
    Q_RCD = 2*speye(N) + Q_RCD; 

    Q_GSCD = sparse(N,N);
    for k=1:K
        Q_GSCD = Q_GSCD + 2*(La_GSCD(k))*(speye(N) - M_bar{k}) ;
    end
    Q_GSCD = 2*speye(N) + Q_GSCD; 

    Q_GD = sparse(N,N);
    for k=1:K
        Q_GD = Q_GD + 2*(La_GD(k))*(speye(N) - M_bar{k}) ;
    end
    Q_GD = 2*speye(N) + Q_GD; 
    Q_GD_eig_max = eigs(Q_GD,1); %MAX EIGENVALUE Q_GD
    step_GD = 1/Q_GD_eig_max;

    %Accuracy at the beginning
    [~,index] = max(Y,[],1,'linear'); 
    [comm,~] = ind2sub([m,N],index);
    comm(I)=[]; %not consider labeled nodes
    matrix_acc(1,r,1) = ((N-sum(c))-wrong(c_exp,comm))/(N-sum(c));
    matrix_acc(1,r,2) = matrix_acc(1,r,1);
    matrix_acc(1,r,3) = matrix_acc(1,r,1);
    matrix_acc(1,r,4) = matrix_acc(1,r,1);

    %Variables at the beginning
    matrix_z_CCD=sparse(zeros(m,N)); 
    matrix_z_RCD=sparse(zeros(m,N)); 
    matrix_z_GSCD=sparse(zeros(m,N)); 
    matrix_z_GD=sparse(zeros(m,N)); 

    %save for computational cost - incremental computations
    %CCD
    BDX_CCD = cell(1,K);
    for k=1:K
        BDX_CCD{k} = BD{k}*matrix_z_CCD'; 
    end
    %RCD
    BDX_RCD = cell(1,K);
    for k=1:K
        BDX_RCD{k} = BD{k}*matrix_z_RCD'; 
    end
    %GSCD
    BDX_GSCD = cell(1,K);
    for k=1:K
        BDX_GSCD{k} = BD{k}*matrix_z_GSCD'; 
    end
    %CCD
    BDX_GD = cell(1,K);
    for k=1:K
        BDX_GD{k} = BD{k}*matrix_z_GD'; 
    end 

    %function at the beginning
    f_CCD_2 = sum(vecnorm(matrix_z_CCD - Y,2,2).^2) +  fun_obj(BDX_CCD,La_CCD,p,K); 
    matrix_fun(1,r,1) = f_CCD_2; 
    f_RCD_2 = sum(vecnorm(matrix_z_RCD - Y,2,2).^2) + fun_obj(BDX_RCD,La_RCD,p,K); 
    matrix_fun(1,r,2) = f_RCD_2;
    f_GSCD_2 = sum(vecnorm(matrix_z_GSCD - Y,2,2).^2) + fun_obj(BDX_GSCD,La_GSCD,p,K); 
    matrix_fun(1,r,3) = f_GSCD_2;
    matrix_fun(1,r,4) = sum(vecnorm(matrix_z_GD - Y,2,2).^2) + fun_obj(BDX_GD,La_GD,p,K);
    
    %gradient 
    %CCD
    grad_CCD_1 = 2*(matrix_z_CCD-Y); 
    LpX_CCD = cell(1,K); %p-Laplacian
    grad_CCD_2 = sparse(m,N); 
    for k=1:K
        LpX_CCD{k} = DB{k}*phi_p(BDX_CCD{k},p); 
        grad_CCD_2 = grad_CCD_2 + La_CCD(k)*p*(LpX_CCD{k})'; 
    end 
    grad_CCD =  grad_CCD_1 + grad_CCD_2; %divide gradient in two parts for incremental computation
                                         %grad_CCD_1 data fitting %grad_CCD_2 p-Laplacian
    %RCD
    grad_RCD_1 = 2*(matrix_z_RCD-Y);
    LpX_RCD = cell(1,K); %p-Laplacian
    grad_RCD_2 = sparse(m,N); 
    for k=1:K
        LpX_RCD{k} = DB{k}*phi_p(BDX_RCD{k},p); 
        grad_RCD_2 = grad_RCD_2 + La_RCD(k)*p*(LpX_RCD{k})'; 
    end 
    grad_RCD =  grad_RCD_1 + grad_RCD_2; %divide gradient in two parts for incremental computation
                                         %grad_CGS_1 data fitting %grad_CGS_2 p-Laplacian
    %GSCD
    grad_GSCD_1 = 2*(matrix_z_GSCD-Y);
    LpX_GSCD = cell(1,K); %p-Laplacian
    grad_GSCD_2 = sparse(m,N); 
    for k=1:K
        LpX_GSCD{k} = DB{k}*phi_p(BDX_GSCD{k},p); 
        grad_GSCD_2 = grad_GSCD_2 + La_GSCD(k)*p*(LpX_GSCD{k})'; 
    end 
    grad_GSCD =  grad_GSCD_1 + grad_GSCD_2; %divide gradient in two parts for incremental computation
                                         %grad_CDAS_1 data fitting %grad_CDAS_2 p-Laplacian
    %GD
    LpX_GD = cell(1,K); %p-Laplacian
    grad_GD = 2*(matrix_z_GSCD-Y)';
    for k=1:K
        LpX_GD{k} = DB{k}*phi_p(BDX_GD{k},p); 
        grad_GD = grad_GD + La_GD(k)*p*LpX_GD{k}; 
    end

    for rr=0:(runs2-1)

        % CCD METHOD
        %CDD each epoch lenght N permute the indices
        mod_rrN = mod(rr,N);
        if mod_rrN==0
           perm = randperm(N);       
        end 

        j_t = perm(mod_rrN+1); %Coordinate selection 
        matrix_z_CCD_j_t = matrix_z_CCD(:,j_t); 
        grad_CCD_j_t = grad_CCD(:,j_t);  
        step_CCD= 1/(2*Q_CCD(j_t,j_t));
        %minimization
        zz = matrix_z_CCD_j_t - step_CCD.*grad_CCD_j_t;
        %updates
        for k=1:K
            BDX_CCD{k} = BDX_CCD{k} + repmat(BD{k}(:,j_t),1,m).*(zz' - matrix_z_CCD_j_t') ; 
        end  
        grad_CCD_1(:,j_t) =  grad_CCD_1(:,j_t) + 2*(zz-matrix_z_CCD_j_t);
        grad_CCD_2 = sparse(m,N); 
        for k=1:K 
            grad_CCD_2 = grad_CCD_2 + La_CCD(k)*p*(DB{k}*BDX_CCD{k})'; %cause p=2
        end 
        grad_CCD = grad_CCD_1 + grad_CCD_2;  

        matrix_z_CCD(:,j_t) = zz;
        matrix_fun(rr+2,r,1) = sum(vecnorm(matrix_z_CCD - Y,2,2).^2) + fun_obj(BDX_CCD,La_CCD,p,K);
        

        % RCD METHOD
        j_t = randi(N); %Coordinate selection
        matrix_z_RCD_j_t = matrix_z_RCD(:,j_t); 
        grad_RCD_j_t = grad_RCD(:,j_t);  
        step_RCD= 1/(2*Q_RCD(j_t,j_t));
        tic;
        %minimization
        zz = matrix_z_RCD_j_t - step_RCD.*grad_RCD_j_t;
        %updates
        for k=1:K
            BDX_RCD{k} = BDX_RCD{k} + repmat(BD{k}(:,j_t),1,m).*(zz' - matrix_z_RCD_j_t') ; 
        end
        grad_RCD_1(:,j_t) =  grad_RCD_1(:,j_t) + 2*(zz-matrix_z_RCD_j_t);
        grad_RCD_2 = sparse(m,N); 
        for k=1:K 
            grad_RCD_2 = grad_RCD_2 + La_RCD(k)*p*(DB{k}*BDX_RCD{k})'; 
        end 
        grad_RCD = grad_RCD_1 + grad_RCD_2;   
        matrix_z_RCD(:,j_t) = zz;
        matrix_fun(rr+2,r,2) = sum(vecnorm(matrix_z_RCD - Y,2,2).^2) + fun_obj(BDX_RCD,La_RCD,p,K);

        % GSCD METHOD
        [~,j_t] = max(abs(grad_GSCD),[],2); %Coordinate selection
        matrix_z_GSCD_j_t = diag(matrix_z_GSCD(:,j_t));
        grad_GSCD_j_t = diag(grad_GSCD(:,j_t));
        step_GSCD= 1./(2*diag(Q_GSCD(j_t,j_t)));
        %minimization
        zz = matrix_z_GSCD_j_t - step_GSCD.*grad_GSCD_j_t;  
        %updates
        for k=1:K
            BDX_GSCD{k} = BDX_GSCD{k} + BD{k}(:,j_t).*(zz' - matrix_z_GSCD_j_t') ; 
        end
        lInd = sub2ind(size(grad_GSCD_1),1:m,j_t');
        grad_GSCD_1(lInd) =  grad_GSCD_1(lInd) + 2*(zz-matrix_z_GSCD_j_t)';
        grad_GSCD_2 = sparse(m,N); 
        for k=1:K
            grad_GSCD_2 = grad_GSCD_2 + La_GSCD(k)*p*(DB{k}*BDX_GSCD{k})'; 
        end 
        grad_GSCD = grad_GSCD_1 + grad_GSCD_2; 
        lInd = sub2ind(size(matrix_z_GSCD),1:m,j_t');
        matrix_z_GSCD(lInd) = zz;
        matrix_fun(rr+2,r,3) = sum(vecnorm(matrix_z_GSCD - Y,2,2).^2) + fun_obj(BDX_GSCD,La_GSCD,p,K);

        % GD METHOD
        if rr<=IT %less iterations 
            %minimization   
            zz = matrix_z_GD - (step_GD.*grad_GD)';
            %updates
            BDX_GD = cell(1,K);
            for k=1:K
                BDX_GD{k} = BD{k}*zz'; 
            end
            %gradient
            LpX_GD = cell(1,K); %p-Laplacian 
            grad_GD = 2*(zz-Y)'; 
            for k=1:K
                LpX_GD{k} = DB{k}*BDX_GD{k}; 
                grad_GD = grad_GD + La_GD(k)*p*LpX_GD{k}; 
            end
            %function
            matrix_fun(rr+2,r,4) = sum(vecnorm(zz - Y,2,2).^2) + fun_obj(BDX_GD,La_GD,p,K); 
            matrix_z_GD = zz;
        end
        
        % CCD
        %Communities partition at this iteration
        [~,index] = max(matrix_z_CCD,[],1,'linear'); 
        [comm,~] = ind2sub([m,N],index);
        %Accuracy
        comm(I)=[]; %not consider labeled nodes
        matrix_acc(rr+2,r,1) = ((N-sum(c))-wrong(c_exp,comm))/(N-sum(c));

        % RCD
        [~,index] = max(matrix_z_RCD,[],1,'linear'); 
        [comm,~] = ind2sub([m,N],index);
        %Accuracy
        comm(I)=[]; %not consider labeled nodes
        matrix_acc(rr+2,r,2) = ((N-sum(c))-wrong(c_exp,comm))/(N-sum(c));
    
        % GSCD
        %Communities partition at this iteration;
        [~,index] = max(matrix_z_GSCD,[],1,'linear'); 
        [comm,~] = ind2sub([m,N],index);
        %Accuracy
        comm(I)=[]; %not consider labeled nodes
        matrix_acc(rr+2,r,3) = ((N-sum(c))-wrong(c_exp,comm))/(N-sum(c));

        if rr<=IT 
            % GD
            %Communities partition at this iteration
            [~,index] = max(matrix_z_GD,[],1,'linear'); 
            [comm,~] = ind2sub([m,N],index);
            %Accuracy
            comm(I)=[]; %not consider labeled nodes
            matrix_acc(rr+2,r,4) = ((N-sum(c))-wrong(c_exp,comm))/(N-sum(c));
        end  

    end
end

%Average Accuracy
mean_matrix_acc=full(mean(matrix_acc,2));
%Average Function
mean_matrix_fun=full(mean(matrix_fun,2));

%Print on file 
dataset_name = num2str(N) + "_" + num2str(m) + "_" + num2str(p) + "_" + num2str(p_in) + "_" + num2str(rho) + "_" + num2str(NL) + "_" + num2str(K);
file_name1 = fullfile('C:\Users\Sara\Desktop\Codes - Copia\Results\ART');
file_name = "efficiency_plot_ART_acc_" + dataset_name + ".csv";
file_name2 = fullfile (file_name1,file_name);
T = array2table(cat(2,mean_matrix_acc(:,:,1),mean_matrix_acc(:,:,2),mean_matrix_acc(:,:,3),mean_matrix_acc(:,:,4))); %from matrix to table
writetable(T,file_name2) %csv file
file_name = "efficiency_plot_ART_fun_" + dataset_name + ".csv";
file_name2 = fullfile (file_name1,file_name);
T = array2table(cat(2,mean_matrix_fun(:,:,1),mean_matrix_fun(:,:,2),mean_matrix_fun(:,:,3),mean_matrix_fun(:,:,4))); %from matrix to table
writetable(T,file_name2) %csv file
%plot
plot_eff_ART(dataset_name,N,m,p_in,NL,runs2,qq_fun,K,p,rho)

end