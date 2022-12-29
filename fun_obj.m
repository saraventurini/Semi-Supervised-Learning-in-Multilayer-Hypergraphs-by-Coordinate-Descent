% objective function initialization - just p-Laplacian term 

function f = fun_obj(BDX,La,p,K)
    f=0;
    for k=1:K
        f = f + La(k)*sum((vecnorm(BDX{k},p,1)).^p);
    end
end

%{
function f = fun_obj(K,X,La,LpX)
    f=0;
    for k=1:K
        %f = f + La(k)*sum(diag(X*LpX{k})); 
        f = f + La(k)*sum(sum(X.*LpX{k}',2)); 
    end
end
%}


%{
function f = fun_obj(M0,N,La,p,K,x,delta_M)
    f=0;
    for k=1:K
        M0_k = M0{k};
        delta_M_k = delta_M(:,k); 
        for i=1:N
            for j=1:N
                f = f + M0_k(i,j)*(abs((x(i)/delta_M_k(i))-(x(j)/delta_M_k(j))))^p;
                %if M0_k(i,j)*(abs((x(i)/delta_M_k(i))-(x(j)/delta_M_k(j))))^p ~= 0
                %    keyboard
                %end
            end
        end
        %f = La(k)*f^(1/p); 
        f = La(k)*f;
    end
end
%}