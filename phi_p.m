%function used to define the p-Laplacian 

function f = phi_p(Y,p)
    f =  abs(Y).^(p-1) .* sign(Y);
end