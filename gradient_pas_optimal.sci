function [x,it]= Gradient_pas_optimal(A,b,rho,tol,itmax,x0)
// A CHANGER
    it = 0;
    x = x0;
    d = b - A*x;
   // while(it < itmax & norm(xn-x) > tol) do
   while(it < itmax) do
        xn = x + rho*d;
        d = b - A*xn;
        rho = (d'*d)/((A*d)'*d);
        it = it + 1;
        if norm(xn-x) < tol then return;end 
        x = xn;
    end;
    mprintf('Arret sur nombre maximum d''itérations %d',it)
endfunction

//norm(d)^2 = d'*d
