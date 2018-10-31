function y=f(x,A,b)
  y = (1/2)* x'*A*x - b'*x
endfunction
 
function y=df(x,A,b)
  y = (A+A')*x/2  - b;
endfunction

function [c,xn,fxn]=gradient(A,b,x0,rho,stop,nmax) 
// une méthode de gradient 
  c=[];
  xn = x0;
  //nmax = 100
  for i=1:nmax 
    xnp1 = xn - rho*df(xn,A,b);
    fxn=f(xn,A,b);
    c=[c,fxn]; //+ add history
    if norm(xnp1-xn) < stop then return;end 
    xn = xnp1;
  end
  mprintf('Arret sur nombre maximum d''itérations %d',nmax) 
endfunction


dom = 5 
n = 10
A = MatDiagDom(n, dom);
A = 1/2*(A + inv(A));

function [x,it]= Gradient_pas_constant(A,b,rho,tol,itmax,x0)

    it = 0;
    x = x0;
    d = b - A*x;
   // while(it < itmax & norm(xn-x) > tol) do
   while(it < itmax) do
        xn = x + rho*d;
        d = b - A*xn;
        it = it + 1;
        if norm(xn-x) < tol then return;end 
        x = xn;
    end;
    mprintf('Arret sur nombre maximum d''itérations %d',it)
endfunction

//rho = 2/(lamda_max+lamda_min) 0.1
//inermax = 1000
B = [2,-1;-1,2]
b = [1;1]
x0 = zeros(b)
//error:Inconsistent row/column dimensions.

// B MARCHE MAIS A MARCHE PAS

x = ones(n,1)
v = -1*ones(n-1,1)
A = diag(n)+diag(v,-1)+diag(v,1)

// grandient converge si et seulement si rho < 2/lamda.max
// rho会在2/lamda.max之间收敛 横坐标是rho
// 作业：理论上找出rho的最优值
// 先去n=5 然后取n=10
// 写出GRANDIENT CONJUGE 和GRADIANT A PAS OPTIMAL
// 比较这三个algo按着这些数据
// 接下来做NON LINERAIRE, QUADRATIQUE : ALGO NEWTON,
// 两页长的rapport, examen之前交

//note
H = MatHilbert(10,10)
Xe = rand(10,1)
b = H* Xe
H*x = b; => x= b/H
err = norm (X - Xe)


