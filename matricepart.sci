function H = MatHilbert(n,m) 
// cette fonction renvoie
// une matrice de Hilbert
[sortie, entree]=argn(0);
if entree==1 then, m=n; end;
H=zeros(n,m);               
for i=1:n               
    for j=1:m 
        H(i,j)=1/(i+j-1);
    end;
end;
endfunction

//objectif est de fariquer une fonction x = Graidient(A,b, rho), s est sol 
// Ax* - b = 0 gradientf(x*) = 0
// gradient a pas constant
// rappel: formulaire
// nb iteration = 100
// algorithme:





function A=MatDiagDom(n,dom)
// renvoie une matrice carrée a diagonale strictement dominante
// dom determine l'importance de la dominance de la diagonale
A=rand(n,n);
[sortie, entree]=argn(0);  
if entree==1 then dom=1, end; // valeur par defaut
A=A-diag(diag(A));
A=A+diag(sum(abs(A),'c'))+dom*eye(A)
endfunction


function A=MatReguliereS(n)
A=triu(rand(n,n));
// on renforce la diagonale de A
A=A+norm(A,'inf')*eye(A);
endfunction

function A=MatSymetrique(n)
A=rand(n,n);
A=A+A';
endfunction


function A=MatHasard(m,n,p) 
[sortie, entree]=argn(0);
select entree
   case 1 then 
      m=n;A=rand(m,n);
   case 2 then
      A=rand(m,n);
   else
      A=rand(m,n);A=p*(2*A-1);     
end;
endfunction


function A=MatReguliere(n)
A=rand(n,n);
// on renforce la diagonale de A
A=A+norm(A,'inf')
endfunction


function A=MatReguliereI(n)
A=tril(rand(n,n));
// on renforce la diagonale de A
A=A+norm(A,'inf')*eye(A);




function A=MatSdp(n)
A=MatSymetrique(n);
[D,P]=bdiag(A);D=abs(D);
D=D+norm(D)*eye(D);
A=P*D*inv(P);
endfunction


function A=MatHasardBin(m,n) 
A=rand(m,n);
A(A<0.5)=0;A(A>=0.5)=1;
endfunction




function A= MatRang(m,n,r)
s=min(m,n);S=max(m,n);
if r>min(m,n) then
   printf('Le rang ne peut etre plus grand que %i ',s) 
   error('Erreur dans la fonction MatRang.') 
else
   A=MatReguliere(s);
   if m>=n then 
         A=[A; ones(S-s,s)]; 
         for k=r+1:s,A(:,k)=rand()*A(:,r);end;
   else
         A=[A  ones(s,S-s)]; 
         for k=r+1:s,A(k,:)=rand()*A(r,:);end;
   end;
end
endfunction


function H = MatHilbert(n,m) 
// cette fonction renvoie
// une matrice de Hilbert
[sortie, entree]=argn(0);
if entree==1 then, m=n; end;
H=zeros(n,m);               
for i=1:n               
    for j=1:m 
        H(i,j)=1/(i+j-1);
    end;
end;
endfunction


function H=House(v)
// Matrice elementaire de Householder
[n,m]=size(v);
if m~=1 then
   error('envoyer un vecteur')
else
   H=eye(n,n);
   n=norm(v);
   if n>1.e-10 then
      H=H -2*v*v'/n/n;
   end;
end;
endfunction

