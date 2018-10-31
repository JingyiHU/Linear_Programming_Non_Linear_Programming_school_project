
// Exercice 1.1 r -> V(r) *************************************************
r = 0.8:0.03:3;
plot(r, r.^(-12)-2*r.^(-6),'-')

// Exeercice 1.2 **********************************************************
// fonction de calcul de J
function J = lennardjones(x)
    N = length(x)/3;
    X = matrix(x,3,N);
    J = 0;
    r = 0;
    for i = 1:N
        for j = i+1:N
           r = norm(X(:,j)-X(:,i));
           V = r.^(-12)-2*r.^(-6);
           J = J + V;
        end
    end
endfunction

// Exemple de Test de fonction "lennardjones"

// N = 2
test2 = [1;0;0;0;0;0]; 
J_test2 = lennardjones(test2);
// Résultat = -1 ( = Jopt | N = 2)

// N = 3 
test3 = [1;0;0;0;1;0;0;0;1];
J_test3 = lennardjones(test3);
// Résultat = - 0.703125

// Exercice 1.3 ********************************************************** 

//focntion d'optimisation "CostLennardjones"
function [f,g,ind] = CostLennardjones(x, ind)
    f = lennardjones(x);
    g = numderivative(lennardjones,x)
endfunction

// Test de la fonction "CostLennardjones" de N = 3 
x = [1;0;0;0;1;0;0;0;1];
[fopt, xopt]=optim(CostLennardjones,x,"qn")

//xopt  = 
//    0.8047358  
//    0.0976273  
//    0.0976273  
//    0.0975944  
//    0.8046994  
//    0.0975926  
//    0.0975944  
//    0.0975926  
//    0.8046994  
// fopt  = - 3. 

// Exercice 1.4  help optim **********************************************

// Exercice 1.5 ********************************************************** 

// N = 4
x = [0;0;0;1;0;0;0;1;0;0;0;1];
[fopt, xopt]=optim(CostLennardjones,x,"qn")

// xopt  =
//  - 0.1039383  
//  - 0.1042170  
//  - 0.1025922  
//    0.8391068  
//    0.1309040  
//    0.1327461  
//    0.1330052  
//    0.8390099  
//    0.1301751  
//    0.1317074  
//    0.1341823  
//    0.8395525  
// fopt  = - 6.  
 
// Exercice 1.7 **********************************************************

// N = 8
x = [0;0;0;1;0;0;1;1;0;0;1;0;0;0;1;1;0;1;1;1;1;0;1;1];
[fopt, xopt]=optim(CostLennardjones,x,"qn")
[fopt, xopt]=optim(CostLennardjones,xopt,"qn")  

// Mêmes résultats pour ces deux tentatives : 
// xopt  =
//  - 0.1842307  
//  - 0.0676354  
//    0.0104576  
//    0.7873812  
//    0.0951828  
//    0.1396228  
//    1.0698784  
//    0.9943452  
//  - 0.1749013  
//    0.1424512  
//    0.8594965  
//    0.1551252  
//    0.1532971  
//    0.1990052  
//    0.9061141  
//    1.0878121  
//  - 0.0977627  
//    1.0668426  
//    0.9117640  
//    0.8412758  
//    0.7940118  
//    0.0214335  
//    1.1660131  
//    1.0924753  
// fopt  = - 18.976056  

// Exercice 1.8 **********************************************************
  //[fopt, xopt]=optim(CostLennardjones,xopt,"qn")
  // Même résultat => peut étre un min global, soit la fonction coinsée dans un min local => re éssayer pls points initicaux  
 x = [0;0.5;0;1;0;0;1;1;0;0;1;0;0;0;1;1;0;1;1;1;1;0;1;1];
 [fopt, xopt]=optim(CostLennardjones,x,"qn");
 // Loop lourd !!!!!!!!
 
// Exercice 1.9 ********************************************************** 

//Fonction d'optimisation avec pénalisation 

 function J = lennardjonesPenalise(x)
    N = length(x)/3;
    X = matrix(x,3,N);
    J = 0;
    r = 0;
    lambda = 10;  // On fixe lambda ici dans au sein de la fonction 
    for i = 1:N
        for j = i+1:N
           r = norm(X(:,j)-X(:,i));
           V = r.^(-12)-2*r.^(-6);
           J = J + V;
        end
    end
    J = J + lambda*(norm(X(:,1)-[0;0;0]))^2
endfunction

function [f,g,ind] = CostLennardjonesPenalise(x, ind)
    f = lennardjonesPenalise(x);
    g = numderivative(lennardjonesPenalise,x)
endfunction

// Exercice 1.5 **********************************************************

// Refait de la Question 5 : N = 4
x = [0;0;0;1;0;0;0;1;0;0;0;1];
[fopt, xopt]=optim(CostLennardjonesPenalise,x,"qn")
// Résultat : xopt[1] assez proche de 0, J reste identique = - 6

// N = 3 (X1 de départ != 0)
x = [1;0;0;0;1;0;0;0;1];
[fopt, xopt]=optim(CostLennardjonesPenalise,x,"qn") 
// Résultat : xopt[1] assez proche de 0, J reste identique = - 3


// N = 8
x = [0;0;0;1;0;0;1;1;0;0;1;0;0;0;1;1;0;1;1;1;1;0;1;1];
[fopt, xopt]=optim(CostLennardjonesPenalise,x,"qn") 

