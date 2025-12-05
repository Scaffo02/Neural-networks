
% definiamo 4 set di basi:

% Base computazioneale uguale per Q e T
zero = [1 0]';
uno = [0 1]';

%Base R
a = pi/4;
R_p = 1/sqrt(2)*(zero + 1i*exp(1i*a)*uno);
R_m = 1/sqrt(2)*(zero - 1i*exp(1i*a)*uno);

%Base s
S_p = 1/sqrt(2)*(zero + uno);
S_m = 1/sqrt(2)*(zero - uno);

%----------------------------------------------

%Generiamo il data set
n = 100000;
label = zeros(n,1);
C_QS = zeros(n,1);
C_RS = zeros(n,1);
C_RT = zeros(n,1);
C_QT = zeros(n,1);


for i = 1 : n

    
    p = rand();
    psi = 1/sqrt(2)*(kron(zero,zero) + kron(uno,uno));


    rho = p*psi*(psi') + (1-p)*eye(4)/4;
    rho_d_QS = ((kron(zero,S_p)')*rho*kron(zero,S_p)).*kron(zero,S_p)*kron(zero,S_p)' + ((kron(zero,S_m)')*rho*kron(zero,S_m)).*kron(zero,S_m)*kron(zero,S_m)' + ((kron(uno,S_p)')*rho*kron(uno,S_p)).*kron(uno,S_p)*kron(uno,S_p)' + ((kron(uno,S_m)')*rho*kron(uno,S_m)).*kron(uno,S_m)*kron(uno,S_m)' ;                                                                       
    rho_d_RS = ((kron(R_p,S_p)')*rho*kron(R_p,S_p)).*kron(R_p,S_p)*kron(R_p,S_p)' + ((kron(R_p,S_m)')*rho*kron(R_p,S_m)).*kron(R_p,S_m)*kron(R_p,S_m)' + ((kron(R_m,S_p)')*rho*kron(R_m,S_p)).*kron(R_m,S_p)*kron(R_m,S_p)' + ((kron(R_m,S_m)')*rho*kron(R_m,S_m)).*kron(R_m,S_m)*kron(R_m,S_m)' ;       
    rho_d_RT = ((kron(R_p,zero)')*rho*kron(R_p,zero)).*kron(R_p,zero)*kron(R_p,zero)' + ((kron(R_p,uno)')*rho*kron(R_p,uno)).*kron(R_p,uno)*kron(R_p,uno)' + ((kron(R_m,zero)')*rho*kron(R_m,zero)).*kron(R_m,zero)*kron(R_m,zero)' + ((kron(R_m,uno)')*rho*kron(R_m,uno)).*kron(R_m,uno)*kron(R_m,uno)' ;
    rho_d_QT = ((kron(zero,zero)')*rho*kron(zero,zero)).*kron(zero,zero)*kron(zero,zero)' + ((kron(zero,uno)')*rho*kron(zero,uno)).*kron(zero,uno)*kron(zero,uno)' + ((kron(uno,zero)')*rho*kron(uno,zero)).*kron(uno,zero)*kron(uno,zero)' + ((kron(uno,uno)')*rho*kron(uno,uno)).*kron(uno,uno)*kron(uno,uno)' ;
    
    C_QS(i,1) = Entropy(rho_d_QS)-Entropy(rho);
    C_RS(i,1) = Entropy(rho_d_RS)-Entropy(rho);
    C_RT(i,1) = Entropy(rho_d_RT)-Entropy(rho);
    C_QT(i,1) = Entropy(rho_d_QT)-Entropy(rho);

    % valutiamo con il criterio PPT se lo stato è entangled o separabile
    %se lo stato è entangled il risultato sarà 0, se separabile sarà 1. Per
    %essere consistenti con l'articolo poniamo però 1 per gli stati
    %entangled e 0 per quelli seprarabili

    if p > 1/3
    label(i,1) = 1;
    else
    label(i,1) = 0;
    end
    
end


dati = [C_QS, C_RS, C_RT, C_QT, label];
dim = size(dati,1);
dataset = dati(randperm(dim),:);
save('test_set2.mat',"dataset");
disp("Programma terminato")
