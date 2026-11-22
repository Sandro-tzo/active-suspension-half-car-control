function K = H2L1_Regulation(A,B1,B2,C,D31,D32,lambda)

% definiamo la LMI per la risoluzione del problema H2L1

%% Dimensions
    n = size(A)*[1;0];
    m = size(B1)*[0;1];
    mw = size(B2)*[0;1];
    p = size(C)*[1;0];

    yalmip('clear');


    %% Variables
    gamma1 = sdpvar(1,1);
    X = sdpvar(n,n);
    Y = sdpvar(m,n);
    mu=sdpvar(1,1);
    
    gamma2 = sdpvar(1,1);
    Q=sdpvar(p,p);
 
a=0.3;
b=0.7;

%%imposto i vincoli LMI L1
F1=(X>=0);
F2=(gamma1>=0);   
F3=([(A*X+B1*Y)+(A*X+B1*Y)'+lambda*X   B2;
             B2'           -mu*eye(mw)]<=0);
F4=([ lambda*X      zeros(n,mw)   (C*X+D32*Y)';
       zeros(mw,n)  (gamma1-mu)*eye(mw)   D31';
    (C*X+D32*Y)    D31     gamma1*eye(p)]>=0);
F5=(mu>=0);

%%imposto i vincoli LMI H2
F6=(gamma2>=0);   
F7=([(A*X+B1*Y)+(A*X+B1*Y)'   B2;
             B2'           -eye(mw)]<=0);
F8=([ Q      (C*X+D32*Y);
    (C*X+D32*Y)'    X]>=0);
F9=(trace(Q)<=gamma2);
F10=(Q>=0);
%%unisco i vincoli
F=F1+F2+F3+F4+F5+F6+F7+F8+F9+F10;
%%risolvo le LMI
%% Solving LMI
    opts=sdpsettings;
    opts.solver='lmilab';
    optimize(F,a*gamma1+b*gamma2,opts)
    K = double(Y)/(double(X));

end