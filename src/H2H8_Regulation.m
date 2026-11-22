function K = H2H8_Regulation(A,B,W,C,D,Dw)

%% Dimensions
    n = size(A)*[1;0];
    m = size(B)*[0;1];
    mw =size(W)*[0;1];
    p = size(C)*[1;0];

    gamma = sdpvar(1,1);
    gamma2 = sdpvar(1,1);

    %% Variables
    X = sdpvar(n,n);
    Y = sdpvar(m,n);

a=0.3;
b=0.7;


    Y = sdpvar(m,n);
    Q = sdpvar(p,p);
%%imposto i vincoli LMI
%% Vincoli Hinf

F1=([(A*X+B*Y)+(A*X+B*Y)'        W            (C*X+D*Y)';
                W'             -gamma*eye(mw)        Dw';
             (C*X+D*Y)               Dw      -gamma*eye(p)]<=0);
F2=(X>=0);
F3=(gamma>=0);
%%% Vincoli H2
F4=(gamma2>=0);   
F5=([(A*X+B*Y)+(A*X+B*Y)'   W;
             W'           -eye(mw)]<=0);
F6=([    Q    (C*X+D*Y);
    (C*X+D*Y)'    X]>=0);
F7=(trace(Q)<=gamma2);
F8=(Q>=0);

F=[F1,F2,F3,F4,F5,F6,F7,F8];

%%risolvo le LMI
opts = sdpsettings;
opts.solver = 'lmilab';
optimize(F,a*gamma+b*gamma2,opts)
K=double(Y)/(double(X));


end
