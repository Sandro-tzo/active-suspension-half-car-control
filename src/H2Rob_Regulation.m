function K = H2Rob_Regulation(A1,A2,A3,A4,B1,B2,B3,B4,W1,W2,W3,W4,C,D)
%% Dimensions
    n = size(A1)*[1;0];
    m = size(B1)*[0;1];
    mw = size(W1)*[0;1];
    p = size(C)*[1;0];

    gamma = sdpvar(1,1);

    %% Variables
    X = sdpvar(n,n,'symmetric');
    Y = sdpvar(m,n);
    Q = sdpvar(p,p,'symmetric');

F1=(X>=0);
F2=(gamma>=0);   

F3=([(A1*X+B1*Y)+(A1*X+B1*Y)'   W1;
             W1'           -eye(mw)]<=0);
F4=([(A2*X+B2*Y)+(A2*X+B2*Y)'   W2;
             W2'           -eye(mw)]<=0);
F5=([(A3*X+B3*Y)+(A3*X+B3*Y)'   W3;
             W3'           -eye(mw)]<=0);
F6=([(A4*X+B4*Y)+(A4*X+B4*Y)'   W4;
             W4'           -eye(mw)]<=0);
F7=([    Q    (C*X+D*Y);
    (C*X+D*Y)'    X]>=0);
F8=(trace(Q)<=gamma);

%%unisco i vincoli
F=[F1,F2,F3,F4,F5,F6,F7,F8];
%%risolvo le LMI
opts=sdpsettings;
    opts.solver='lmilab';
    optimize(F,gamma,opts)
    K = double(Y)/(double(X));

end
