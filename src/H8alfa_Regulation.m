function K = H8alfa_Regulation(A,B,W,C,D,Dw)

alfa=57;
teta=pi/4;
r=110;

%% Dimensions
    n = size(A)*[1;0];
    m = size(B)*[0;1];
    mw =size(W)*[0;1];
    p = size(C)*[1;0];

    gamma = sdpvar(1,1);

    %% Variables
    X = sdpvar(n,n);
    Y = sdpvar(m,n);
    Q = sdpvar(p,p);
%%imposto i vincoli x Regioni LMI
S1=([2*alfa*X+A*X+X*A'+B*Y+Y'*B']<=0);
S2=([sin(teta)*(A*X+B*Y+X*A'+Y'*B') cos(teta)*(-A*X-B*Y+X*A'+Y'*B');
    cos(teta)*(A*X+B*Y-X*A'-Y'*B') sin(teta)*(A*X+B*Y+X*A'+Y'*B')]<=0);
S3=([-r*X  A;
      A' -r*X]<=0);
%
% Vincoli Hinfinito
F1=([(A*X+B*Y)+(A*X+B*Y)'        W            (C*X+D*Y)';
                W'             -gamma*eye(mw)        Dw';
             (C*X+D*Y)               Dw      -gamma*eye(p)]<=0);
F2=(X>=0);
F3=(gamma>=0);

F = [F1, F2, F3, S1, S2, S3];

opts = sdpsettings;
    opts.solver = 'lmilab';
    optimize(F,gamma,opts)
K=double(Y)/(double(X));

end


