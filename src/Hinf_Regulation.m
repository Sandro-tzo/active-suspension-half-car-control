function K = Hinf_Regulation( A,B,Bw,C,D,Dw )

    %% Dimensions
    n = size(A)*[1;0];
    m = size(B)*[0;1];
    mw =size(Bw)*[0;1];
    p = size(C)*[1;0];

    gamma = sdpvar(1,1);

    %% Variables
    X = sdpvar(n,n);
    Y = sdpvar(m,n);
    %% LMI
    
    V1= [(A*X+B*Y)+(A*X+B*Y)'        Bw            (C*X+D*Y)';
                Bw'             -gamma*eye(mw)        Dw';
             (C*X+D*Y)               Dw          -gamma*eye(p)]<=0;
    V2 = X>=0;
    V3 = gamma>=0;

    V = [V1,V2,V3];

    %% Solving LMI
    opts=sdpsettings;
    opts.solver='lmilab';
    optimize(V,gamma,opts)
    K = double(Y)/(double(X));

end

