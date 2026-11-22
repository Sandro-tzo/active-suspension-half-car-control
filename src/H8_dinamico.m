 function [ Ak, Bk, Ck, Dk ] = H8_dinamico(A,B1, B2, C1, D11,D12, C, D)

    nx = size(A)*[1;0];
    nu = size(B1)*[0;1];
    nw = size(B2)*[0;1];
    nz = size(C1)*[1;0];
    np =  size(C)*[1;0];

 % Variabili LMI
 X=sdpvar(nx,nx);
 Y=sdpvar(nx,nx);
 Ack=sdpvar(nx,nx);
 Bck=sdpvar(nx,np);
 Cck=sdpvar(nu,nx);
 Dck=sdpvar(nu,np);
 gamma=sdpvar(1,1);
 
 %vincoli LMI
 V11=[ (A*Y+B1*Cck)+(A*Y+B1*Cck)'  Ack'+(A+B1*Dck*C);
 Ack+(A+B1*Dck*C)'   (X*A+Bck*C)+(X*A+Bck*C)'];

 V12=[(B2+B1*Dck*D)  (C1*Y+D12*Cck)';
      (X*B2+Bck*D)   (C1+D12*Dck*C)'];
 
 V21=[(B2+B1*Dck*D)'  (X*B2+Bck*D)';
     (C1*Y+D12*Cck)   (C1+D12*Dck*C)];

 D1s=(D11+D12*Dck*D);
 V22=[-gamma*eye(nw,nw)  D1s';
      D1s               -gamma*eye(nz,nz)];
 %
 V1=[V11   V12;
     V21   V22]<=0;

 V2=[Y    eye(nx);
     eye(nx)  X]>=0;

 V3=gamma>=0;

 V4=(X>=0);

 V=[V1,V2,V3,V4];
 
 % Ottimizzazione
 opts=sdpsettings;
    opts.solver='lmilab';
 optimize(V, gamma, opts)

% Risultati
 X1=double(X);
 Y1=double(Y);
 Ack1=double(Ack);
 Bck1=double(Bck);
 Cck1=double(Cck);
 Dck1=double(Dck);
 %
 [N, MT]=lu(eye(nx)-X1*Y1);

 % mappatura inversa
 Dk=Dck1;
 Ck=(Cck1- Dk*C*Y1)*inv(MT);
 Bk=inv(N)*(Bck1- X1*B1*Dk);
 Ak=inv(N)*[ Ack1- N*Bk*C*Y1-X1*B1*Ck*(MT)- X1*(A+B1*Dk*C)*Y1 ]*inv(MT);

 end