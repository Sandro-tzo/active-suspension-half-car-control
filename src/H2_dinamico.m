 function [ Ak, Bk, Ck, Dk ] = H2_dinamico(A,B1, B2, C2,D21, D22, C, D)

    nx = size(A)*[1;0];
    nu = size(B1)*[0;1];
    nw = size(B2)*[0;1];
    nz = size(C2)*[1;0];
    np =  size(C)*[1;0];

 % Variabili LMI
 X=sdpvar(nx,nx);
 Y=sdpvar(nx,nx);
 Q=sdpvar(nz,nz);
 Ack=sdpvar(nx,nx);
 Bck=sdpvar(nx,np);
 Cck=sdpvar(nu,nx);
 Dck=sdpvar(nu,np);
 gamma=sdpvar(1,1);
 
 %vincoli LMI

V21=(Ack+(A+B1*Dck*C)');
V31=(B2+B1*Dck*D)';
V32=(X*B2+Bck*D)';

 %
 V1=[(A*Y+B1*Cck)+(A*Y+B1*Cck)', V21', V31';
     V21,   (X*A+Bck*C)+(X*A+Bck*C)', V32';
     V31,   V32     -eye(nw, nw)]<=0;

 V2=[Q ,  (C2*Y+D22*Cck), (C2+D22*Dck*C);
     (C2*Y+D22*Cck)',  Y,    eye(nx) ;
     (C2+D22*Dck*C)',   eye(nx),  X]>=0;

 V3=gamma>=0;

 V4=(trace(Q)<=gamma);

 V5=(X>=0);

 V6=(Q>=0);

 V=[V1,V2,V3,V4,V5,V6];
 
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

 if (D21+D22*Dck1*D)~=0
     disp('Error : matrix D2 is not 0');
     return
 end

 % mappatura inversa
 Dk=Dck1;
 Ck=(Cck1- Dk*C*Y1)*inv(MT);
 Bk=inv(N)*(Bck1- X1*B1*Dk);
 Ak=inv(N)*[ Ack1- N*Bk*C*Y1-X1*B1*Ck*(MT)- X1*(A+B1*Dk*C)*Y1 ]*inv(MT);



 end