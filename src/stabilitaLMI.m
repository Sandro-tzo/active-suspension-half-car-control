function stabilitaLMI(A)

yalmip('clear');

n = size(A)*[1;0];
X = sdpvar(n, n);
X = (X + X')/2; %  simmetria


v1= ([-X zeros(n,n);
    zeros(n,n)  X*A+A'*X]<=0);

options = sdpsettings('solver', 'sedumi','verbose',1); 
diagnostics = optimize(v1, [], options);


if diagnostics.problem == 0
    disp('Il sistema è as stabile');
    disp('Matrice X trovata:');
    disp(value(X));
elseif diagnostics.problem == 1
    disp('il sistema non è as stabile');
else
    disp('Errore durante la risoluzione del problema.');
end

end
