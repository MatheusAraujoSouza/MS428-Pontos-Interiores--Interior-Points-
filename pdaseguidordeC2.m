%Bruno Trevisan Ra:168170
%Matheus Araujo Souza Ra:184145  MS428


function m2=pdaseguidordeC2(A,b,c)
tau=0.99995;
epsi=1*10^-8;


%*******************pontos iniciais***************************************

%calculando o ponto inicial x=A'(AA)^-1*b que satisfaz a relação primal
%Ax=b
[m,n] = size(A);
[R,p]=chol(A*A');
y=R'\(R\b);
x = A' * y;
E2 = 100;
E1 = min([max(x), E2, norm(b,1)/(E2*norm(A,1))]);
E3 = 1 + norm(c,1);

%determinando o valor inicial de x
x = max(x, E1); % valor inicial de x
%determinando o valor inicial de y
y = zeros(m, 1); %valor inicial de y
%vamos colocar um valor provisorio z
%vamos satisfazer as seguintes condições estabelecidas em aula

%********************calculo do valor de z*********************************
z = c; % valor provisorio de z

%fazendo essa substituição todos os calculos efetuados com z irão ser
%realizados em c
% para z > 0
sub = find(z >0);
z(sub,1) = z(sub, 1) + E3;
%para z<=0
%aqui vou criar o vetor com elementos menores ou iguais a zero
%com esse comando eu atribuo zero para as posições que são maiores que ele
% e também guardo a exata posição 
aux = min(z, 0);
%-z se z<=-E3
sub = find(aux <= -E3); % z <= -E3
aux(sub,1) = -z(sub,1);
%E3 se -E3< z <=0
sub = find(aux > -E3); % z >= -E3
aux(sub ,1) = E3;
z=max(aux, z); % valor inicial de z

%************************começando as interações*******************8*******
i=1;
k=0; %contador de interações 
e=ones(n,1);
while i>epsi
if k==1000 % Evita loop Eterno
    fprintf(1,' \n');
    fprintf(1,'Número máximo de iterações: %d \n',k);
    break
end    
    
    
rp = b - A*x;
rd = c - A'*y - z;
gama = trace(diag(x)*diag(z));
sigma=1/sqrt(n);
if gama < 1
   sigma=gama/n;
else
    sigma=1/sqrt(n);
end
mi=(sigma*gama)/n;
rc = mi*e - diag(x)*diag(z)*e;
%*****vamos começar o calculo da direção afim escala******************* 
D = diag(x.^-1)*diag(z);
[R,p] = chol(A*diag(diag(D).^-1)*A');
dy = R\(R'\(rp + A*diag(diag(D).^-1)*rd - A*diag(z.^-1)*rc));
dx = diag(diag(D).^-1)*(A'*dy - rd + diag(x.^-1)*rc);
dz = diag(x.^-1)*(rc - diag(z)*dx);
%***********tamanho do passo calculo dos alpha*************************** 
xdx = -x./dx;
zdz = -z./dz;
xdx(xdx<0) = [];
zdz(zdz<0) = [];
rhop = min(xdx);
rhod = min(zdz);
alphap = min([tau*rhop, 1]);
alphad = min([tau*rhod, 1]);
%**************atualizando as posições********************************** 
x = x + alphap*dx;
y = y + alphad*dy;
z = z + alphad*dz;
%****************8***condições de otimalidade****************************
Fp = norm(b-A*x)/(norm(b)+1) ;
Fd = norm(c-A'*y - z)/(norm(c)+1);
Fa = abs(c'*x - b'*y)/(abs(c'*x)+abs(b'*y)+1);
F = [Fp; Fd; Fa];
i = norm(F, 1);
k = k+1;
end
x
y
z


end



