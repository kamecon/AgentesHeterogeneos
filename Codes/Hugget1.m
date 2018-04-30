function g = Hugget1(x)
global TAU ALFA BETA SIGMA DELTA GRID NK THETA S Q value_good value_bad P Pk joint
%global indica que se comparten dichas variables para varias funciones
%****************************************************************************
%CONSTRUIMOS MATRICES INICIALES
%Esto se realiza para que MATLAB utilice menos memoria al construir las
%matrices que vienen a continuación
consumo_alto=zeros(NK,NK);
consumo_bajo=zeros(NK,NK);
U_alto=zeros(NK,NK);
U_bajo=zeros(NK,NK);
Tv1=zeros(1,NK);
Tv2=zeros(1,NK);
v1=zeros(1,NK);
v2=zeros(1,NK);
policy_good_Kap=zeros(1,NK);
policy_bad_Kap=zeros(1,NK);
%****************************************************************************
%PARAMETROS DE LA ITERACION
%Numero maximo de iteraciones
MAXIT=2500;
MAXIT2=15;
%Nivel de tolerancia
TOL2=1*10^(-8);
TOL=1*10^(-4);
ln=1;
distance1=1;
distance2=1;
r=x+1;
%****************************************************************************
%CALCULO DEL CONSUMO Y LA UTILIDAD
I=ones(1,NK);
consumo_alto=THETA(1)+(r*(I'*GRID))-GRID'*I;
consumo_bajo=THETA(2)+(r*(I'*GRID))-GRID'*I;
 if SIGMA==1
U_alto=log(consumo_alto);
U_bajo=log(consumo_bajo);
else
   U_alto=1/(1-SIGMA)*(consumo_alto.^(1-SIGMA));
   U_bajo=1/(1-SIGMA)*(consumo_bajo.^(1-SIGMA));
end;
for i=1:NK
   for j=1:NK
      if consumo_alto(i,j)<=0;
         U_alto(i,j)=-realmax;
      end;
      if consumo_bajo(i,j)<=0;
         U_bajo(i,j)=-realmax;
      end;
   end;
end;
%****************************************************************************
%FUNCION VALOR INICIAL;
v1=zeros(1,NK);
v2=zeros(1,NK);
o=ones(NK,1);
%Funcion valor en el estado bueno
Tv1=max(U_alto+(o*BETA*v1)');		
%Funcion valor en el estado malo
Tv2=max(U_bajo+(o*BETA*v2)');				
v1=Tv1;
v2=Tv2;
%****************************************************************************
%ITERACION DE LA FUNCION VALOR;
n=1;
while distance1>TOL2 distance2>TOL2 & n<MAXIT;
   %Valor presente funcion valor período siguiente en el estado bueno hoy
   B_alto=BETA*(o*(Q(1,1)*v1+Q(1,2)*v2));			
   %Valor presente funcion valor período siguiente en el estado malo hoy
   B_bajo=BETA*(o*(Q(2,1)*v1+Q(2,2)*v2));			
  	Tv1=max(U_alto+(B_alto)');
   Tv2=max(U_bajo+(B_bajo)');
   distance1=norm(abs(Tv1-v1));
   distance2=norm(abs(Tv2-v2));
   v1=Tv1;
   v2=Tv2;
   n=n+1;
end;
if n>=MAXIT
     disp(['Advertencia: Numero maximo de iteraciones alcanzado ',num2str(n)])
  end;
[value_good,PGKINDEX]=max(U_alto+(B_alto)');
[value_bad,PBKINDEX]=max(U_bajo+(B_bajo)');
PGK=GRID(PGKINDEX);
PBK=GRID(PBKINDEX);
%****************************************************************************
%FUNCION DE POLITICA PARA EL CONSUMO
%POL_CON_GOOD=THETA(1)+(r*GRID)-policy_good_Kap;
%POL_CON_BAD=THETA(2)+(r*GRID)-policy_bad_Kap;
%****************************************************************************
%CONSTRUCCION DE LA MATRIZ DE TRANCISION CONJUNTA
g1=sparse(NK,NK);
g2=sparse(NK,NK);
for i=1:NK
g1(i,PGKINDEX(i))=1;
g2(i,PBKINDEX(i))=1;
end;
trans_joint=[Q(1,1)*g1 Q(1,2)*g1; Q(2,1)*g2 Q(2,2)*g2]';
%****************************************************************************
%DISTRIBUCION ESTACIONARIA DE ACTIVOS
%Distribución conjunta inicial 
P=(1/(S*NK))*ones(S*NK,1);
distance_distr=1;
while distance_distr>TOL;
   P1=trans_joint*P;
   distance_distr=norm(abs(P1-P));
   P=P1;
end
%Se localizan los puntos para los cuales la masa en puntos de la
%distribucion son menores a una tolerancia
II=find(P<=TOL);
%Se crea un vector con estos puntos para realizar la suma
PP=zeros(1,length(II));
for i=1:length(II)
    PP(i)=P(II(i));
end
%Se reparte la suma de los puntos obtenidos anteriormente entre el resto de
%la distribucion: length(P)-length(II)
Pp=zeros(S*NK,1);
Rep=sum(PP)/(length(P)-length(II));
for i=1:S*NK
    if P(i)<=TOL;
        Pp(i)=0;
    else
        Pp(i)=P(i)+Rep;
    end
end
joint=zeros(NK,S);
joint(:)=Pp;
Pk=sum(joint');
%Calculamos la oferta de activos medio como SUM F(a,s)*g(a,s)
PFK=[PGK' PBK'];
kk=PFK(:);
Mean_K=P'*kk;
g=Mean_K;