clear
t1=clock;
global TAU ALFA BETA SIGMA DELTA GRID NK THETA S Q value_good value_bad P Pk joint
disp(['****************************************************************************'])
disp(['************OBTENCION DE LA DISTRIBUCION ESTACIONARIA DE *******************'])
disp(['*************ACTIVOS EN UN MODELO DE AGENTES HETEROGENEOS*******************'])
disp(['****************CON PRECIOS ENDOGENOS SIN PRODUCCION************************'])
disp(['----------------------------------------------------------------------------'])
disp(['*******************Codigo elaborado por: KAMAL ROMERO***********************'])
disp(['***************cualquier comentario por favor escribir a********************'])
disp(['**************************karomero@ucm.es***********************************'])
disp(['****************************************************************************'])
%PARAMETROS DEL MODELO;
%Parametro tecnologico
TAU=1;					
%Share del capital en el producto
ALFA=0.36;	        		
%Factor subjetivo de descuento
BETA=0.96;
%Aversión al riesgo (curvatura de la función de utilidad)
SIGMA=1.5;
%Depreciacion
DELTA=0.08;
%****************************************************************************
%PROCESO ESTOCASTICO;
%Distribución de Probabilidad de las unididades de eficiencia del trabajo
THETA=[1 0.25];							
%Número de estados
S=length(THETA);
%Probabilidad de quedarse en el estado b (estando en el estado b)
Pbb=0.8;					
%Probabilidad de transitar al estado m (estando en el estado b)
Pbm=1-Pbb;				
%Probabilidad de quedarse en el estado m (estando en el estado m)
Pmm=0.5;					
%Probabilidad de transitar al estado b (estando en el estado m)
Pmb=1-Pmm;				
%Matriz de transición
Q=[Pbb Pbm; Pmb Pmm];
%****************************************************************************
%GENERACION DEL GRID
%Nivel minimo de capital en el grid
MINK=-2; 
%Nivel maximo de capital en el grid                      
MAXK=10;
%Tamaño del incremento en el grid
INK=0.1;                      
GRID=[MINK:INK:MAXK];
%Número de celdas en el grid
NK=length(GRID);
%****************************************************************************
%TASA DE INTERES DE EQUILIBRIO
disp(['INICIANDO LA LOCALIZACION DE LA TASA DE INTERES DE EQUILIBRIO CON '])
disp(['                      LA FUNCION FZERO '])
%Se busca la tasa de interes de equilibrio como aquella que vacia el 
%mercado de activos, o lo que es lo mismo, activos agregados igual a cero
R=fzero('Hugget1',0.03,optimset('disp','iter'))
%Una vez localizada la tasa de interes mostramos los resultados
timer=etime(clock,t1);
disp(['Encontrado una tasa de interes que vacia el mercado de activos a '])
disp(['un nivel de ',num2str(R),' en un tiempo de ',num2str(timer),' segundos'])
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
%****************************************************************************
%CALCULO DEL CONSUMO Y LA UTILIDAD
I=ones(1,NK);
consumo_alto=THETA(1)+((1+R)*(I'*GRID))-GRID'*I;
consumo_bajo=THETA(2)+((1+R)*(I'*GRID))-GRID'*I;
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
POL_CON_GOOD=THETA(1)+((1+R)*GRID)-PGK;
POL_CON_BAD=THETA(2)+((1+R)*GRID)-PBK;
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
%Calculamos el consumo medio como SUM F(a,s)*c(a,s)
PFC=[POL_CON_GOOD' POL_CON_BAD'];
cc=PFC(:);
Mean_C=P'*cc;
%Calculamos la oferta de activos medio como SUM F(a,s)*g(a,s)
PFK=[PGK' PBK'];
kk=PFK(:);
Mean_K=P'*kk;
%Se establecen los limites del grafico de la distribucion de activos
for i=1:length(Pk)
if Pk(i)==0
    break
    ii=i;
end
end
GG=find(GRID>0);
gg=GG(1)-1;
indebt=sum(Pk(1:gg));
disp(['=================================================================='])
disp(['Medida de Agentes Endeudados '           ,num2str(indebt)])
disp(['Medida de Agentes Borrowed Constrained ' ,num2str(Pk(1))])
disp(['=================================================================='])
figure
plot(GRID,Pk)
title(['Distribucion Estacionaria de Activos r= ',num2str(R),' B= ',num2str(BETA),' SIGMA= ',num2str(SIGMA),''])
xlabel('Nivel de Activos');
ylabel('% de Agentes');
axis([GRID(1) GRID(i) min(Pk) max(Pk)])
figure
plot(Pp)
title(['Distribucion Estacionaria Condicionada de Activos r= ',num2str(R),' B= ',num2str(BETA),' SIGMA= ',num2str(SIGMA),''])
xlabel('Nivel de Activos');
ylabel('% de Agentes');
