t1=clock;
global TAU ALFA BETA SIGMA DELTA GRID NK THETA S Q value_good value_bad P Pk joint
%global indica que se comparten dichas variables para varias funciones
disp(['****************************************************************************'])
disp(['************OBTENCION DE LA DISTRIBUCION ESTACIONARIA DE *******************'])
disp(['*************ACTIVOS EN UN MODELO DE AGENTES HETEROGENEOS*******************'])
disp(['************************CON PRECIOS EXOGENOS********************************'])
disp(['----------------------------------------------------------------------------'])
disp(['*******************Codigo elaborado por: KAMAL ROMERO***********************'])
disp(['***************cualquier comentario por favor escribir a********************'])
disp(['**************************karomero@ucm.es***********************************'])
disp(['****************************************************************************'])
%PARAMETROS DEL MODELO;
op=menu('Elección de Parámetros','Predeterminado','Personalizado');
if op==1
%Parametro tecnologico
TAU=1;					
%Share del capital en el producto
ALFA=0.36;	        		
%Factor subjetivo de descuento
BETA=0.96;
%Aversión al riesgo (curvatura de la función de utilidad)
SIGMA=2;
%Depreciacion
DELTA=0.08;
%Tasa de interes
r=0.03;
%****************************************************************************
%PROCESO ESTOCASTICO;
%Distribución de Probabilidad de las unididades de eficiencia del trabajo
THETA=[1 0.5];							
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
elseif op==2
%Factor subjetivo de descuento
BETA=input('Introduzca el valor del factor subjetivo de descuento   ');
%Aversión al riesgo (curvatura de la función de utilidad)
SIGMA=input('Introduzca el valor del parametro de aversion al riesgo   ');
%****************************************************************************
%PROCESO ESTOCASTICO;
%Distribución de Probabilidad de las unididades de eficiencia del trabajo
THETA=[1 0.5];							
%Número de estados
S=length(THETA);
%Probabilidad de quedarse en el estado b (estando en el estado b)
Pbb=input('Introduzca la probabilidad de quedarse empleado (estando empleado)   ');
%Probabilidad de transitar al estado m (estando en el estado b)
Pbm=1-Pbb;				
%Probabilidad de quedarse en el estado m (estando en el estado m)
Pmm=input('Introduzca la probabilidad de quedarse desempleado (estando desempleado)   ');
%Probabilidad de transitar al estado b (estando en el estado m)
Pmb=1-Pmm;				
%Matriz de transición
Q=[Pbb Pbm; Pmb Pmm];
end
%****************************************************************************
%GENERACION DEL GRID
%Nivel minimo de capital en el grid
MINK=0.1; 
%Nivel maximo de capital en el grid                      
MAXK=20;
%Tamaño del incremento en el grid
INK=0.1;                      
GRID=[MINK:INK:MAXK];
%Número de celdas en el grid
NK=length(GRID);
disp(['=================================================================='])
disp(['		Valores de los Parametros	'])
disp(['Tasa de Interes                             ',num2str(r)])
disp(['Factor subjetivo de descuento               ',num2str(BETA)])
disp(['Grado de aversión al riesgo                 ',num2str(SIGMA)])
disp(['Share del capital en el PIB                 ',num2str(ALFA)])
disp(['Tasa de depreciacion                        ',num2str(DELTA)])
disp(['Parametro tecnologico                       ',num2str(TAU)])
disp(['=================================================================='])
disp(['		Proceso Estocastico	  '])
disp(['Numero de estados ',num2str(S)])
disp(['Valores de las unidades de eficiencia del trabajo '])
disp([''])
disp([      THETA ])
disp(['Matriz de Trancision'])
disp([''])
disp([      Q ])
disp(['=================================================================='])
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
toler=0.0001;
metric=100;
liter=1;
maxiter=20;
ln=1;
distance1=1;
distance2=1;
%Nivel de activos inicial
Astart=5;
disp(['=================================================================='])
disp(['		         Iterando en el Nivel de Activos	'])
disp(['=================================================================='])
disp(['Iteracion   Distancia  Capital nuevo  Capital viejo   '])
while (metric>toler) & (liter <= maxiter);
%****************************************************************************
%CALCULO DEL CONSUMO Y LA UTILIDAD
I=ones(1,NK);
% A=(((r+DELTA)/(ALFA*TAU))*(Mean_L^(ALFA-1)))^(1/(ALFA-1));
% w=(1-ALFA)*TAU*(A^ALFA)*(Mean_L^(-ALFA));
consumo_alto=THETA(1)+((1+r)*(I'*GRID))-GRID'*I;
consumo_bajo=THETA(2)+((1+r)*(I'*GRID))-GRID'*I;
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
%GOLDEN RULE SEARCH BRACKETS
for l=1:2
if l==1
for i=1:NK
   if i==1
      ax=PGK(i);
      bx=PGK(i);
      cx=PGK(i+1);
   elseif i==NK
      ax=PGK(i-1);
      bx=PGK(i);
      cx=PGK(i);
   else
      ax=PGK(i-1);
      bx=PGK(i);
      cx=PGK(i+1);
   end
   policy_good_Kap(i)=golden(i,l,ax,bx,cx);
end
else
   for i=1:NK
   if i==1
      ax=PBK(i);
      bx=PBK(i);
      cx=PBK(i+1);
   elseif i==NK
      ax=PBK(i-1);
      bx=PBK(i);
      cx=PBK(i);
   else
      ax=PBK(i-1);
      bx=PBK(i);
      cx=PBK(i+1);
   end
   policy_bad_Kap(i)=golden(i,l,ax,bx,cx);
end
end
end
%****************************************************************************
%FUNCION DE POLITICA PARA EL CONSUMO
POL_CON_GOOD=THETA(1)+((1+r)*GRID)-policy_good_Kap;
POL_CON_BAD=THETA(2)+((1+r)*GRID)-policy_bad_Kap;
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
%Calculamos el nivel de activos medio como SUM F(a,s)*a
PFKK=[GRID' GRID'];
kkk=PFKK(:);
Mean_A=P'*kkk;
%Calculamos la oferta de activos medio como SUM F(a,s)*g(a,s)
PFK=[policy_good_Kap' policy_bad_Kap'];
kk=PFK(:);
Mean_K=P'*kk;
%Calculamos el consumo medio como SUM F(a,s)*c(a,s)
PFC=[POL_CON_GOOD' POL_CON_BAD'];
cc=PFC(:);
Mean_C=P'*cc;
%Se establecen los limites del grafico de la distribucion de activos
for i=1:length(Pk)
if Pk(i)==0
    break
    ii=i;
end
end
Aold=Astart;
g=0.75;
Anew=g*Mean_K+(1-g)*Aold;
metric=abs((Aold-Mean_K)/Aold);
Astart=Anew;
liter=liter+1;
disp([   liter          metric          Anew                  Aold   ])
end;
timer=etime(clock,t1);
disp(['Encontrado un punto fijo de la distribucion de activos a '])
disp(['un nivel medio de ',num2str(Mean_K),' en un tiempo de ',num2str(timer),' segundos'])
disp(['=================================================================='])
disp(['		Valores de los Agregados	'])
disp(['Capital                         ',num2str(Mean_K)])
disp(['Consumo                         ',num2str(Mean_C)])
%disp(['Demanda Agregada                ',num2str(DA)])
%disp(['PIB neto de depreciacion        ',num2str(OA)])
disp(['=================================================================='])
figure
plot(GRID,Pk)
title(['Distribucion Estacionaria de Activos r= ',num2str(r),' B= ',num2str(BETA),' SIGMA= ',num2str(SIGMA),''])
xlabel('Nivel de Activos');
ylabel('% de Agentes');
axis([GRID(1) GRID(i) min(Pk) max(Pk)])
figure
plot(Pp)
title(['Distribucion Estacionaria Condicionada de Activos r= ',num2str(r),' B= ',num2str(BETA),' SIGMA= ',num2str(SIGMA),''])
xlabel('Nivel de Activos');
ylabel('% de Agentes');
%DISTRIBUCION ESTACIONARIA DEL INGRESO
%Definimos un vector de ingreso (1+r)a + w
Income=[((1+r)*GRID)+THETA(1) ((1+r)*GRID)+THETA(2)];
%Ordenamos el nivel de ingreso en orden ascendente y guardamos su posición en 
%el vector original
[IN index]=sort(Income(:));
%Ordenamos las probabilidades
pjoint = joint(:);
fraction=pjoint(index);
figure
plot(IN,fraction);
title(['Distribucion Estacionaria del Ingreso r= ',num2str(r),' B= ',num2str(BETA),' SIGMA= ',num2str(SIGMA),''])
xlabel('Nivel de Ingreso');
ylabel('% de Agentes');
