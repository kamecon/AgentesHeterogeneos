clear
t1=clock;
global TAU ALFA BETA SIGMA DELTA GRID NK THETA S Q value_good value_bad P Pk joint Mean_L POL_CON_GOOD POL_CON_BAD policy_good_Kap policy_bad_Kap Pp
disp(['****************************************************************************'])
disp(['************OBTENCION DE LA DISTRIBUCION ESTACIONARIA DE *******************'])
disp(['*************ACTIVOS EN UN MODELO DE AGENTES HETEROGENEOS*******************'])
disp(['****************CON PRECIOS ENDOGENOS CON PRODUCCION************************'])
disp(['----------------------------------------------------------------------------'])
disp(['*******************Codigo elaborado por: KAMAL ROMERO***********************'])
disp(['***************cualquier comentario por favor escribir a********************'])
disp(['**************************karomero@ucm.es***********************************'])
disp(['****************************************************************************'])
%global indica que se comparten dichas variables para varias funciones
%PARAMETROS DEL MODELO;
%Parametro tecnologico
TAU=1;					
%Share del capital en el producto
ALFA=0.36;	        		
%Factor subjetivo de descuento
BETA=0.96;
%Aversión al riesgo (curvatura de la función de utilidad)
SIGMA=2;
%Depreciacion
DELTA=0.03;
%****************************************************************************
%PROCESO ESTOCASTICO;
%Distribución de Probabilidad de las unididades de eficiencia del trabajo
THETA=[1 0.05];							
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
%TRABAJO AGREGADO
%calculamos la matriz de transicion ergodica para poder calcular la 
%distribucion inicial del trabajo. Sargent-Ljungqvist pag.367
%Matriz de Transicion Inicial
p=Q;
%Distribucion Inicial de las unidades de eficiencia
g=[0.5,0.5];
distance=1;
toll=1*10^(-8);
while distance>toll
p1=p*p;
g1=g*p;
distance=norm(abs(p1-p));
p=p1;
g=g1;
end
%g es la distribucion invariante de las unidades de eficiencia entre 
%los trabajadores
%****************************************************************************
%TRABAJO AGREGADO EXOGENO
Mean_L=THETA*g';
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
%****************************************************************************
disp(['=================================================================='])
disp(['		Valores de los Parametros	'])
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
disp(['Distribucion invariante de las unidades de eficiencia'])
disp([''])
disp([      g ])
disp(['=================================================================='])
%****************************************************************************
%CAPITAL DE EQUILIBRIO
disp(['INICIANDO LA LOCALIZACION DE LA TASA DE INTERES DE EQUILIBRIO CON '])
disp(['                      LA FUNCION FZERO '])
%Se busca la tasa de interes de equilibrio como aquella que vacia el 
%mercado de activos, o lo que es lo mismo, activos agregados igual a cero
R=fzero('aiyagari1',0.059,optimset('disp','iter'))
%****************************************************************************
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
%Calculamos el nivel de producto  y demanda agregada implicado por las
%reglas de decision
%Oferta agregada Y=A*[(K^ALFA)*(L^(1-ALFA))] + (1-DELTA)*K
OA=TAU*(Mean_A^ALFA)*(Mean_L^(1-ALFA))+((1-DELTA)*Mean_A);
%Demanda agregada SUM F(a,s)*g(a,s) + SUM F(a,s)*c(a,s)
DA=Mean_C+Mean_K;
K=Mean_L*(ALFA/(R+DELTA)).^(1/(1-ALFA));
w=(1-ALFA)*TAU*(K^ALFA)*(Mean_L^(-ALFA));
timer=etime(clock,t1);
disp(['Encontrado una tasa de interes que reproduce el nivel de capital ',num2str(Mean_K)])
disp(['a un nivel de ',num2str(R),' en un tiempo de ',num2str(timer),' segundos'])
disp(['=================================================================='])
disp(['		Valores de los Agregados	'])
disp(['Tasa de Interes                 ',num2str(R)])
disp(['Salario                         ',num2str(w)])
disp(['Capital                         ',num2str(Mean_K)])
disp(['Consumo                         ',num2str(Mean_C)])
disp(['Demanda Agregada                ',num2str(DA)])
disp(['PIB neto de depreciacion        ',num2str(OA)])
disp(['=================================================================='])
r=R;
figure
plot(GRID,Pk)
title(['Distribucion Estacionaria de Activos r= ',num2str(r),' B= ',num2str(BETA),' SIGMA= ',num2str(SIGMA),''])
xlabel('Nivel de Activos');
ylabel('% de Agentes');
figure
plot(Pp)
title(['Distribucion Estacionaria Condicionada de Activos r= ',num2str(r),' B= ',num2str(BETA),' SIGMA= ',num2str(SIGMA),''])
xlabel('Nivel de Activos');
ylabel('% de Agentes');
%DISTRIBUCION ESTACIONARIA DEL INGRESO
%Definimos un vector de ingreso (1+r)a + w
Income=[((1+r)*GRID)+w*THETA(1) ((1+r)*GRID)+w*THETA(2)];
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