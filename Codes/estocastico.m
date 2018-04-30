clear
%La funcion menu permite que puedas escoger entre opciones
op=menu('Elección de Parámetros','Prestablecido','Personalizado');
if op==1
%Parametro tecnologico
A=1;					
%Share del capital en el producto
ALFA=0.35;	        		
%Factor subjetivo de descuento
BETA=0.95;                    
%Aversión al riesgo (curvatura de la función de utilidad)
SIGMA=2;
%Tasa de depreciacion   	
DELTA=0.1;   				
elseif op==2
%El input permite la introduccion de datos por pantalla
disp(['A continuación introducirá el valor de los parámetros--presione enter para continuar--']);pause;
BETA=input('Introduzca el valor del factor subjetivo de descuento   ');
SIGMA=input('Introduzca el valor de la aversión al riesgo   ');
ALFA=input('Introduzca el valor del share del capital   ');		
DELTA=input('Introduzca el valor de la tasa de depreciacion   ');
A=input('Introduzca el valor del parametro tecnologico   ');
end
disp(['=================================================================='])
disp(['	ITERACION DE LA FUNCION VALOR PARA EL MODELO DE CRECIMIENTO'])
disp([' 		          NEOCLASICO ESTOCASTICO		   '])
disp(['=================================================================='])
disp(['Elaborado por: Kamal Romero'])
disp(['Comentarios por favor enviar a karomero@ucm.es'])
disp(['--presione enter para continuar--']);pause;
%****************************************************************************
%PROCESO ESTOCASTICO;
%Distribución de Probabilidad del parámetro tecnológico
THETA=[1.2 0.8];							
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
%ESTADO ESTACIONARIO;
%Capital estado estacionario determinista
Kss=((1-(1-DELTA)*BETA)/(BETA*ALFA*A))^(1/(ALFA-1));	
%Capital estado estacionario estado bueno
Kssg=((1-(1-DELTA)*BETA)/(BETA*ALFA*THETA(1)))^(1/(ALFA-1));	
%Capital estado estacionario estado malo
Kssb=((1-(1-DELTA)*BETA)/(BETA*ALFA*THETA(2)))^(1/(ALFA-1));	
%Producto estado estacionario estocástico
Yss=A*Kss^ALFA;			
%Consumo estado estacionario estocástico
Css=Yss-DELTA*Kss;
disp(['=================================================================='])
disp(['		Valores de Estado Estacionario		'])
disp(['Capital estado determinista ',num2str(Kss)])
disp(['Capital estado bueno ',num2str(Kssg)])
disp(['Capital estado malo ',num2str(Kssb)])
disp(['=================================================================='])
%****************************************************************************
%GENERACION DEL GRID;
%Nivel minimo de capital en el grid
MINK=0.75*Kssb;                    
%Nivel maximo de capital en el grid
MAXK=1.2*Kssg;                     
%Número de celdas en el grid
NK=input('Introduzca el tamaño del grid (número impar)  ');
%Tamaño del incremento en el grid
INK=((MAXK-MINK)/(NK-1));              
GRID=[MINK:INK:MAXK];
%****************************************************************************
%CALCULO DEL CONSUMO Y LA INVERSION;
I=ones(1,NK);
inversion=GRID'*I-((1-DELTA)*(I'*GRID));
consumo_alto=I'*((THETA(1)*A*GRID).^ALFA)-inversion;
consumo_bajo=I'*((THETA(2)*A*GRID).^ALFA)-inversion;
%****************************************************************************
%CALCULO DE LA UTILIDAD;
U=zeros(NK,NK);
if SIGMA==1
   U_alto=log(consumo_alto);
   U_bajo=log(consumo_bajo);
else
   U_alto=(1/(1-SIGMA))*(consumo_alto.^(1-SIGMA));
   U_bajo=(1/(1-SIGMA))*(consumo_bajo.^(1-SIGMA));
end
%Se le asigna un valor extremadamente bajo a los consumos no factibles
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
tic
figure(1)
n=1;
%Numero maximo de iteraciones
MAXIT=2500;                  
%Nivel de tolerancia
TOL=1*10^(-8);                  
distance1=1;
distance2=1;
while distance1>TOL & distance2>TOL & n<MAXIT;
   %Valor presente funcion valor período siguiente en el estado bueno hoy
   B_alto=BETA*(o*(Q(1,1)*v1+Q(1,2)*v2));			
   %Valor presente funcion valor período siguiente en el estado malo hoy
   B_bajo=BETA*(o*(Q(2,1)*v1+Q(2,2)*v2));			
   Tv1=max(U_alto+(B_alto)');  
   Tv2=max(U_bajo+(B_bajo)');
   V=[Tv1;Tv2];
   distance1=norm(abs(Tv1-v1));
   distance2=norm(abs(Tv2-v2));
   plot(GRID,V)
   xlabel('Capital en t+1')
   ylabel('Función Valor')
   title(['Función Valor Óptima alcanzada luego de ',num2str(n),' iteraciones'])
   pause(0.01)
   v1=Tv1;
   v2=Tv2;
   n=n+1;
      end;
toc
if n>=MAXIT
    disp(['Advertencia: Numero maximo de iteraciones alcanzado ',num2str(n)])
else
    disp(['El algoritmo ha convergido en '...
            num2str(n),' iteraciones'])
    disp(['Operacion realizada en '...
         num2str(toc),' segundos'])
end
%****************************************************************************
%CALCULO FUNCIONES DE POLITICA   
   [value_good,j]=max(U_alto+(B_alto)');
   policy_good_Kap=GRID(j);
   policy_good_Cons=(THETA(1)*A*GRID.^ALFA)+((1-DELTA)*GRID)-(policy_good_Kap);
   [value_bad,jj]=max(U_bajo+(B_bajo)');
   policy_bad_Kap=GRID(jj);
   policy_bad_Cons=(THETA(2)*A*GRID.^ALFA)+((1-DELTA)*GRID)-(policy_bad_Kap);
   PK=[policy_good_Kap;policy_bad_Kap];
%GRAFICOS DE RESULTADOS;
disp('- - Pausa - - Pulse una tecla para ver gráfico 2'); pause;
figure(2)
   plot(GRID, PK, GRID, GRID,':')
   axis tight
   xlabel('Capital en t')
   ylabel('Capital en t+1')
   title(['Regla de Desicion del Capital (Ambos estados)'])
   text(Kss,Kss,'**  Estado Estacionario Determinista')
   text(Kssb,Kssb,'**  Estado Estacionario Estado Malo')
   text(Kssg,Kssg,'**  Estado Estacionario Estado Bueno')
   disp('- - Pausa - - Pulse una tecla para ver gráfico 3'); pause;
figure(3)
   plot(GRID, policy_good_Cons,'-', GRID, policy_bad_Cons,':')
   axis tight
   xlabel('Capital en t')
   ylabel('Consumo en t')
   title(['Regla de Desicion del Consumo'])
   legend('Estado Bueno','Estado Malo',0)