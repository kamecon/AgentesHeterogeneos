clear
%La funcion menu permite que puedas escoger entre opciones
op=menu('Elección de Parámetros','Prestablecido','Personalizado');
if op==1
%Parametro tecnologico
A=10;					
%Share del capital en el producto
ALFA=0.35;	        		
%Factor subjetivo de descuento
BETA=0.95;                    
%Aversión al riesgo (curvatura de la función de utilidad)
SIGMA=2;
%Tasa de depreciacion   	
DELTA=0.06;   				
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
disp([' 		          NEOCLASICO CON DEPRECIACION		   '])
disp(['=================================================================='])
disp(['Elaborado por: Kamal Romero'])
disp(['Comentarios por favor enviar a karomero@ucm.es'])
disp(['--presione enter para continuar--']);pause;
%************************************************************************************************
%ESTADO ESTACIONARIO;
%Capital estado estacionario
Kss=((1-(1-DELTA)*BETA)/(BETA*ALFA*A))^(1/(ALFA-1));	
%Producto estado estacionario
Yss=A*Kss^ALFA;			
%Consumo estado estacionario			
Css=Yss-DELTA*Kss;					
disp(['=================================================================='])
disp(['		Valores de Estado Estacionario		'])
disp(['Capital ',num2str(Kss)])
disp(['Consumo ',num2str(Css)])
disp(['Producto ',num2str(Yss)])
disp(['=================================================================='])
%************************************************************************************************
%GENERACION DEL GRID;
%Nivel minimo de capital en el grid
MINK=0.95*Kss;                     
%Nivel maximo de capital en el grid
MAXK=1.01*Kss;                     
%Número de celdas en el grid
NK=input('Introduzca el tamaño del grid (número impar)  ');			    	
%Tamaño del incremento en el grid
INK=((MAXK-MINK)/(NK-1));              
GRID=[MINK:INK:MAXK];
%************************************************************************************************
%CALCULO DEL CONSUMO Y LA INVERSION;
I=ones(1,NK);
inversion=GRID'*I-((1-DELTA)*(I'*GRID));
consumo=I'*(A*GRID.^ALFA)-inversion;
%************************************************************************************************
%CALCULO DE LA UTILIDAD;
U=zeros(NK,NK);
if SIGMA==1
   U=log(consumo);
else
   U=(1/(1-SIGMA))*(consumo.^(1-SIGMA));
end
%Le asignamos un valor extremamente negativo a los consumos no factibles,
%de manera tal que el agente nunca los escoga como optimos
for i=1:NK
    for j=1:NK
       if consumo(i,j)<=0;
          U(i,j)=-realmax;
       end;
    end;
 end;
%************************************************************************************************
%FUNCION VALOR INICIAL;
v1=zeros(1,NK);
o=ones(NK,1);
Tv1=max(U+(o*BETA*v1)');		
%************************************************************************************************
%ITERACION DE LA FUNCION VALOR;
tic
n=1;
%Numero maximo de iteraciones
MAXIT=2500;                  
%Nivel de tolerancia
TOL=1*10^(-3);                  
distance1=1;
while distance1>TOL & n<MAXIT;
   Tv1=max(U+(o*BETA*v1)');  
   distance1=norm(abs(Tv1-v1));
   plot(GRID,Tv1)
   grid
   xlabel('Capital en t+1')
   ylabel('Función Valor')
   title(['Función Valor Óptima alcanzada luego de ',num2str(n),' iteraciones con error '...
         num2str(norm(abs(Tv1-v1)))])
   pause(0.01)
   v1=Tv1;
   n=n+1;
      end;
toc
%************************************************************************************************
%GRAFICOS DE RESULTADOS;
	[value_function,j]=max(U+(o*BETA*v1)');
   policy_function_Kap=GRID(j);
   policy_function_Cons=(A*GRID.^ALFA)+(1-DELTA)*GRID-(policy_function_Kap);
   if n>=MAXIT
     disp(['Advertencia: Numero maximo de iteraciones alcanzado ',num2str(n)])
 else
   disp(['El algoritmo ha convergido en '...
            num2str(n),' iteraciones'])
   disp(['Operacion realizada en '...
         num2str(toc),' segundos'])
   figure(1)
plot(GRID,Tv1)
   grid
   xlabel('Capital en t+1')
   ylabel('Función Valor')
   title(['Función Valor Óptima alcanzada luego de ',num2str(n),' iteraciones'])
text(mean(GRID),mean(Tv1),'     Presione enter para ver grafico 2')
end;
pause;
   figure(2)
%El ':' hace gráficos con líneas discontinuas
   plot(GRID,policy_function_Kap, GRID, GRID,':') 
   xlabel('Capital en t')
   ylabel('Capital en t+1')
axis tight
   title(['Regla de Desicion del Capital'])
   text(Kss,Kss,'**   Capital Estado Estacionario')
   text(Kss-2.5,Kss-2.5,'   Presione enter para ver grafico 3')
pause;
    figure(3)
   plot(GRID,policy_function_Cons)
   xlabel('Capital en t')
   ylabel('Consumo en t')
axis tight
	title(['Regla de Desicion del Consumo'])
