function [a_optimo]=golden(i,l,ax,bx,cx)
%*---------------------------------------------------------------------*
%golden(i,l,ax,bx,cx,w,r)
%Función que calcula el nivel de activos óptimo en cada intervalo de
%la función valor interpolada
%Inputs
%i: posición en el grid original del activo en el período actual
%l: estado de la naturaleza actual
%ax<bx<cx: entorno del grid original dentro del cual se evalúa bellman
%vn= función valor en el estado n
%w y r: salrio y tipo de interés respectivamente
%Output
%a_optimo: nivel de activo el período siguiente fuera del grid original
%*---------------------------------------------------------------------*

global value_good value_bad w r
%*********************************************************************
%PARAMETROS DE LA ITERACION
%Tolerancia
tol=1*10^(-8);
%Golden Ratios
r1=0.61803399;
r2=1-r1;
%Limites
%En todo momento trabajaremos sobre 4 puntos x0,x1,x3,x4
x0=ax;
x3=cx;
%*********************************************************************
%BRACKET INICIAL
%Aca definimos los otros dos puntos x1 y x2
if abs(cx-bx)<=abs(bx-ax)
x1=bx;
x2=x1+r2*(cx-bx);
else
x2=bx;
x1=x2-r2*(bx-ax);
end
%*********************************************************************
%EVALUACION INICIAL
ii=i;
ll=l;
f1=bellman(ll,ii,x1);
f2=bellman(ll,ii,x2);
%*********************************************************************
%GOLDEN SECTION SEARCH
while abs(x3-x0)>=tol*(abs(x1)+abs(x2))
if f2>f1
x0=x1;
x1=x2;
x2=r1*x1+r2*x2;
f1=f2;
f2=bellman(ll,ii,x2);
else
x3=x2;
x2=x1;
x1=r1*x2+r2*x0;
f2=f1;
f1=bellman(ll,ii,x1);
end
end
if f1>=f2
a_optimo=x1;
else
a_optimo=x2;
end
