function [f]=bellman(l,i,x)
%*---------------------------------------------------------------------*
%bellman(x,i,l)
%Función que calcula el lado derecho de la ecuación de Bellman para cada x
%inputs:
%x: valor del activo el período siguiente
%i: posición en el grid original del activo en el período actual
%l: estado de la naturaleza actual
%Output
%f: Valor del lado derecho de la ecuación de Bellman
%*---------------------------------------------------------------------*

global THETA Q GRID SIGMA ALFA BETA value_good value_bad w r
%**************************************************************
c=THETA(l)*w+((1+r)*GRID(i))-x;
if c<=0
   U=-10000;
else
   U=1/(1-SIGMA)*(c^(1-SIGMA));
end
v1=value_good;
v2=value_bad;
V1=interp1(GRID,v1,x);
V2=interp1(GRID,v2,x);
f=U+(BETA*(Q(1,1)*V1+Q(1,2)*V2));
