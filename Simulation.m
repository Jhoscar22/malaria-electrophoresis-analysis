%Meshgrid
m=3;
dm=0.5;
[xx, yy, zz]=meshgrid(-m:dm:m);

%Constantes
k=8.998e9;
q=1e-6*80;  %microCoulomb
%q = q*1000
%k=1;
%q=1;

%Celula
r = 0.2; %Radio de la esfera
%dc = 2*r; %distancia entre las cargas
dc = 2*r;
qC = 1e-6*1; %Carga de la celula de un microCoulonb por polo

%/////////////////////////////////////////////////////////////////
% Seccion de creacion automatizada de puntos a base de parametros de
%  Distribucion

%Parametros ajustables
Num_Cargas = 5; %Entero impar (funciona con par, pero las cargas negativas toman valores de .25 y -75)
Dist_placas = 3; %Entero par o impar
Dist_pos = 1; %Distancia entre cargas positivas
Dist_neg = 0.5; %Dist_pos/2; %Ditancia entre cargas negativas, por default es la mitad entre positivas
x_ini = 0; %Cualquier valor

%Cargas positivas
A = x_ini*ones(1, Num_Cargas);
B = Dist_placas/2*ones(1, Num_Cargas);
B= -B; %adecuacion de Kevin
C = zeros(1, Num_Cargas);
C_ini = -(Num_Cargas-1)/2*Dist_pos;
cont=1;
for i=C_ini:Dist_pos:-C_ini
    C(1,cont)=i;
    cont = cont + 1;
end

%Cargas negativas
D = x_ini*ones(1, Num_Cargas);
E = Dist_placas/2*ones(1, Num_Cargas);
F = zeros(1, Num_Cargas);
F_ini = -(Num_Cargas-1)/2*Dist_neg;
cont=1;
for i=F_ini:Dist_neg:-F_ini
    F(1,cont)=i;
    cont = cont + 1;
end
%/////////////////////////////////////////////////////////////////

%Cargas positivas
%A = [0 0 0 0];
%B = [-1.5 -1.5 -1.5 -1.5];
%C = [-1 0 1 2];

%Cargas negativas
%D = [0 0 0 0 0];
%E = [1.5 1.5 1.5 1.5 1.5];
%F = [-1 -0.5 0 0.5 1];
%D = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%E = [1.5 1.5 1.5 1.5 1.5 2 2 2 2 2 2.5 2.5 2.5 2.5 2.5];
%F = [-1 -0.5 0 0.5 1 -1 -0.5 0 0.5 1 -1 -0.5 0 0.5 1];
%D = [0 0 0 0 0 0 0 0 0 ]
%E = [1.5 1.75 1.75 1.75 2 2 2 2 2]
%F = [0 -0.25 0 0.25 -0.5 -0.25 0 0.25 0.5]
%D = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%E = [1.5 1.75 1.75 1.75 2 2 2 2 2 2.25 2.25 2.25 2.25 2.25 2.25 2.25 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5];
%F = [0 -0.5 0 0.5 -1 -0.5 0 0.5 1 -1.5 -1 -0.5 0 0.5 1 1.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2];

%Acumuladores de cargas positivas
ExPosi = zeros(size(xx));
EyPosi = zeros(size(yy));
EzPosi = zeros(size(zz));

%Acumuladores de cargas negativas
ExNega = zeros(size(xx));
EyNega = zeros(size(yy));
EzNega = zeros(size(zz));

%Creacion de Campo Electrico de cargas positivas
i=1;
while i<=size(A,2)
    
    x1=xx-A(i);
    y1=yy-B(i);
    z1=zz-C(i);
    
    
    %Asignacion del campo resultante
    ExPosi = ExPosi + k*q.*x1./(x1.^2+y1.^2+z1.^2).^1.5;  % acomodada la función de otra manera                       
    EyPosi = EyPosi + k*q.*y1./(x1.^2+y1.^2+z1.^2).^1.5;
    EzPosi = EzPosi + k*q.*z1./(x1.^2+y1.^2+z1.^2).^1.5;
    
    %dibujando las cargas
    plot3(A(i),B(i),C(i),'r +') %dibujar la carga
    hold on
    plot3(A(i),B(i),C(i),'r o','MarkerSize',15) %dibujar el signo de la carga
    hold on
    
    i=i+1;
end

%Creacion de Campo Electrico de cargas positivas
i=1;
while i<=size(D,2)
    
    x1=xx-D(i);
    y1=yy-E(i);
    z1=zz-F(i);
    
    
    %Asignacion del campo resultante
    ExNega = ExNega + k*q.*x1./(x1.^2+y1.^2+z1.^2).^1.5;  % acomodada la función de otra manera                       
    EyNega = EyNega + k*q.*y1./(x1.^2+y1.^2+z1.^2).^1.5;
    EzNega = EzNega + k*q.*z1./(x1.^2+y1.^2+z1.^2).^1.5;
    
    %dibujando las cargas
    plot3(D(i),E(i),F(i),'b _') %dibujar la carga
    hold on
    plot3(D(i),E(i),F(i),'b o','MarkerSize',15) %dibujar el signo de la carga
    hold on
    
    i=i+1;
end

%Calculo del campo electrico resultante
Fx =ExPosi-ExNega;
Fy =EyPosi-EyNega;
Fz =EzPosi-EzNega;

%Modificacion a vectores unitarios
Exx=Fx./(sqrt(Fx.^2+Fy.^2+Fz.^2));
Eyy=Fy./(sqrt(Fx.^2+Fy.^2+Fz.^2));
Ezz=Fz./(sqrt(Fx.^2+Fy.^2+Fz.^2));

%Graficacion y recorrido de la celula
quiver3(xx,yy,zz,Fx,Fy,Fz,0.5,'b')
xlabel('x');
ylabel('y');
zlabel('z');
    
%vista del grafico
view([120 15])

%Tiempo para poder poner la pantalla completa
pause(.1);

%distancia de recorrido
dr=3.5; %comenzara en (0, 0, dr)

%Declaracion de la esfera
[X, Y, Z] = sphere(10); %Esta comienza en (0i, 0j, 0k)

axis equal
hold on %No recuerdo que hace este aqui

X2 = X * r; %\
Y2 = Y * r; % |>Se adecua la esfera a su radio
Z2 = Z * r; %/

%Posicion de la celula
xc=0;
yc=0;
zc=dr;

%Posicion de la carga positiva
xcp = xc;
ycp = yc+(dc/2);
zcp = zc;

%Posicion de la carga negativa
xcn = xc;
ycn = yc-(dc/2);
zcn = zc;

%esferita de prueba para posicion de cargas
X3 = X * 0.1; %\
Y3 = Y * 0.1; % |>Se adecua la esfera a su radio
Z3 = Z * 0.1; %/

h=-dr;
while h<=dr
    
    %Punto de origen de las cargas positivas y negativas
    px = [xc xc];
    py = [yc+(dc/2) yc-(dc/2)];
    pz = [zc zc];
    Carga = [qC -qC];

    %Campo resultante iniciado en 0
    Vect_Result_celX = 0;
    Vect_Result_celY = 0;
    Vect_Result_celZ = 0;

    %surf(X, Y, Z)
    axis equal
    
    hold on
    surf(X2+xc, Y2+yc, Z2+zc)
    
    %Cargas de prueba
    hold on
    surf(X3+xcp, Y3+ycp, Z3+zcp) %positiva
    hold on
    surf(X3+xcn, Y3+ycn, Z3+zcn) %negativa
    
    h = h+0.25; %Con cada iteracion avanza una unidad hacia abajo en z
    hold on
    
    %Graficacion y recorrido de la celula
    quiver3(xx,yy,zz,Fx,Fy,Fz,0.5,'b')
    
    %vista del grafico
    view([120 15])
    
    %Intervalos de tiempo para el movimiento de la celula
    xlabel('x');
    ylabel('y');
    zlabel('z');
    pause(0.1);
    
    %Calculo del vector resultante sobre la celula
    
    %como las cargas estan a lo largo del eje y, a lo largo de 
    %este se movera la celula
    
    j=1;
    while j<=2 %size(px,2)

        %variables utilizadas para acumular el campo resultante
        Ex1 = 0;
        Ey1 = 0;
        Ez1 = 0;
        Ex2 = 0;
        Ey2 = 0;
        Ez2 = 0;
        %Creacion de Campo Electrico de cargas positivas
        i=1;
        while i<=size(A,2)

            %Calculo de catetos para distancia
            x=px(j)-A(i);
            y=py(j)-B(i);
            z=pz(j)-C(i);

            %calculo de vector resultante de carga postiva
            r=sqrt(x^2+y^2+z^2);

            %Asignacion del campo resultante
            Ex1 = Ex1 + ((k*q*Carga(j).*x)./(r.^1.5));
            Ey1 = Ey1 + ((k*q*Carga(j).*y)./(r.^1.5));
            Ez1 = Ez1 + ((k*q*Carga(j).*z)./(r.^1.5));

            %dibujando las cargas
            %plot3(A(i),B(i),C(i),'r +') %dibujar la carga
            %hold on
            %plot3(A(i),B(i),C(i),'r o','MarkerSize',15) %dibujar el signo de la carga
            %hold on

            i=i+1;
        end

        %Creacion de Campo Electrico de cargas negativa
        i=1;
        while i<=size(D,2)

            %Calculo de catetos para distancia
            x=px(j)-D(i);
            y=py(j)-E(i);
            z=pz(j)-F(i);

            %calculo de vector resultante de carga postiva
            r=sqrt(x^2+y^2+z^2);

            %Asignacion del campo resultante
            Ex2 = Ex2 + ((k*(-q)*Carga(j).*x)./(r.^1.5));
            Ey2 = Ey2 + ((k*(-q)*Carga(j).*y)./(r.^1.5));
            Ez2 = Ez2 + ((k*(-q)*Carga(j).*z)./(r.^1.5));

            i=i+1;
        end

        Ex=Ex1+Ex2 %Campo resultante de una carga
        Ey=Ey1+Ey2
        Ez=Ez1+Ez2

        Vect_Result_celX = Vect_Result_celX + Ex;
        Vect_Result_celY = Vect_Result_celY + Ey;
        Vect_Result_celZ = Vect_Result_celZ + Ez;

        j=j+1;
    end
    
    %if Vect_Result_celX>0 || Vect_Result_celY>0 || Vect_Result_celZ>0
    %    %Modificacion a vectores unitarios
    %    Vect_Result_celX=Vect_Result_celX./(sqrt(Vect_Result_celX.^2+Vect_Result_celY.^2+Vect_Result_celZ.^2));
    %    Vect_Result_celY=Vect_Result_celY./(sqrt(Vect_Result_celX.^2+Vect_Result_celY.^2+Vect_Result_celZ.^2));
    %    Vect_Result_celZ=Vect_Result_celZ./(sqrt(Vect_Result_celX.^2+Vect_Result_celY.^2+Vect_Result_celZ.^2));
    %end
    
    %Posicion de la celula
    xc = xc + Vect_Result_celX;
    yc = yc + Vect_Result_celY;
    zc = zc + Vect_Result_celZ -0.25;
    %zc = -h;

    %Posicion de la carga positiva
    xcp = xc;
    ycp = yc+(dc/2);
    zcp = zc;

    %Posicion de la carga negativa
    xcn = xc;
    ycn = yc-(dc/2);
    zcn = zc;

    hold off
end

%quiver3(xx,yy,zz,Ex,Ey,Ez,0.95)

%dibujando la carga
%hold on
%plot3(0,0,0,'r +') %dibujar la carga
%plot3(0,0,0,'r o','MarkerSize',15) %dibujar el signo de la carga
%hold off