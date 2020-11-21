%% RAZON DE REDFIELD
clear all
close all
clc

%cargamos las variables 

load 'SALT.mat'
load 'DIC.mat'
load 'NH4.mat'
load 'NO3.mat'

time1=datevec(time);
D=size(DIC);

%% 
%[DIC(1,1,1) NH4(1,1,1) NO3(1,1,1)]

%Razon de Redfield
%los valores estan en mmol/m^3 que equivale a milimol por metro cubico
%pasemos los milimol/m^3 a milimol/litro dividiendo por mil que es la
%molaridad https://sciencing.com/how-to-convert-mg-to-mmoll-13666361.html

DIC_f=DIC./1000; % mmol/l
NH4_f=NH4./1000; % mmol/l
NO3_f=NO3./1000; % mmol/l


%pasamos a mg  a partir de sus masa molares
%notar que el DIC corresponde a la suma de CT = [CO2*] + [HCO3] + [CO32]
% y  ( [CO2*] = [CO2] + [H2CO3])

c=(((12+16*2)) + ((12+1*2+3*16))) + ((12+1+3*16)) + ((12+16*3));  % [CO2] + [H2CO3] + [HCO3] + [CO3]
n1= (14+4*1);  
n2= (14+3*16);

DIC_ff=DIC_f.*c; %DIC en mg/L
NH4_ff=NH4_f.*n1; %NH4 en mg/L
NO3_ff=NO3_f.*n2; %NO3 en mg/L

%extraemos la cantidad(o porcentaje) de carbono y nitrogeno que hay en cada
%compuesto

%porcentajes de N y C
C=((12/(12+16*2)))  +(12/(12+1*2+3*16)) + (12/(12+1+3*16)) + (12/(12+16*3));  % [C]/[CO2] + [C]/[H2CO3] + [C]/[HCO3] + [C]/[CO3]
N1= 14/(14+4*1);  %que corresponde al porcentaje de nitrogeno en el amonio
N2= 14/(14+3*16);  %que corresponde al porcentaje de nitrogeno en el nitrato (wikipedia jeje)

DIC_fff=DIC_ff.*(C); % mg/L de carbono, no estoy seguro si esta bien hacer el 1-C
NH4_fff=NH4_ff.*N1; %mg/L de nitrogeno
NO3_fff=NO3_ff.*N2; %mg/L de nitrogeno

%volvemos a mmol

DIC_ffff=DIC_fff./12;   %nmoles de carbono en DIC
NH4_ffff=NH4_fff./14;   %mmoles de nitrogeno en NH4
NO3_ffff=NO3_fff./14;   %mmoles de nitrogeno del NO3


R1=DIC_ffff./(NH4_ffff+NO3_ffff);



%tomamos una lat
for i=1:D(3)
    R1_5=squeeze(R1(:,115,:)); %5
    R1_10=squeeze(R1(:,96,:)); %10
    R1_15=squeeze(R1(:,77,:)); %15
    R1_25=squeeze(R1(:,49,:)); %25
    R1_35=squeeze(R1(:,30,:)); %35
    R1_45=squeeze(R1(:,11,:)); %45
end


figure()
pcolor(lon(:,1)-360,time,R1_5')
datetick('y')
axis([lon(1,1)-360 lon(61,1)-360 715885 731945])
% caxis([-3 3])
 shading interp
% colormap(cmap)
colorbar;
% hold on
% contour(lon(:,1)-360,time,NO3_a_5',[0 0],'k','Linewidth',1)
 xlabel('Longitud °W')
 ylabel('Tiempo')
title('Razón de Redfield en 5°S','FontSize', 22)


figure()
pcolor(lon(:,1)-360,time,R1_10')
datetick('y')
axis([lon(1,1)-360 lon(64,1)-360 715885 731945])
% caxis([-3 3])
 shading interp
% colormap(cmap)
colorbar;
% hold on
% contour(lon(:,1)-360,time,NO3_a_5',[0 0],'k','Linewidth',1)
 xlabel('Longitud °W')
 ylabel('Tiempo')
title('Razón de Redfield en 10°S','FontSize', 22)


figure()
pcolor(lon(:,1)-360,time,R1_15')
datetick('y')
axis([lon(1,1)-360 lon(66,1)-360 715885 731945])
% caxis([-3 3])
 shading interp
% colormap(cmap)
colorbar;
% hold on
% contour(lon(:,1)-360,time,NO3_a_5',[0 0],'k','Linewidth',1)
 xlabel('Longitud °W')
 ylabel('Tiempo')
title('Razón de Redfield en 15°S','FontSize', 22)


figure()
pcolor(lon(:,1)-360,time,R1_25')
datetick('y')
axis([lon(1,1)-360 lon(71,1)-360 715885 731945])
% caxis([-3 3])
 shading interp
% colormap(cmap)
colorbar;
% hold on
% contour(lon(:,1)-360,time,NO3_a_5',[0 0],'k','Linewidth',1)
 xlabel('Longitud °W')
 ylabel('Tiempo')
title('Razón de Redfield en 25°S','FontSize', 22)


figure()
pcolor(lon(:,1)-360,time,R1_35')
datetick('y')
axis([lon(1,1)-360 lon(69,1)-360 715885 731945])
% caxis([-3 3])
 shading interp
% colormap(cmap)
colorbar;
% hold on
% contour(lon(:,1)-360,time,NO3_a_5',[0 0],'k','Linewidth',1)
 xlabel('Longitud °W')
 ylabel('Tiempo')
title('Razón de Redfield en 35°S','FontSize', 22)


figure()
pcolor(lon(:,1)-360,time,R1_45')
datetick('y')
axis([lon(1,1)-360 lon(68,1)-360 715885 731945])
% caxis([-3 3])
 shading interp
% colormap(cmap)
colorbar;
% hold on
% contour(lon(:,1)-360,time,NO3_a_5',[0 0],'k','Linewidth',1)
 xlabel('Longitud °W')
 ylabel('Tiempo')
title('Razón de Redfield en 45°S','FontSize', 22)


%% Anomalias de la razon 
R_s=106/16;

R_a=R1-R_s;



%tomamos una lat 
for i=1:D(3)
    R_a_5=squeeze(R_a(:,115,:)); %5
    R_a_10=squeeze(R_a(:,96,:)); %10
    R_a_15=squeeze(R_a(:,77,:)); %15
    R_a_25=squeeze(R_a(:,49,:)); %25
    R_a_35=squeeze(R_a(:,30,:)); %35
    R_a_45=squeeze(R_a(:,11,:)); %45
end


figure()
pcolor(lon(:,1)-360,time,R_a_5')
datetick('y')
axis([lon(1,1)-360 lon(61,1)-360 715885 731945])
% caxis([-3 3])
 shading interp
% colormap(cmap)
colorbar;
% hold on
% contour(lon(:,1)-360,time,NO3_a_5',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Razón de Redfield en 5°S','FontSize', 22)


figure()
pcolor(lon(:,1)-360,time,R_a_10')
datetick('y')
axis([lon(1,1)-360 lon(64,1)-360 715885 731945])
% caxis([-3 3])
 shading interp
% colormap(cmap)
colorbar;
% hold on
% contour(lon(:,1)-360,time,NO3_a_5',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Razón de Redfield en 10°S','FontSize', 22)


figure()
pcolor(lon(:,1)-360,time,R_a_15')
datetick('y')
axis([lon(1,1)-360 lon(66,1)-360 715885 731945])
% caxis([-3 3])
 shading interp
% colormap(cmap)
colorbar;
% hold on
% contour(lon(:,1)-360,time,NO3_a_5',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Razón de Redfield en 15°S','FontSize', 22)


figure()
pcolor(lon(:,1)-360,time,R_a_25')
datetick('y')
axis([lon(1,1)-360 lon(71,1)-360 715885 731945])
% caxis([-3 3])
 shading interp
% colormap(cmap)
colorbar;
% hold on
% contour(lon(:,1)-360,time,NO3_a_5',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Razón de Redfield en 25°S','FontSize', 22)




figure()
pcolor(lon(:,1)-360,time,R_a_35')
datetick('y')
axis([lon(1,1)-360 lon(69,1)-360 715885 731945])
% caxis([-3 3])
 shading interp
% colormap(cmap)
colorbar;
% hold on
% contour(lon(:,1)-360,time,NO3_a_5',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Razón de Redfield en 5°S','FontSize', 22)


figure()
pcolor(lon(:,1)-360,time,R_a_45')
datetick('y')
axis([lon(1,1)-360 lon(68,1)-360 715885 731945])
% caxis([-3 3])
 shading interp
% colormap(cmap)
colorbar;
% hold on
% contour(lon(:,1)-360,time,NO3_a_5',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Razón de Redfield en 45°S','FontSize', 22)




%% anomalias climatologicas

%le quitamos el ciclo anual y estandarizamos 
R_ac=cicloanual_estan(R1,time);


VAR = R_ac; 
[X,Y,T] = size(VAR);
% Filtraje de la variable
VAR_ft  = zeros(X,Y,T);   % Mejora numerica [Variable filtrada]
P       = 11;               % Periodo
% Ciclo iterativo
for i = 1:X
    for j = 1:Y
        VAR_ft(i,j,:) = lanczosfilter(squeeze(VAR(i,j,:)),1,1/P,[],'low');
    end
end

%segunda pasada
VAR1 = VAR_ft; 
[X,Y,T] = size(VAR);
% Filtraje de la variable
VAR_ft1  = zeros(X,Y,T);   % Mejora numerica [Variable filtrada]
P       = 11;               % Periodo
% Ciclo iterativo
for i = 1:X
    for j = 1:Y
        VAR_ft1(i,j,:) = lanczosfilter(squeeze(VAR1(i,j,:)),1,1/P,[],'low');
    end
end

R_ac=VAR_ft1;

%tomamos una lat
for i=1:D(3)
    R_ac_5=squeeze(R_ac(:,115,:)); %5
    R_ac_10=squeeze(R_ac(:,96,:)); %10
    R_ac_15=squeeze(R_ac(:,77,:)); %15
    R_ac_25=squeeze(R_ac(:,49,:)); %25
    R_ac_35=squeeze(R_ac(:,30,:)); %35
    R_ac_45=squeeze(R_ac(:,11,:)); %45
end
cmap = colormap_cpt('Optimus_Prime.cpt',265);


figure()
pcolor(lon(:,1)-360,time,R_ac_5')
datetick('y')
axis([lon(1,1)-360 lon(61,1)-360 715885 731945])
% caxis([-3 3])
 shading interp
colormap(cmap)
colorbar;
% hold on
% contour(lon(:,1)-360,time,NO3_a_5',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías climatologica de la razón de Redfield en 5°S','FontSize', 15)


figure()
pcolor(lon(:,1)-360,time,R_ac_10')
datetick('y')
axis([lon(1,1)-360 lon(64,1)-360 715885 731945])
% caxis([-3 3])
 shading interp
colormap(cmap)
colorbar;
% hold on
% contour(lon(:,1)-360,time,NO3_a_5',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías climatologica de la razón de Redfield en 10°S','FontSize', 15)

figure()
pcolor(lon(:,1)-360,time,R_ac_15')
datetick('y')
axis([lon(1,1)-360 lon(66,1)-360 715885 731945])
% caxis([-3 3])
 shading interp
colormap(cmap)
colorbar;
% hold on
% contour(lon(:,1)-360,time,NO3_a_5',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías climatologica de la razón de Redfield en 15°S','FontSize', 15)

figure()
pcolor(lon(:,1)-360,time,R_ac_25')
datetick('y')
axis([lon(1,1)-360 lon(71,1)-360 715885 731945])
% caxis([-3 3])
 shading interp
colormap(cmap)
colorbar;
% hold on
% contour(lon(:,1)-360,time,NO3_a_5',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías climatologica de la razón de Redfield en 25°S','FontSize', 15)


figure()
pcolor(lon(:,1)-360,time,R_ac_35')
datetick('y')
axis([lon(1,1)-360 lon(69,1)-360 715885 731945])
% caxis([-3 3])
 shading interp
colormap(cmap)
colorbar;
% hold on
% contour(lon(:,1)-360,time,NO3_a_5',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías climatologica de la razón de Redfield en 35°S','FontSize', 15)

figure()
pcolor(lon(:,1)-360,time,R_ac_45')
datetick('y')
axis([lon(1,1)-360 lon(68,1)-360 715885 731945])
% caxis([-3 3])
 shading interp
colormap(cmap)
colorbar;
% hold on
% contour(lon(:,1)-360,time,NO3_a_5',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías climatologica de la razón de Redfield en 45°S','FontSize', 15)

