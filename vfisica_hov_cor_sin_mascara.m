%% HOVMOLLERS de SST, NIVEL DEL MAR, TEMPERATUR Y OXIGENO Y CORRELACION CON ESTOS Y EL ONI
%% PERO SIN LA MASCARA, ES DECIR VAN A TENER TODAS LAS LONGITUDES EN LATITUDES MAS ALTAS
%% ADEMAS SERAN OTRAS 3 LATITUDES RESPECTO A EL AVANCE ANTERIOR ES DECIR 5, 10, 15, 25, 35 Y 40 °S
%% y de las variables filtradas de variabilidad interanual
clear all
close all
clc

%% oxigeno 

load 'ONI.mat'
load 'o2.mat'
load 'sam.mat'
time1=datevec(time);

D=size(o2);
%le quitamos el ciclo anual y estandarizamos 
o2_a=cicloanual_estan(o2,time);


%filtro de 11 meses usando el metodo coseno-lanczos
%primera pasada
VAR = o2_a; 
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

o2_af=VAR_ft1;



%tomamos una lat 
o2_a_5=squeeze(o2_af(:,115,:)); %5
o2_a_10=squeeze(o2_af(:,96,:)); %10
o2_a_15=squeeze(o2_af(:,77,:)); %15
o2_a_25=squeeze(o2_af(:,49,:)); %25
o2_a_35=squeeze(o2_af(:,30,:)); %35
o2_a_45=squeeze(o2_af(:,11,:)); %45

cmap = colormap_cpt('Optimus_Prime.cpt',265);


o2_a_5=o2_a_5';
o2_a_10=o2_a_10';
o2_a_15=o2_a_15';
o2_a_25=o2_a_25';
o2_a_35=o2_a_35';
o2_a_45=o2_a_45';

figure()
pcolor(lon(:,1)-360,time,o2_a_5)
datetick('y')
axis([lon(1,1)-360 lon(61,1)-360 715885 731945])
caxis([-2.5 2.5])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,o2_a_5,[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías de oxigeno en 5°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,o2_a_10)
datetick('y')
axis([lon(1,1)-360 lon(64,1)-360 715885 731945])
caxis([-2.5 2.5])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,o2_a_10,[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías de oxigeno en 10°S','FontSize', 22)


figure()
pcolor(lon(:,1)-360,time,o2_a_15)
datetick('y')
axis([lon(1,1)-360 lon(66,1)-360 715885 731945])
caxis([-2.5 2.5])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,o2_a_15,[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías de oxigeno en 15°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,o2_a_25)
datetick('y')
axis([lon(1,1)-360 lon(71,1)-360 715885 731945])
caxis([-2.5 2.5])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,o2_a_25,[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías de oxigeno en 25°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,o2_a_35)
datetick('y')
axis([lon(1,1)-360 lon(69,1)-360 715885 731945])
caxis([-3 3])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,o2_a_35,[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías de oxigeno en 35°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,o2_a_45)
datetick('y')
axis([lon(1,1)-360 lon(68,1)-360 715885 731945])
caxis([-3 3])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,o2_a_45,[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías de oxigeno en 45°S','FontSize', 22)


% correlacion entre ONI y anomalias de oxigeno

ddd=size(o2_a_5);

for i=1:ddd(2)
    o2_a_5_2(:,i)=mediamovil(o2_a_5(:,i),3);
    o2_a_10_2(:,i)=mediamovil(o2_a_10(:,i),3);
    o2_a_15_2(:,i)=mediamovil(o2_a_15(:,i),3);
    o2_a_25_2(:,i)=mediamovil(o2_a_25(:,i),3);
    o2_a_35_2(:,i)=mediamovil(o2_a_35(:,i),3);
    o2_a_45_2(:,i)=mediamovil(o2_a_45(:,i),3);
    c5_oni(:,i)=corr(o2_a_5_2(:,i),ONI);
    c10_oni(:,i)=corr(o2_a_10_2(:,i),ONI);
    c15_oni(:,i)=corr(o2_a_15_2(:,i),ONI);
    c25_oni(:,i)=corr(o2_a_25_2(:,i),ONI);
    c35_oni(:,i)=corr(o2_a_35_2(:,i),ONI);
    c45_oni(:,i)=corr(o2_a_45_2(:,i),ONI);
    c5_sam(:,i)=corr(o2_a_5_2(:,i),sam);
    c10_sam(:,i)=corr(o2_a_10_2(:,i),sam);
    c15_sam(:,i)=corr(o2_a_15_2(:,i),sam);
    c25_sam(:,i)=corr(o2_a_25_2(:,i),sam);
    c35_sam(:,i)=corr(o2_a_35_2(:,i),sam);
    c45_sam(:,i)=corr(o2_a_45_2(:,i),sam);
end

figure()
plot(lon(:,1)-360,c5_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(61,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y Oxigeno a 5°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c5_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(61,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y Oxigeno a 5°S','FontSize', 24)



figure()
plot(lon(:,1)-360,c10_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(64,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y Oxigeno a 10°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c10_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(64,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y Oxigeno a 10°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c15_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(66,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y Oxigeno a 15°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c15_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(66,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y Oxigeno a 15°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c25_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(71,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y Oxigeno a 25°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c25_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(71,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y Oxigeno a 25°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c35_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(69,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y Oxigeno a 35°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c35_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(69,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y Oxigeno a 35°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c45_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(68,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y Oxigeno a 45°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c45_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(68,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y Oxigeno a 45°S','FontSize', 24)


%ANTICICLON
antic=load('ASPS.txt');
%solo el año que nos importa

antic=antic(1:34,:); %desde la columna 2 a la 13 son los meses y la 14 es el promedio anual

anticiclon=antic(:,2:13)';
anticiclon=anticiclon(:);

o2_a_5a=o2_a_5(121:end,:);
o2_a_10a=o2_a_10(121:end,:);
o2_a_15a=o2_a_15(121:end,:);
o2_a_25a=o2_a_25(121:end,:);
o2_a_35a=o2_a_35(121:end,:);
o2_a_45a=o2_a_45(121:end,:);


for i=1:ddd(2)
    o2_a_5_2a(:,i)=mediamovil(o2_a_5a(:,i),2);
    o2_a_10_2a(:,i)=mediamovil(o2_a_10a(:,i),2);
    o2_a_15_2a(:,i)=mediamovil(o2_a_15a(:,i),2);
    o2_a_25_2a(:,i)=mediamovil(o2_a_25a(:,i),2);
    o2_a_35_2a(:,i)=mediamovil(o2_a_35a(:,i),2);
    o2_a_45_2a(:,i)=mediamovil(o2_a_45a(:,i),2);
    c5_anti(:,i)=corr(o2_a_5_2a(:,i),anticiclon);
    c10_anti(:,i)=corr(o2_a_10_2a(:,i),anticiclon);
    c15_anti(:,i)=corr(o2_a_15_2a(:,i),anticiclon);
    c25_anti(:,i)=corr(o2_a_25_2a(:,i),anticiclon);
    c35_anti(:,i)=corr(o2_a_35_2a(:,i),anticiclon);
    c45_anti(:,i)=corr(o2_a_45_2a(:,i),anticiclon);
end

figure()
plot(lon(:,1)-360,c5_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(61,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y Oxigeno a 5°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c10_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(64,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y Oxigeno a 10°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c15_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(66,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y Oxigeno a 15°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c25_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(71,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y Oxigeno a 25°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c35_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(69,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y Oxigeno a 35°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c45_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(68,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y Oxigeno a 45°S','FontSize', 24)




%giff
timee=datevec(time);
timee=timee(:,1:2);

figure(1)
filename= 'ggifOf.gif' ;
set(gcf,'Renderer','zbuffer')
for i=1:D(3)
    pcolor(lon,lat,o2_af(:,:,i))
    colorbar
    caxis([-2.5 2.5])
    shading interp
    colormap(cmap)
    c = colorbar;
    c.Label.String = '[mmol/m^3]';
    hold on
    contour(lon,lat,o2_af(:,:,i),[0 0],'k','Linewidth',1)
    title(['anomalías mensuales de Oxigeno filtrado a 1100m en ZMO ' mat2str(timee(i,:))])
    drawnow
    frame=getframe(1);
    im=frame2im(frame);
    [imind,cm]=rgb2ind(im,256);
    if i==1
        imwrite(imind,cm,filename,'gif','LoopCount',inf)
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append')
    end
end



%% temperatura
clear all
clc

load 'ONI.mat'
load 'TEMP.mat'
load 'sam.mat'
time1=datevec(time);
D=size(TEMP);

%le quitamos el ciclo anual y estandarizamos 
TEMP_a=cicloanual_estan(TEMP,time);


%filtro de 11 meses usando el metodo coseno-lanczos
%primera pasada
VAR = TEMP_a; 
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

TEMP_af=VAR_ft1;


%tomamos una lat 
for i=1:D(3)
    TEMP_a_5=squeeze(TEMP_af(:,115,:)); %5
    TEMP_a_10=squeeze(TEMP_af(:,96,:)); %10
    TEMP_a_15=squeeze(TEMP_af(:,77,:)); %15
    TEMP_a_25=squeeze(TEMP_af(:,49,:)); %25
    TEMP_a_35=squeeze(TEMP_af(:,30,:)); %35
    TEMP_a_45=squeeze(TEMP_af(:,11,:)); %45
end


cmap = colormap_cpt('Optimus_Prime.cpt',265);



figure()
pcolor(lon(:,1)-360,time,TEMP_a_5')
datetick('y')
axis([lon(1,1)-360 lon(61,1)-360 715885 731945])
caxis([-2 2])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[°C]';
hold on
contour(lon(:,1)-360,time,TEMP_a_5',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías de temperatura a 1100m en 5°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,TEMP_a_10')
datetick('y')
axis([lon(1,1)-360 lon(64,1)-360 715885 731945])
caxis([-2 2])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[°C]';
hold on
contour(lon(:,1)-360,time,TEMP_a_10',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías de temperatura a 1100m en 10°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,TEMP_a_15')
datetick('y')
axis([lon(1,1)-360 lon(66,1)-360 715885 731945])
caxis([-2.5 2.5])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[°C]';
hold on
contour(lon(:,1)-360,time,TEMP_a_15',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías de temperatura a 1100m en 15°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,TEMP_a_25')
datetick('y')
axis([lon(1,1)-360 lon(71,1)-360 715885 731945])
caxis([-2.5 2.5])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[°C]';
hold on
contour(lon(:,1)-360,time,TEMP_a_25',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías de temperatura a 1100m en 25°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,TEMP_a_35')
datetick('y')
axis([lon(1,1)-360 lon(69,1)-360 715885 731945])
caxis([-3 3])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[°C]';
hold on
contour(lon(:,1)-360,time,TEMP_a_35',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías de temperatura a 1100m en 35°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,TEMP_a_45')
datetick('y')
axis([lon(1,1)-360 lon(68,1)-360 715885 731945])
caxis([-3 3])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[°C]';
hold on
contour(lon(:,1)-360,time,TEMP_a_45',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías de temperatura a 1100m en 45°S','FontSize', 22)


% correlacion entre ONI y anomalias de oxigeno

TEMP_a_5=TEMP_a_5';
TEMP_a_10=TEMP_a_10';
TEMP_a_15=TEMP_a_15';
TEMP_a_25=TEMP_a_25';
TEMP_a_35=TEMP_a_35';
TEMP_a_45=TEMP_a_45';
ddd=size(TEMP_a_5);

for i=1:ddd(2)
    TEMP_a_5_2(:,i)=mediamovil(TEMP_a_5(:,i),3);
    TEMP_a_10_2(:,i)=mediamovil(TEMP_a_10(:,i),3);
    TEMP_a_15_2(:,i)=mediamovil(TEMP_a_15(:,i),3);
    TEMP_a_25_2(:,i)=mediamovil(TEMP_a_25(:,i),3);
    TEMP_a_35_2(:,i)=mediamovil(TEMP_a_35(:,i),3);
    TEMP_a_45_2(:,i)=mediamovil(TEMP_a_45(:,i),3);
    c5_oni(:,i)=corr(TEMP_a_5_2(:,i),ONI);
    c10_oni(:,i)=corr(TEMP_a_10_2(:,i),ONI);
    c15_oni(:,i)=corr(TEMP_a_15_2(:,i),ONI);
    c25_oni(:,i)=corr(TEMP_a_25_2(:,i),ONI);
    c35_oni(:,i)=corr(TEMP_a_35_2(:,i),ONI);
    c45_oni(:,i)=corr(TEMP_a_45_2(:,i),ONI);
    c5_sam(:,i)=corr(TEMP_a_5_2(:,i),sam);
    c10_sam(:,i)=corr(TEMP_a_10_2(:,i),sam);
    c15_sam(:,i)=corr(TEMP_a_15_2(:,i),sam);
    c25_sam(:,i)=corr(TEMP_a_25_2(:,i),sam);
    c35_sam(:,i)=corr(TEMP_a_35_2(:,i),sam);
    c45_sam(:,i)=corr(TEMP_a_45_2(:,i),sam);
end

figure()
plot(lon(:,1)-360,c5_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(61,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y temperatura a 1100m 5°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c5_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(61,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y temperatura a 1100m 5°S','FontSize', 24)


figure()
plot(lon(:,1)-360,c10_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(64,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y temperatura a 1100m 10°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c10_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(64,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y temperatura a 1100m 10°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c15_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(66,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y temperatura a 1100m 15°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c15_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(66,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y temperatura a 1100m 15°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c25_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(71,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y temperatura a 1100m 25°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c25_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(71,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y temperatura a 1100m 25°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c35_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(69,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y temperatura a 1100m 35°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c35_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(69,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y temperatura a 1100m 35°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c45_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(68,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y temperatura a 1100m 45°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c45_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(68,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y temperatura a 1100m 45°S','FontSize', 24)

%ANTICICLON
antic=load('ASPS.txt');
%solo el año que nos importa

antic=antic(1:34,:); %desde la columna 2 a la 13 son los meses y la 14 es el promedio anual

anticiclon=antic(:,2:13)';
anticiclon=anticiclon(:);


TEMP_a_5a=TEMP_a_5(121:end,:);
TEMP_a_10a=TEMP_a_10(121:end,:);
TEMP_a_15a=TEMP_a_15(121:end,:);
TEMP_a_25a=TEMP_a_25(121:end,:);
TEMP_a_35a=TEMP_a_35(121:end,:);
TEMP_a_45a=TEMP_a_45(121:end,:);

for i=1:ddd(2)
    TEMP_a_5_2a(:,i)=mediamovil(TEMP_a_5a(:,i),2);
    TEMP_a_10_2a(:,i)=mediamovil(TEMP_a_10a(:,i),2);
    TEMP_a_15_2a(:,i)=mediamovil(TEMP_a_15a(:,i),2);
    TEMP_a_25_2a(:,i)=mediamovil(TEMP_a_25a(:,i),2);
    TEMP_a_35_2a(:,i)=mediamovil(TEMP_a_35a(:,i),2);
    TEMP_a_45_2a(:,i)=mediamovil(TEMP_a_45a(:,i),2);
    c5_anti(:,i)=corr(TEMP_a_5_2a(:,i),anticiclon);
    c10_anti(:,i)=corr(TEMP_a_10_2a(:,i),anticiclon);
    c15_anti(:,i)=corr(TEMP_a_15_2a(:,i),anticiclon);
    c25_anti(:,i)=corr(TEMP_a_25_2a(:,i),anticiclon);
    c35_anti(:,i)=corr(TEMP_a_35_2a(:,i),anticiclon);
    c45_anti(:,i)=corr(TEMP_a_45_2a(:,i),anticiclon);
end

figure()
plot(lon(:,1)-360,c5_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(61,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y temperatura a 1100m 5°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c10_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(64,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y temperatura a 1100m 10°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c15_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(66,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y temperatura a 1100m 15°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c25_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(71,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y temperatura a 1100m 25°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c35_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(69,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y temperatura a 1100m 35°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c45_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(68,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y temperatura a 1100m 45°S','FontSize', 24)



%giff
timee=datevec(time);
timee=timee(:,1:2);

figure(1)
filename= 'ggifTf.gif' ;
set(gcf,'Renderer','zbuffer')
for i=1:D(3)
    pcolor(lon,lat,TEMP_af(:,:,i))
    caxis([-2.5 2.5])
    shading interp
    colormap(cmap)
    c = colorbar;
    c.Label.String = '[°C]';
    hold on
    contour(lon,lat,TEMP_af(:,:,i),[0 0],'k','Linewidth',1)
    title(['anomalías mensuales de temperatura filtrada a 1100m en ZMO ' mat2str(timee(i,:))])
    drawnow
    frame=getframe(1);
    im=frame2im(frame);
    [imind,cm]=rgb2ind(im,256);
    if i==1
        imwrite(imind,cm,filename,'gif','LoopCount',inf)
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append')
    end
end


%% SST filtrado
clear all
clc

load 'ONI.mat'
load 'SST1.mat'
load 'sam.mat'
time1=datevec(time);
D=size(SST1);



%le quitamos el ciclo anual y estandarizamos 
SST_a=cicloanual_estan(SST1,time);

%filtro de 11 meses usando el metodo coseno-lanczos
%primera pasada
VAR = SST_a; 
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

SST_af=VAR_ft1;



%tomamos una lat 
for i=1:D(3)
    SST_a_5=squeeze(SST_af(:,115,:)); %5
    SST_a_10=squeeze(SST_af(:,96,:)); %10
    SST_a_15=squeeze(SST_af(:,77,:)); %15
    SST_a_25=squeeze(SST_af(:,49,:)); %25
    SST_a_35=squeeze(SST_af(:,30,:)); %35
    SST_a_45=squeeze(SST_af(:,11,:)); %45
end
cmap = colormap_cpt('Optimus_Prime.cpt',265);

figure()
pcolor(lon(:,1)-360,time,SST_a_5')
datetick('y')
axis([lon(1,1)-360 lon(61,1)-360 715885 731945])
caxis([-2 2])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[°C]';
hold on
contour(lon(:,1)-360,time,SST_a_5',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías de SST en 5°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,SST_a_10')
datetick('y')
axis([lon(1,1)-360 lon(64,1)-360 715885 731945])
caxis([-2 2])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[°C]';
hold on
contour(lon(:,1)-360,time,SST_a_10',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías de SST en 10°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,SST_a_15')
datetick('y')
axis([lon(1,1)-360 lon(66,1)-360 715885 731945])
caxis([-2 2])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[°C]';
hold on
contour(lon(:,1)-360,time,SST_a_15',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías de SST en 15°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,SST_a_25')
datetick('y')
axis([lon(1,1)-360 lon(71,1)-360 715885 731945])
caxis([-1.5 1.5])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[°C]';
hold on
contour(lon(:,1)-360,time,SST_a_25',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías de SST en 25°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,SST_a_35')
datetick('y')
axis([lon(1,1)-360 lon(69,1)-360 715885 731945])
caxis([-1 1])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[°C]';
hold on
contour(lon(:,1)-360,time,SST_a_35',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías de SST en 35°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,SST_a_45')
datetick('y')
axis([lon(1,1)-360 lon(68,1)-360 715885 731945])
caxis([-1 1])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[°C]';
hold on
contour(lon(:,1)-360,time,SST_a_45',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías de SST en 45°S','FontSize', 22)

% correlacion entre ONI y sam y anomalias de sst

SST_a_5=SST_a_5';
SST_a_10=SST_a_10';
SST_a_15=SST_a_15';
SST_a_25=SST_a_25';
SST_a_35=SST_a_35';
SST_a_45=SST_a_45';
ddd=size(SST_a_5);

for i=1:ddd(2)
    SST_a_5_2(:,i)=mediamovil(SST_a_5(:,i),3);
    SST_a_10_2(:,i)=mediamovil(SST_a_10(:,i),3);
    SST_a_15_2(:,i)=mediamovil(SST_a_15(:,i),3);
    SST_a_25_2(:,i)=mediamovil(SST_a_25(:,i),3);
    SST_a_35_2(:,i)=mediamovil(SST_a_35(:,i),3);
    SST_a_45_2(:,i)=mediamovil(SST_a_45(:,i),3);
    c5_oni(:,i)=corr(SST_a_5_2(:,i),ONI);
    c10_oni(:,i)=corr(SST_a_10_2(:,i),ONI);
    c15_oni(:,i)=corr(SST_a_15_2(:,i),ONI);
    c25_oni(:,i)=corr(SST_a_25_2(:,i),ONI);
    c35_oni(:,i)=corr(SST_a_35_2(:,i),ONI);
    c45_oni(:,i)=corr(SST_a_45_2(:,i),ONI);
    c5_sam(:,i)=corr(SST_a_5_2(:,i),sam);
    c10_sam(:,i)=corr(SST_a_10_2(:,i),sam);
    c15_sam(:,i)=corr(SST_a_15_2(:,i),sam);
    c25_sam(:,i)=corr(SST_a_25_2(:,i),sam);
    c35_sam(:,i)=corr(SST_a_35_2(:,i),sam);
    c45_sam(:,i)=corr(SST_a_45_2(:,i),sam);
end

figure()
plot(lon(:,1)-360,c5_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(61,1)-360 0.3 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y SST a 5°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c5_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(61,1)-360 -1 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y SST a 5°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c10_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(64,1)-360 0.3 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y SST a 10°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c10_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(64,1)-360 -1 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y SST a 10°S','FontSize', 24)


figure()
plot(lon(:,1)-360,c15_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(66,1)-360 0.3 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y SST a 15°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c15_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(66,1)-360 -1 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y SST a 15°S','FontSize', 24)


figure()
plot(lon(:,1)-360,c25_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(71,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y SST a 25°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c25_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(71,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y SST a 25°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c35_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(69,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y SST a 35°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c35_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(69,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y SST a 35°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c45_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(68,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y SST a 45°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c45_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(68,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y SST a 45°S','FontSize', 24)


%ANTICICLON

antic=load('ASPS.txt');
%solo el año que nos importa

antic=antic(1:34,:); %desde la columna 2 a la 13 son los meses y la 14 es el promedio anual

anticiclon=antic(:,2:13)';
anticiclon=anticiclon(:);

SST_a_5a=SST_a_5(121:end,:);
SST_a_10a=SST_a_10(121:end,:);
SST_a_15a=SST_a_15(121:end,:);
SST_a_25a=SST_a_25(121:end,:);
SST_a_35a=SST_a_35(121:end,:);
SST_a_45a=SST_a_45(121:end,:);

for i=1:ddd(2)
    SST_a_5_2a(:,i)=mediamovil(SST_a_5a(:,i),2);
    SST_a_10_2a(:,i)=mediamovil(SST_a_10a(:,i),2);
    SST_a_15_2a(:,i)=mediamovil(SST_a_15a(:,i),2);
    SST_a_25_2a(:,i)=mediamovil(SST_a_25a(:,i),2);
    SST_a_35_2a(:,i)=mediamovil(SST_a_35a(:,i),2);
    SST_a_45_2a(:,i)=mediamovil(SST_a_45a(:,i),2);
    c5_anti(:,i)=corr(SST_a_5_2a(:,i),anticiclon);
    c10_anti(:,i)=corr(SST_a_10_2a(:,i),anticiclon);
    c15_anti(:,i)=corr(SST_a_15_2a(:,i),anticiclon);
    c25_anti(:,i)=corr(SST_a_25_2a(:,i),anticiclon);
    c35_anti(:,i)=corr(SST_a_35_2a(:,i),anticiclon);
    c45_anti(:,i)=corr(SST_a_45_2a(:,i),anticiclon);
end

figure()
plot(lon(:,1)-360,c5_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(61,1)-360 -1 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y SST a 5°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c10_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(64,1)-360 -1 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y SST a 10°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c15_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(66,1)-360 -1 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y SST a 15°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c25_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(71,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y SST a 25°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c35_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(69,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y SST a 35°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c45_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(68,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y SST a 45°S','FontSize', 24)



%giff
timee=datevec(time);
timee=timee(:,1:2);

figure(1)
filename= 'ggifSSTf.gif' ;
set(gcf,'Renderer','zbuffer')
for i=1:D(3)
    pcolor(lon,lat,SST_af(:,:,i))
    colorbar
    caxis([-1.5 1.5])
    shading interp
    colormap(cmap)
    c = colorbar;
    c.Label.String = '[°C]';
    hold on
    contour(lon,lat,SST_af(:,:,i),[0 0],'k','Linewidth',1)
    title(['anomalías mensuales de SST filtrado ' mat2str(timee(i,:))])
    drawnow
    frame=getframe(1);
    im=frame2im(frame);
    [imind,cm]=rgb2ind(im,256);
    if i==1
        imwrite(imind,cm,filename,'gif','LoopCount',inf)
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append')
    end
end



%% niVEL DEL MAR SSH

clear all
clc

load 'ONI.mat'
load 'SSH.mat'
load 'sam.mat'
time1=datevec(time);
D=size(SSH);


%le quitamos el ciclo anual y estandarizamos 
SSH_a=cicloanual_estan(SSH,time);


%filtro de 11 meses usando el metodo coseno-lanczos
%primera pasada
VAR = SSH_a; 
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

SSH_af=VAR_ft1;


%tomamos una lat 
for i=1:D(3)
    SSH_a_5=squeeze(SSH_af(:,115,:)); %5
    SSH_a_10=squeeze(SSH_af(:,96,:)); %10
    SSH_a_15=squeeze(SSH_af(:,77,:)); %15
    SSH_a_25=squeeze(SSH_af(:,49,:)); %25
    SSH_a_35=squeeze(SSH_af(:,30,:)); %35
    SSH_a_45=squeeze(SSH_af(:,11,:)); %45
end
cmap = colormap_cpt('Optimus_Prime.cpt',265);


figure()
pcolor(lon(:,1)-360,time,SSH_a_5')
datetick('y')
axis([lon(1,1)-360 lon(61,1)-360 715885 731945])
caxis([-3 3])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[°cm]';
hold on
contour(lon(:,1)-360,time,SSH_a_5',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías del nivel del mar en 5°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,SSH_a_10')
datetick('y')
axis([lon(1,1)-360 lon(64,1)-360 715885 731945])
caxis([-3 3])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[°cm]';
hold on
contour(lon(:,1)-360,time,SSH_a_10',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías del nivel del mar en 10°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,SSH_a_15')
datetick('y')
axis([lon(1,1)-360 lon(66,1)-360 715885 731945])
caxis([-3 3])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[°cm]';
hold on
contour(lon(:,1)-360,time,SSH_a_15',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías del nivel del mar en 15°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,SSH_a_25')
datetick('y')
axis([lon(1,1)-360 lon(71,1)-360 715885 731945])
caxis([-3 3])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[°cm]';
hold on
contour(lon(:,1)-360,time,SSH_a_25',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías del nivel del mar en 25°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,SSH_a_35')
datetick('y')
axis([lon(1,1)-360 lon(69,1)-360 715885 731945])
caxis([-2.5 2.5])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[°cm]';
hold on
contour(lon(:,1)-360,time,SSH_a_35',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías del nivel del mar en 35°S','FontSize', 22)


figure()
pcolor(lon(:,1)-360,time,SSH_a_45')
datetick('y')
axis([lon(1,1)-360 lon(68,1)-360 715885 731945])
caxis([-2.5 2.5])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[°cm]';
hold on
contour(lon(:,1)-360,time,SSH_a_45',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Hovmoller de anomalías del nivel del mar en 45°S','FontSize', 22)



% correlacion entre ONI y sam y anomalias de SSH

SSH_a_5=SSH_a_5';
SSH_a_10=SSH_a_10';
SSH_a_15=SSH_a_15';
SSH_a_25=SSH_a_25';
SSH_a_35=SSH_a_35';
SSH_a_45=SSH_a_45';
ddd=size(SSH_a_5);

for i=1:ddd(2)
    SSH_a_5_2(:,i)=mediamovil(SSH_a_5(:,i),3);
    SSH_a_10_2(:,i)=mediamovil(SSH_a_10(:,i),3);
    SSH_a_15_2(:,i)=mediamovil(SSH_a_15(:,i),3);
    SSH_a_25_2(:,i)=mediamovil(SSH_a_25(:,i),3);
    SSH_a_35_2(:,i)=mediamovil(SSH_a_35(:,i),3);
    SSH_a_45_2(:,i)=mediamovil(SSH_a_45(:,i),3);
    c5_oni(:,i)=corr(SSH_a_5_2(:,i),ONI);
    c10_oni(:,i)=corr(SSH_a_10_2(:,i),ONI);
    c15_oni(:,i)=corr(SSH_a_15_2(:,i),ONI);
    c25_oni(:,i)=corr(SSH_a_25_2(:,i),ONI);
    c35_oni(:,i)=corr(SSH_a_35_2(:,i),ONI);
    c45_oni(:,i)=corr(SSH_a_45_2(:,i),ONI);
    c5_sam(:,i)=corr(SSH_a_5_2(:,i),sam);
    c10_sam(:,i)=corr(SSH_a_10_2(:,i),sam);
    c15_sam(:,i)=corr(SSH_a_15_2(:,i),sam);
    c25_sam(:,i)=corr(SSH_a_25_2(:,i),sam);
    c35_sam(:,i)=corr(SSH_a_35_2(:,i),sam);
    c45_sam(:,i)=corr(SSH_a_45_2(:,i),sam);
end

figure()
plot(lon(:,1)-360,c5_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(61,1)-360 -0.4 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y nivel del mar a 5°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c5_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(61,1)-360 -0.4 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y nivel del mar a 5°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c10_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(64,1)-360 -0.4 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y nivel del mar a 10°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c10_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(64,1)-360 -0.4 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y nivel del mar a 10°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c15_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(66,1)-360 -0.4 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y nivel del mar a 15°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c15_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(66,1)-360 -0.4 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y nivel del mar a 15°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c25_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(71,1)-360 -0.55 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y nivel del mar a 25°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c25_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(71,1)-360 -0.55 0.5])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y nivel del mar a 25°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c35_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(69,1)-360 -0.55 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y nivel del mar a 35°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c35_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(69,1)-360 -0.55 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y nivel del mar a 35°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c45_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(68,1)-360 -0.55 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y nivel del mar a 45°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c45_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(68,1)-360 -0.55 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y nivel del mar a 45°S','FontSize', 24)


%ANTICICLON

antic=load('ASPS.txt');
%solo el año que nos importa

antic=antic(1:34,:); %desde la columna 2 a la 13 son los meses y la 14 es el promedio anual

anticiclon=antic(:,2:13)';
anticiclon=anticiclon(:);


SSH_a_5a=SSH_a_5(121:end,:);
SSH_a_10a=SSH_a_10(121:end,:);
SSH_a_15a=SSH_a_15(121:end,:);
SSH_a_25a=SSH_a_25(121:end,:);
SSH_a_35a=SSH_a_35(121:end,:);
SSH_a_45a=SSH_a_45(121:end,:);

for i=1:ddd(2)
    SSH_a_5_2a(:,i)=mediamovil(SSH_a_5a(:,i),2);
    SSH_a_10_2a(:,i)=mediamovil(SSH_a_10a(:,i),2);
    SSH_a_15_2a(:,i)=mediamovil(SSH_a_15a(:,i),2);
    SSH_a_25_2a(:,i)=mediamovil(SSH_a_25a(:,i),2);
    SSH_a_35_2a(:,i)=mediamovil(SSH_a_35a(:,i),2);
    SSH_a_45_2a(:,i)=mediamovil(SSH_a_45a(:,i),2);
    c5_anti(:,i)=corr(SSH_a_5_2a(:,i),anticiclon);
    c10_anti(:,i)=corr(SSH_a_10_2a(:,i),anticiclon);
    c15_anti(:,i)=corr(SSH_a_15_2a(:,i),anticiclon);
    c25_anti(:,i)=corr(SSH_a_25_2a(:,i),anticiclon);
    c35_anti(:,i)=corr(SSH_a_35_2a(:,i),anticiclon);
    c45_anti(:,i)=corr(SSH_a_45_2a(:,i),anticiclon);
end

figure()
plot(lon(:,1)-360,c5_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(61,1)-360 -0.4 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y nivel del mar a 5°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c10_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(64,1)-360 -0.4 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y nivel del mar a 10°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c15_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(66,1)-360 -0.4 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y nivel del mar a 15°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c25_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(71,1)-360 -0.55 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y nivel del mar a 25°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c35_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(69,1)-360 -0.55 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y nivel del mar a 35°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c45_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(68,1)-360 -0.55 1])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y nivel del mar a 45°S','FontSize', 24)


%giff
timee=datevec(time);
timee=timee(:,1:2);


figure(1)
filename= 'ggifSSHf.gif' ;
set(gcf,'Renderer','zbuffer')
for i=1:D(3)
    pcolor(lon,lat,SSH_af(:,:,i))
    colorbar
    caxis([-2.5 2.5])
    shading interp
    colormap(cmap)
    c = colorbar;
    c.Label.String = '[°cm]';
    hold on
    contour(lon,lat,SSH_af(:,:,i),[0 0],'k','Linewidth',1)
    title(['anomalías mensuales del nivel del mar filtrado' mat2str(timee(i,:))])
    drawnow
    frame=getframe(1);
    im=frame2im(frame);
    [imind,cm]=rgb2ind(im,256);
    if i==1
        imwrite(imind,cm,filename,'gif','LoopCount',inf)
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append')
    end
end

