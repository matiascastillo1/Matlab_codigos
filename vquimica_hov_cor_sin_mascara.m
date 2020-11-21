%% HOVMOLLERS de DIC, NH4 y NO3 .Y CORRELACION CON ESTOS Y EL ONI
%% PERO SIN LA MASCARA, ES DECIR VAN A TENER TODAS LAS LONGITUDES EN LATITUDES MAS ALTAS
%% ADEMAS SERAN OTRAS 3 LATITUDES RESPECTO A EL AVANCE ANTERIOR ES DECIR 5, 10, 15, 25, 35 Y 40 °S
%% y de las variables filtradas de variabilidad interanual

clear all
close all
clc

%% DIC
load 'SALT.mat'
load 'ONI.mat'
load 'DIC.mat'
load 'sam.mat'
time1=datevec(time);
D=size(DIC);

% Normalizamos el DIC con respecto a la salinidad usando la salanidad de
% referencia de 35 PSU, la formula se extrajo del paper (Friss. K, et al.
% 2003)
for i=1:D(1)
    for j=1:D(2)
        for k=1:D(3)
            DIC_n(i,j,k)=(DIC(i,j,k)*35)/SALT(i,j,k);
        end
    end
end


%le quitamos el ciclo anual y estandarizamos 
DIC_a=cicloanual_estan(DIC_n,time);

%filtro de 11 meses usando el metodo coseno-lanczos
%primera pasada
VAR = DIC_a; 
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

DIC_af=VAR_ft1;


%tomamos una lat 
DIC_a_5=squeeze(DIC_af(:,115,:)); %5
DIC_a_10=squeeze(DIC_af(:,96,:)); %10
DIC_a_15=squeeze(DIC_af(:,77,:)); %15
DIC_a_25=squeeze(DIC_af(:,49,:)); %25
DIC_a_35=squeeze(DIC_af(:,30,:)); %35
DIC_a_45=squeeze(DIC_af(:,11,:)); %45

cmap = colormap_cpt('Optimus_Prime.cpt',265);

figure()
pcolor(lon(:,1)-360,time,DIC_a_5')
datetick('y')
axis([lon(1,1)-360 lon(61,1)-360 715885 731945])
caxis([-2.5 2.5])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,DIC_a_5',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías de DIC en 5°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,DIC_a_10')
datetick('y')
axis([lon(1,1)-360 lon(64,1)-360 715885 731945])
caxis([-2.5 2.5])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,DIC_a_10',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías de DIC en 10°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,DIC_a_15')
datetick('y')
axis([lon(1,1)-360 lon(66,1)-360 715885 731945])
caxis([-3 3])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,DIC_a_15',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías de DIC en 15°S','FontSize', 22)


figure()
pcolor(lon(:,1)-360,time,DIC_a_25')
datetick('y')
axis([lon(1,1)-360 lon(71,1)-360 715885 731945])
caxis([-2.5 2.5])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,DIC_a_25',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías de DIC en 25°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,DIC_a_35')
datetick('y')
axis([lon(1,1)-360 lon(69,1)-360 715885 731945])
caxis([-2.5 2.5])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,DIC_a_35',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías de DIC en 35°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,DIC_a_45')
datetick('y')
axis([lon(1,1)-360 lon(68,1)-360 715885 731945])
caxis([-2.5 2.5])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,DIC_a_45',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías de DIC en 45°S','FontSize', 22)

% correlacion entre ONI y sam y anomalias de DIC

DIC_a_5=DIC_a_5';
DIC_a_10=DIC_a_10';
DIC_a_15=DIC_a_15';
DIC_a_25=DIC_a_25';
DIC_a_35=DIC_a_35';
DIC_a_45=DIC_a_45';
ddd=size(DIC_a_5);

for i=1:ddd(2)
    DIC_a_5_2(:,i)=mediamovil(DIC_a_5(:,i),3);
    DIC_a_10_2(:,i)=mediamovil(DIC_a_10(:,i),3);
    DIC_a_15_2(:,i)=mediamovil(DIC_a_15(:,i),3);
    DIC_a_25_2(:,i)=mediamovil(DIC_a_25(:,i),3);
    DIC_a_35_2(:,i)=mediamovil(DIC_a_35(:,i),3);
    DIC_a_45_2(:,i)=mediamovil(DIC_a_45(:,i),3);
    c5_oni(:,i)=corr(DIC_a_5_2(:,i),ONI);
    c5_sam(:,i)=corr(DIC_a_5_2(:,i),sam);
    c10_oni(:,i)=corr(DIC_a_10_2(:,i),ONI);
    c10_sam(:,i)=corr(DIC_a_10_2(:,i),sam);
    c15_oni(:,i)=corr(DIC_a_15_2(:,i),ONI);
    c15_sam(:,i)=corr(DIC_a_15_2(:,i),sam);
    c25_oni(:,i)=corr(DIC_a_25_2(:,i),ONI);
    c25_sam(:,i)=corr(DIC_a_25_2(:,i),sam);
    c35_oni(:,i)=corr(DIC_a_35_2(:,i),ONI);
    c35_sam(:,i)=corr(DIC_a_35_2(:,i),sam);
    c45_oni(:,i)=corr(DIC_a_45_2(:,i),ONI);
    c45_sam(:,i)=corr(DIC_a_45_2(:,i),sam);
end

figure()
plot(lon(:,1)-360,c5_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(61,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y DIC a 5°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c5_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(61,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y DIC a 5°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c10_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(64,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y DIC a 10°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c10_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(64,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y DIC a 10°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c15_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(66,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y DIC a 15°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c15_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(66,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y DIC a 15°S','FontSize', 24)


figure()
plot(lon(:,1)-360,c25_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(71,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y DIC a 25°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c25_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(71,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y DIC a 25°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c35_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(69,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y DIC a 35°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c35_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(69,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y DIC a 35°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c45_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(68,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y DIC a 45°S','FontSize', 24)


figure()
plot(lon(:,1)-360,c45_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(68,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y DIC a 45°S','FontSize', 24)



%cor anticiclon
antic=load('ASPS.txt');
%solo el año que nos importa

antic=antic(1:34,:); %desde la columna 2 a la 13 son los meses y la 14 es el promedio anual

anticiclon=antic(:,2:13)';
anticiclon=anticiclon(:);

DIC_a_5a=DIC_a_5(121:end,:);
DIC_a_10a=DIC_a_10(121:end,:);
DIC_a_15a=DIC_a_15(121:end,:);
DIC_a_25a=DIC_a_25(121:end,:);
DIC_a_35a=DIC_a_35(121:end,:);
DIC_a_45a=DIC_a_45(121:end,:);

for i=1:ddd(2)
    DIC_a_5_2a(:,i)=mediamovil(DIC_a_5a(:,i),2);
    DIC_a_10_2a(:,i)=mediamovil(DIC_a_10a(:,i),2);
    DIC_a_15_2a(:,i)=mediamovil(DIC_a_15a(:,i),2);
    DIC_a_25_2a(:,i)=mediamovil(DIC_a_25a(:,i),2);
    DIC_a_35_2a(:,i)=mediamovil(DIC_a_35a(:,i),2);
    DIC_a_45_2a(:,i)=mediamovil(DIC_a_45a(:,i),2);
    c5_anti(:,i)=corr(DIC_a_5_2a(:,i),anticiclon);
    c10_anti(:,i)=corr(DIC_a_10_2a(:,i),anticiclon);
    c15_anti(:,i)=corr(DIC_a_15_2a(:,i),anticiclon);
    c25_anti(:,i)=corr(DIC_a_25_2a(:,i),anticiclon);
    c35_anti(:,i)=corr(DIC_a_35_2a(:,i),anticiclon);
    c45_anti(:,i)=corr(DIC_a_45_2a(:,i),anticiclon);
end



figure()
plot(lon(:,1)-360,c5_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(61,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y DIC a 5°S','FontSize', 24)


figure()
plot(lon(:,1)-360,c10_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(64,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y DIC a 10°S','FontSize', 24)


figure()
plot(lon(:,1)-360,c15_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(66,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y DIC a 15°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c25_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(71,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y DIC a 25°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c35_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(69,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y DIC a 35°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c45_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(68,1)-360 -0.75 0.75])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y DIC a 45°S','FontSize', 24)


%GIFF
timee=datevec(time);
timee=timee(:,1:2);

figure(1)
filename= 'ggifDICf.gif' ;
set(gcf,'Renderer','zbuffer')
for i=1:D(3)
    pcolor(lon,lat,DIC_af(:,:,i))
    colorbar
    caxis([-2.5 2.5])
    shading interp
    colormap(cmap)
    c = colorbar;
    c.Label.String = '[mmol/m^3]';
    hold on
    contour(lon,lat,DIC_af(:,:,i),[0 0],'k','Linewidth',1)
    title(['anomalías mensuales de DIC filtrado a 1100m en ZMO ' mat2str(timee(i,:))])
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

%% NH4 Filtrado
clear all
clc 

load 'ONI.mat'
load 'NH4.mat'
load 'sam.mat'
time1=datevec(time);
D=size(NH4);

%le quitamos el ciclo anual y estandarizamos 
NH4_a=cicloanual_estan(NH4,time);


%filtro de 11 meses usando el metodo coseno-lanczos
%primera pasada
VAR = NH4_a; 
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

NH4_af=VAR_ft1;

%tomamos una lat 
for i=1:D(3)
    NH4_a_5=squeeze(NH4_af(:,115,:)); %5
    NH4_a_10=squeeze(NH4_af(:,96,:)); %10
    NH4_a_15=squeeze(NH4_af(:,77,:)); %15
    NH4_a_25=squeeze(NH4_af(:,49,:)); %25
    NH4_a_35=squeeze(NH4_af(:,30,:)); %35
    NH4_a_45=squeeze(NH4_af(:,11,:)); %45
end

cmap = colormap_cpt('Optimus_Prime.cpt',265);


figure()
pcolor(lon(:,1)-360,time,NH4_a_5')
datetick('y')
axis([lon(1,1)-360 lon(61,1)-360 715885 731945])
caxis([-2 2])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,NH4_a_5',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías de NH4 en 5°S','FontSize', 22)


figure()
pcolor(lon(:,1)-360,time,NH4_a_10')
datetick('y')
axis([lon(1,1)-360 lon(64,1)-360 715885 731945])
caxis([-2 2])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,NH4_a_10',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías de NH4 en 10°S','FontSize', 22)


figure()
pcolor(lon(:,1)-360,time,NH4_a_15')
datetick('y')
axis([lon(1,1)-360 lon(66,1)-360 715885 731945])
caxis([-2.5 2.5])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,NH4_a_15',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías de NH4 en 15°S','FontSize', 22)



figure()
pcolor(lon(:,1)-360,time,NH4_a_25')
datetick('y')
axis([lon(1,1)-360 lon(71,1)-360 715885 731945])
caxis([-1.5 1.5])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,NH4_a_25',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías de NH4 en 25°S','FontSize', 22)


figure()
pcolor(lon(:,1)-360,time,NH4_a_35')
datetick('y')
axis([lon(1,1)-360 lon(69,1)-360 715885 731945])
caxis([-1.5 1.5])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,NH4_a_35',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías de NH4 en 35°S','FontSize', 22)



figure()
pcolor(lon(:,1)-360,time,NH4_a_45')
datetick('y')
axis([lon(1,1)-360 lon(68,1)-360 715885 731945])
caxis([-1 1])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,NH4_a_45',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías de NH4 en 45°S','FontSize', 22)

% correlacion entre ONI y SAM y anomalias de amonio 

NH4_a_5=NH4_a_5';
NH4_a_10=NH4_a_10';
NH4_a_15=NH4_a_15';
NH4_a_25=NH4_a_25';
NH4_a_35=NH4_a_35';
NH4_a_45=NH4_a_45';
ddd=size(NH4_a_5);

for i=1:ddd(2)
    NH4_a_5_2(:,i)=mediamovil(NH4_a_5(:,i),3);
    NH4_a_10_2(:,i)=mediamovil(NH4_a_10(:,i),3);
    NH4_a_15_2(:,i)=mediamovil(NH4_a_15(:,i),3);
    NH4_a_25_2(:,i)=mediamovil(NH4_a_25(:,i),3);
    NH4_a_35_2(:,i)=mediamovil(NH4_a_35(:,i),3);
    NH4_a_45_2(:,i)=mediamovil(NH4_a_45(:,i),3);
    c5_oni(:,i)=corr(NH4_a_5_2(:,i),ONI);
    c10_oni(:,i)=corr(NH4_a_10_2(:,i),ONI);
    c15_oni(:,i)=corr(NH4_a_15_2(:,i),ONI);
    c25_oni(:,i)=corr(NH4_a_25_2(:,i),ONI);
    c35_oni(:,i)=corr(NH4_a_35_2(:,i),ONI);
    c45_oni(:,i)=corr(NH4_a_45_2(:,i),ONI);
    c5_sam(:,i)=corr(NH4_a_5_2(:,i),sam);
    c10_sam(:,i)=corr(NH4_a_10_2(:,i),sam);
    c15_sam(:,i)=corr(NH4_a_15_2(:,i),sam);
    c25_sam(:,i)=corr(NH4_a_25_2(:,i),sam);
    c35_sam(:,i)=corr(NH4_a_35_2(:,i),sam);
    c45_sam(:,i)=corr(NH4_a_45_2(:,i),sam);
end

figure()
plot(lon(:,1)-360,c5_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(61,1)-360 -1 0.3])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y NH4 a 5°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c5_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(61,1)-360 -1 0.3])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y NH4 a 5°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c10_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(64,1)-360 -1 0.3])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y NH4 a 10°S','FontSize', 24)


figure()
plot(lon(:,1)-360,c10_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(64,1)-360 -1 0.3])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y NH4 a 10°S','FontSize', 24)



figure()
plot(lon(:,1)-360,c15_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(66,1)-360 -1 0.3])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y NH4 a 15°S','FontSize', 24)


figure()
plot(lon(:,1)-360,c15_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(66,1)-360 -1 0.3])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y NH4 a 15°S','FontSize', 24)



figure()
plot(lon(:,1)-360,c25_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(71,1)-360 -0.6 0.45])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y NH4 a 25°S','FontSize', 24)


figure()
plot(lon(:,1)-360,c25_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(71,1)-360 -0.6 0.45])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y NH4 a 25°S','FontSize', 24)



figure()
plot(lon(:,1)-360,c35_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(69,1)-360 -0.6 0.45])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y NH4 a 35°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c35_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(69,1)-360 -0.6 0.45])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y NH4 a 35°S','FontSize', 24)


figure()
plot(lon(:,1)-360,c45_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(68,1)-360 -0.6 0.45])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y NH4 a 45°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c45_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(68,1)-360 -0.6 0.45])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y NH4 a 45°S','FontSize', 24)

%ANTICICLON

antic=load('ASPS.txt');
%solo el año que nos importa

antic=antic(1:34,:); %desde la columna 2 a la 13 son los meses y la 14 es el promedio anual

anticiclon=antic(:,2:13)';
anticiclon=anticiclon(:);

NH4_a_5a=NH4_a_5(121:end,:);
NH4_a_10a=NH4_a_10(121:end,:);
NH4_a_15a=NH4_a_15(121:end,:);
NH4_a_25a=NH4_a_25(121:end,:);
NH4_a_35a=NH4_a_35(121:end,:);
NH4_a_45a=NH4_a_45(121:end,:);


for i=1:ddd(2)
    NH4_a_5_2a(:,i)=mediamovil(NH4_a_5a(:,i),2);
    NH4_a_10_2a(:,i)=mediamovil(NH4_a_10a(:,i),2);
    NH4_a_15_2a(:,i)=mediamovil(NH4_a_15a(:,i),2);
    NH4_a_25_2a(:,i)=mediamovil(NH4_a_25a(:,i),2);
    NH4_a_35_2a(:,i)=mediamovil(NH4_a_35a(:,i),2);
    NH4_a_45_2a(:,i)=mediamovil(NH4_a_45a(:,i),2);
    c5_anti(:,i)=corr(NH4_a_5_2a(:,i),anticiclon);
    c10_anti(:,i)=corr(NH4_a_10_2a(:,i),anticiclon);
    c15_anti(:,i)=corr(NH4_a_15_2a(:,i),anticiclon);
    c25_anti(:,i)=corr(NH4_a_25_2a(:,i),anticiclon);
    c35_anti(:,i)=corr(NH4_a_35_2a(:,i),anticiclon);
    c45_anti(:,i)=corr(NH4_a_45_2a(:,i),anticiclon);
end

figure()
plot(lon(:,1)-360,c5_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(61,1)-360 -1 0.3])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y NH4 a 5°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c10_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(64,1)-360 -1 0.3])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y NH4 a 10°S','FontSize', 24)


figure()
plot(lon(:,1)-360,c15_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(66,1)-360 -1 0.3])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y NH4 a 15°S','FontSize', 24)


figure()
plot(lon(:,1)-360,c25_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(71,1)-360 -0.6 0.45])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y NH4 a 25°S','FontSize', 24)


figure()
plot(lon(:,1)-360,c35_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(69,1)-360 -0.6 0.45])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y NH4 a 35°S','FontSize', 24)


figure()
plot(lon(:,1)-360,c45_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(68,1)-360 -0.6 0.45])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y NH4 a 45°S','FontSize', 24)





%giff
timee=datevec(time);
timee=timee(:,1:2);

figure(1)
filename= 'ggifNH4f.gif' ;
set(gcf,'Renderer','zbuffer')
for i=1:D(3)
    pcolor(lon,lat,NH4_af(:,:,i))
    colorbar
    caxis([-2 2])
    shading interp
    colormap(cmap)
    c = colorbar;
    c.Label.String = '[mmol/m^3]';
    hold on
    contour(lon,lat,NH4_af(:,:,i),[0 0],'k','Linewidth',1)
    title(['anomalías mensuales de NH4 filtrado a 1100m en ZMO ' mat2str(timee(i,:))])
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


%% NO3 Filtrado
clear all
clc 

load 'ONI.mat'
load 'NO3.mat'
load 'sam.mat'
time1=datevec(time);
D=size(NO3);


%le quitamos el ciclo anual y estandarizamos 
NO3_a=cicloanual_estan(NO3,time);


%filtro de 11 meses usando el metodo coseno-lanczos
%primera pasada
VAR = NO3_a; 
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

NO3_af=VAR_ft1;


%tomamos una lat 
for i=1:D(3)
    NO3_a_5=squeeze(NO3_af(:,115,:)); %5
    NO3_a_10=squeeze(NO3_af(:,96,:)); %10
    NO3_a_15=squeeze(NO3_af(:,77,:)); %15
    NO3_a_25=squeeze(NO3_af(:,49,:)); %25
    NO3_a_35=squeeze(NO3_af(:,30,:)); %35
    NO3_a_45=squeeze(NO3_af(:,11,:)); %45
end

cmap = colormap_cpt('Optimus_Prime.cpt',265);

figure()
pcolor(lon(:,1)-360,time,NO3_a_5')
datetick('y')
axis([lon(1,1)-360 lon(61,1)-360 715885 731945])
caxis([-3 3])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,NO3_a_5',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías de NO3 en 5°S','FontSize', 22)


figure()
pcolor(lon(:,1)-360,time,NO3_a_10')
datetick('y')
axis([lon(1,1)-360 lon(64,1)-360 715885 731945])
caxis([-3 3])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,NO3_a_10',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías de NO3 en 10°S','FontSize', 22)


figure()
pcolor(lon(:,1)-360,time,NO3_a_15')
datetick('y')
axis([lon(1,1)-360 lon(66,1)-360 715885 731945])
caxis([-3 3])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,NO3_a_15',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías de NO3 en 15°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,NO3_a_25')
datetick('y')
axis([lon(1,1)-360 lon(71,1)-360 715885 731945])
caxis([-3 3])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,NO3_a_25',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías de NO3 en 25°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,NO3_a_35')
datetick('y')
axis([lon(1,1)-360 lon(69,1)-360 715885 731945])
caxis([-2.5 2.5])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,NO3_a_35',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías de NO3 en 35°S','FontSize', 22)

figure()
pcolor(lon(:,1)-360,time,NO3_a_45')
datetick('y')
axis([lon(1,1)-360 lon(68,1)-360 715885 731945])
caxis([-2.5 2.5])
shading interp
colormap(cmap)
c = colorbar;
c.Label.String = '[mmol/m^3]';
hold on
contour(lon(:,1)-360,time,NO3_a_45',[0 0],'k','Linewidth',1)
xlabel('Longitud °W')
ylabel('Tiempo')
title('Anomalías de NO3 en 45°S','FontSize', 22)

% correlacion entre ONI y sam  y anomalias de nitrato

NO3_a_5=NO3_a_5';
NO3_a_10=NO3_a_10';
NO3_a_15=NO3_a_15';
NO3_a_25=NO3_a_25';
NO3_a_35=NO3_a_35';
NO3_a_45=NO3_a_45';
ddd=size(NO3_a_5);

for i=1:ddd(2)
    NO3_a_5_2(:,i)=mediamovil(NO3_a_5(:,i),3);
    NO3_a_10_2(:,i)=mediamovil(NO3_a_10(:,i),3);
    NO3_a_15_2(:,i)=mediamovil(NO3_a_15(:,i),3);
    NO3_a_25_2(:,i)=mediamovil(NO3_a_25(:,i),3);
    NO3_a_35_2(:,i)=mediamovil(NO3_a_35(:,i),3);
    NO3_a_45_2(:,i)=mediamovil(NO3_a_45(:,i),3);
    c5_oni(:,i)=corr(NO3_a_5_2(:,i),ONI);
    c10_oni(:,i)=corr(NO3_a_10_2(:,i),ONI);
    c15_oni(:,i)=corr(NO3_a_15_2(:,i),ONI);
    c25_oni(:,i)=corr(NO3_a_25_2(:,i),ONI);
    c35_oni(:,i)=corr(NO3_a_35_2(:,i),ONI);
    c45_oni(:,i)=corr(NO3_a_45_2(:,i),ONI);
    c5_sam(:,i)=corr(NO3_a_5_2(:,i),sam);
    c10_sam(:,i)=corr(NO3_a_10_2(:,i),sam);
    c15_sam(:,i)=corr(NO3_a_15_2(:,i),sam);
    c25_sam(:,i)=corr(NO3_a_25_2(:,i),sam);
    c35_sam(:,i)=corr(NO3_a_35_2(:,i),sam);
    c45_sam(:,i)=corr(NO3_a_45_2(:,i),sam);
end

figure()
plot(lon(:,1)-360,c5_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(61,1)-360 -0.7 0.7])
xlabel('longitud °W')
ylabel('correlacion')
title('Correlacion entre ONI y NO3 a 5°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c5_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(61,1)-360 -0.7 0.7])
xlabel('longitud °W')
ylabel('correlacion')
title('Correlacion entre SAM y NO3 a 5°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c10_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(64,1)-360 -0.7 0.7])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y NO3 a 10°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c10_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(64,1)-360 -0.7 0.7])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y NO3 a 10°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c15_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(66,1)-360 -0.7 0.7])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y NO3 a 15°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c15_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(66,1)-360 -0.7 0.7])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y NO3 a 15°S','FontSize', 24)


figure()
plot(lon(:,1)-360,c25_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(71,1)-360 -0.6 0.6])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y NO3 a 25°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c25_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(71,1)-360 -0.6 0.6])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y NO3 a 25°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c35_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(69,1)-360 -0.7 0.7])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y NO3 a 35°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c35_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(69,1)-360 -0.7 0.7])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y NO3 a 35°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c45_oni,'LineWidth',4)
axis([lon(1,1)-360 lon(68,1)-360 -0.7 0.7])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre ONI y NO3 a 45°S','FontSize', 24)
figure()
plot(lon(:,1)-360,c45_sam,'LineWidth',4)
axis([lon(1,1)-360 lon(68,1)-360 -0.7 0.7])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre SAM y NO3 a 45°S','FontSize', 24)



%ANTICICLON
antic=load('ASPS.txt');
%solo el año que nos importa

antic=antic(1:34,:); %desde la columna 2 a la 13 son los meses y la 14 es el promedio anual

anticiclon=antic(:,2:13)';
anticiclon=anticiclon(:);


NO3_a_5a=NO3_a_5(121:end,:);
NO3_a_10a=NO3_a_10(121:end,:);
NO3_a_15a=NO3_a_15(121:end,:);
NO3_a_25a=NO3_a_25(121:end,:);
NO3_a_35a=NO3_a_35(121:end,:);
NO3_a_45a=NO3_a_45(121:end,:);

for i=1:ddd(2)
    NO3_a_5_2a(:,i)=mediamovil(NO3_a_5a(:,i),2);
    NO3_a_10_2a(:,i)=mediamovil(NO3_a_10a(:,i),2);
    NO3_a_15_2a(:,i)=mediamovil(NO3_a_15a(:,i),2);
    NO3_a_25_2a(:,i)=mediamovil(NO3_a_25a(:,i),2);
    NO3_a_35_2a(:,i)=mediamovil(NO3_a_35a(:,i),2);
    NO3_a_45_2a(:,i)=mediamovil(NO3_a_45a(:,i),2);
    c5_anti(:,i)=corr(NO3_a_5_2a(:,i),anticiclon);
    c10_anti(:,i)=corr(NO3_a_10_2a(:,i),anticiclon);
    c15_anti(:,i)=corr(NO3_a_15_2a(:,i),anticiclon);
    c25_anti(:,i)=corr(NO3_a_25_2a(:,i),anticiclon);
    c35_anti(:,i)=corr(NO3_a_35_2a(:,i),anticiclon);
    c45_anti(:,i)=corr(NO3_a_45_2a(:,i),anticiclon);
end

figure()
plot(lon(:,1)-360,c5_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(61,1)-360 -0.7 0.7])
xlabel('longitud °W')
ylabel('correlacion')
title('Correlacion entre Índice del Anticiclón y NO3 a 5°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c10_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(64,1)-360 -0.7 0.7])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y NO3 a 10°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c15_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(66,1)-360 -0.7 0.7])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y NO3 a 15°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c25_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(71,1)-360 -0.6 0.6])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y NO3 a 25°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c35_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(69,1)-360 -0.7 0.7])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y NO3 a 35°S','FontSize', 24)

figure()
plot(lon(:,1)-360,c45_anti,'LineWidth',4)
axis([lon(1,1)-360 lon(68,1)-360 -0.7 0.7])
xlabel('longitud °W')
ylabel('correlación')
title('Correlación entre Índice del Anticiclón y NO3 a 45°S','FontSize', 24)




%giff
timee=datevec(time);
timee=timee(:,1:2);

figure(1)
filename= 'ggifNO3f.gif' ;
set(gcf,'Renderer','zbuffer')
for i=1:D(3)
    pcolor(lon,lat,NO3_af(:,:,i))
    caxis([-2.5 2.5])
    shading interp
    colormap(cmap)
    c = colorbar;
    c.Label.String = '[mmol/m^3]';
    hold on
    contour(lon,lat,NO3_af(:,:,i),[0 0],'k','Linewidth',1)
    title(['anomalías mensuales de NO3 filtrado a 1100m en ZMO ' mat2str(timee(i,:))])
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

