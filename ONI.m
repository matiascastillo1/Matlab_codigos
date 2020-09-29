%% CALCULAR EL INDICE DEL NIÑO OCEANICO (ONI)

%se parte cargando los datos de temperatura superficial del mar de la
%region conocida como niño 3.4 (5N-5S, 170W-120W)

load 'SST.mat'

%se usa datevec para poder ver el tiempo de mejor manera
time1=datevec(time);
D=size(SST);
%se extrae el ciclo anual para trabajar con anomalias

enero=find(time1(:,2)==1);
febrero=find(time1(:,2)==2);
marzo=find(time1(:,2)==3);
abril=find(time1(:,2)==4);
mayo=find(time1(:,2)==5);
junio=find(time1(:,2)==6);
julio=find(time1(:,2)==7);
agos=find(time1(:,2)==8);
sept=find(time1(:,2)==9);
oct=find(time1(:,2)==10);
nov=find(time1(:,2)==11);
dic=find(time1(:,2)==12);

SST_a=ones(D(1),D(2),D(3));
for i=1:D(1)
    for j=1:D(2)
        SST_a(i,j,enero)=(SST(i,j,enero)-nanmean(SST(i,j,enero)));%./nanstd(SST(i,j,:));
        SST_a(i,j,febrero)=(SST(i,j,febrero)-nanmean(SST(i,j,febrero)));%./nanstd(SST(i,j,:));
        SST_a(i,j,marzo)=(SST(i,j,marzo)-nanmean(SST(i,j,marzo)));%./nanstd(SST(i,j,:));
        SST_a(i,j,abril)=(SST(i,j,abril)-nanmean(SST(i,j,abril)));%./nanstd(SST(i,j,:));
        SST_a(i,j,mayo)=(SST(i,j,mayo)-nanmean(SST(i,j,mayo)));%./nanstd(SST(i,j,:));
        SST_a(i,j,junio)=(SST(i,j,junio)-nanmean(SST(i,j,junio)));%./nanstd(SST(i,j,:));
        SST_a(i,j,julio)=(SST(i,j,julio)-nanmean(SST(i,j,julio)));%./nanstd(SST(i,j,:));
        SST_a(i,j,agos)=(SST(i,j,agos)-nanmean(SST(i,j,agos)));%./nanstd(SST(i,j,:));
        SST_a(i,j,sept)=(SST(i,j,sept)-nanmean(SST(i,j,sept)));%./nanstd(SST(i,j,:));
        SST_a(i,j,oct)=(SST(i,j,oct)-nanmean(SST(i,j,oct)));%./nanstd(SST(i,j,:));
        SST_a(i,j,nov)=(SST(i,j,nov)-nanmean(SST(i,j,nov)));%./nanstd(SST(i,j,:));
        SST_a(i,j,dic)=(SST(i,j,dic)-nanmean(SST(i,j,dic)));%./nanstd(SST(i,j,:));
    end
end

%ahora se promedia espacialmente para asi obtener una serie de tiempo de la
%temperatura superficial en el niño 3.4

for i=1:D(3)
    SST_1(:,i)=reshape(squeeze(SST_a(:,:,i)),D(1)*D(2),1);
end
for i=1:D(3)
    SST_ap(i,1)=mean(SST_1(:,i));
end

%una vez tenemos nuestro vector temporal de anomalias de temperatura de la
%zona niño 3.4 se le aplica una media movil de ventana 3


ONI=mediamovil(SST_ap,3);  %INDICE DEL NIÑO OCEANICO

% Asi, se calcula el indice del niño oceanico (ONI) de acuerdo a la
% definicion empleada de la NOAA
%https://climatedataguide.ucar.edu/climate-data/nino-sst-indices-nino-12-3-34-4-oni-and-tni?fbclid=IwAR2w8a7LBMZs475qRBvx2Zuw6yGnISQao7OEZzLEDQcyrsZfOw09uzn1D00


%%
%Se puede hacer unos pequeños plots
%En el caso del ONI se considera un fenomeno del Niño cuando las anomalias
%son mayores a 0.5C y Niña cuando son menores a -0.5C

figure()
h=line([715916 731915],[0.5 0.5])
h0=line([715916 731915],[-0.5 -0.5])
h1=line([715916 731915],[1.5 1.5])
h2=line([715916 731915],[-1.5 -1.5])
hold on
plot(time(2:end-1),ONI,'k')
set(h,'color','r')
set(h1,'color','r')
set(h0,'color','b')
set(h2,'color','b')
datetick
axis 'tight'
title('Indice del Niño Oceanico (ONI) desde 1960 hasta 2005')
ylabel('°C')
legend('El Niño','La Niña')
