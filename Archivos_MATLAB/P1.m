clear all
%% Carga de datos iniciales
load('OngoingEEG_data.mat')
length = 334848;
f_samp = 1024;
data = transpose(EEGData);

ventana1 = 0.5*f_samp; %0.5[s]
ventana2 = 1.5*f_samp; %1.5[s]
ventana3 = 3*f_samp; %3[s]
%% Gráfico señal original en frecuencia
eje_x = linspace(-length/2,length/2-1,length)*f_samp/length;
magnitud = fft(data)/length;
magnitud_abs = abs([magnitud(length/2:length) magnitud(1:length/2-1)]);

figure(1)
plot(eje_x,magnitud_abs) 
xlim([-40 40])
xlabel('f[Hz]')
grid on

%% ventana 1
prom1 = promediacion(data,ventana1,length);
[values1,f1] = FFT(prom1,ventana1); %f1 es un vector del largo de la cantidad de datos de la ventana
figure(2)
plot(f1,values1)
xlim([-40 40])
xlabel('f[Hz]')
ylabel('Amplitud')
grid on

%% ventana 2
prom2 = promediacion(data,ventana2,length);
[values2,f2] = FFT(prom2,ventana2); %f2 es un vector del largo de la cantidad de datos de la ventana

figure(3)
plot(f2,values2)
xlim([-40 40])
xlabel('f[Hz]')
ylabel('Amplitud')
grid on

%% ventana3
prom3 = promediacion(data,ventana3,length);
[values3,f3] = FFT(prom3,ventana3); %f3 es un vector del largo de la cantidad de datos de la ventana

figure(4)
plot(f3,values3)
xlim([-40 40])
xlabel('f[Hz]')
ylabel('Amplitud')
grid on
%% wavelet señal original en reposo
figure (5)
cwt(data,f_samp)

%% ventana 0.5[s] usando Wavelet
figure(6)
cwt(prom1,f_samp)

%% ventana 1.5[s] usando wavelet
figure(7)
cwt(prom2,f_samp)

%% ventana 3[s] usando wavelet
figure(8)
cwt(prom3,f_samp)

%% funciones
function potencia_media = calculo_psd(Datos_abs,f_sampling,largo)
    potencia_media = (Datos_abs.^2)/(largo*f_sampling);

end

function val_promedio = promediacion(datos,VENTANA,largo)
    total = 0;
    cont = 0;
    j = 1;
    
    while j+VENTANA-1<=largo
        cont = cont+1;
        ventana_i=datos(j:j+VENTANA-1);
        total = ventana_i+total;
        
        j=j+VENTANA;
    end
    val_promedio = total/cont;
end

function [fourier,f] = FFT(DATA,LARGO)
    f = linspace(-LARGO/2,LARGO/2-1,LARGO)*1024/LARGO;
    fourier_aux = fft(DATA)/LARGO;
    fourier = abs([fourier_aux(LARGO/2:LARGO) fourier_aux(1:LARGO/2-1)]);
end