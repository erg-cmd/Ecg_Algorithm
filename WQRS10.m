%% WQRS10 ---
% Implementacion del algoritmo WQRS descripto en "Robust-OpenSource
% Algorithm to detect the onset and duration complexes".
% Se modifico para funcionar con el ECG-KIT

%  usage: [QRS] = WQRS10(ECG,HEADER,PROGRESO,PAYLOAD)
%  param in
%			ECG: Matriz de nsig * nsamp
%			HEADER: descripcion de la señal
%			PROGRESO: ----barrita
%			PAYLOAD: -----?
%  param out
%			QRS:Cell que contiene las posiciones en muestras de las detecciones 
%			MULTI: segun las convenciones multileed que desconocemos....
%/*------------------------------------------------------------------*/
function [QRS,MULTI]=WQRS10(ECG,HEADER,PROGRESO,PAYLOAD)
	ModoDebug = false;
    %ModoDebug = true;
	if(ModoDebug)
		t=[1:HEADER.nsamp]/HEADER.freq; 
    end
    
    PwrLnFrq = 60;						% Frecuencia de la linea AC
	LP2n = 2*HEADER.freq/PwrLnFrq;		% Notch a la frecuencia de linea
	% Un objeto filtro con 
	H = Filtro_FIR(HEADER.freq);
	%Aplicamos el filtro
	ECG_filtered=filter(H,ECG);
    ECG_filtered = 17.7828 .* ECG_filtered;
	[phi,w]=phasedelay(H);
	PhaseDelay = 2*round(phi(5));
	%Corregimos las muestras x el atraso del filtro
	ECG_filtered=circshift(ECG_filtered,[-PhaseDelay 0]);
	
	%ajuste de muestras por delay del filtro
	for k=1:HEADER.nsig
		ECG_filtered(HEADER.nsamp-PhaseDelay:HEADER.nsamp,k) = ECG_filtered(HEADER.nsamp-PhaseDelay-1,k);
	end
	%Filtrado a la frecuencia de linea; Una muestra de relleno
	dy = diff(ECG_filtered,1,1)/LP2n;
	dy = vertcat(dy,dy(end,:));
	% retorna la Longitud de la Curva de cada señal
	LT = CurveLengthTrans(dy,HEADER,ModoDebug);
	minimo = min(LT);	
	% Leads individuales
	QRS = cell(HEADER.nsig,1);
	for n=1:HEADER.nsig
        % Se corrige el corrimiento del zero
		if(minimo(n) < 0) 
            LT(:,n) = LT(:,n) + abs(minimo(n));   
        end
        if(~ModoDebug)
			QRS{n} = BloqueDecision(LT(:,n),HEADER.freq,HEADER.gain(n));
		else
			deteccion = BloqueDecision(LT(:,n),HEADER.freq,HEADER.gain(n));
			QRS{n} = deteccion;		% retorna las detecciones, la posicion es en muestras
			figure(n);
			plot(t,LT(:,n));title(['Length Transform ',num2str(n),'lead']);grid on;hold on;xlabel('Time(sec)');
			plot(t(deteccion),LT(deteccion,n),'c*','Markersize',5);
			for index=1:length(deteccion)
				text(t(deteccion(index)),LT(deteccion(index),n),[num2str(index)]);
			end
		end
	end
	
	% Sumando los leads
	LTsum = sum(LT,2);
	MULTI = BloqueDecision(LTsum,HEADER.freq,HEADER.gain(1));
end


%/*-------------------------------------------------------------------*/
%%FIR_Equiri Returns a discrete-time filter object.
% Equiripple Lowpass filter designed using the FIRPM function.

%/*-------------------------------------------------------------------*/
function Hd = Filtro_FIR(Fs)
	% All frequency values are in Hz.
	%~ Fpass = 14;              % Passband Frequency
	%~ Fstop = 70;              % Stopband Frequency
	%~ Fpass = 10;              % Passband Frequency
	%~ Fstop = 50;              % Stopband Frequency
	%~ Dpass = 0.057501127785;  % Passband Ripple% 
	%~ Dstop = 0.00001;          % Stopband Attenuation
	%~ dens  = 20;              % Density Factor

	Fpass = 16;              % Passband Frequency
	Fstop = 70;              % Stopband Frequency
	Dpass = 0.057501127785;  % Passband Ripple
	Dstop = 0.0001;          % Stopband Attenuation
	dens  = 20;              % Density Factor
	
	% Calculate the order from the parameters using FIRPMORD.
	[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
	
	% Calculate the coefficients using the FIRPM function.
	b  = firpm(N, Fo, Ao, W, {dens});
	Hd = dfilt.dffir(b);
end


%/*-------------------------------------------------------------------*/
%% CurveLengthTrans - Multilead
%
%	Se pretende devolver la matriz longitud de una curva y
%	especificando una ventana (w), número de muestras(N), periodo de
%	muestreo (1/fs = t). Condición 1+w < i < N
%
% 	param in:	
%		muestras: señal/es a procesar
%		nsamples: numero de muestras
%		nsig:	numero de señales
% 	param out:	
%		void
%/*-------------------------------------------------------------------*/

function [ltprom]=CurveLengthTrans(muestras,header,opcion)
	
	MaxQRSw = 0.13;							% tiempo de la ventana en segundos >>0.13
	c = 1.25*header.gain(1)*header.gain(1);	% Factor de escalamiento en segundos >>1.25
	w = round(MaxQRSw*header.freq); 		% ancho de la ventana en muestras
	t = c/header.freq;						% periodo de muestreo o "factor de escalamiento"
	
	lt = zeros(header.nsamp,header.nsig);	% Vector "lt" del tamaño de muestras llenado con 0
	ltprom = zeros(header.nsamp,header.nsig);	% Vector "LTprom" del tamaño de muestras llenado con 0
	delta = zeros(1,header.nsig);	% Acumula la diferencia entre la muestra actual y actual-w
	
	for m=1:header.nsamp
		lt(m,:) = (t+(muestras(m).^2)).^(1/2); 
		if m > w
			delta = delta + lt(m,:) - lt(m-w,:);
			ltprom(m,:) = delta;
        else
                delta = delta + lt(m,:) - lt(1,:);
                ltprom(m,:) = delta;
		end
	end
	
	% las muestras que no llegan a procesarce por el tamaño del 
	% filtro, busco el min de la Lt
	%~ for m=0:w
		%~ lt(end-m,:)=min(lt(1:end-w,:)); 
    %~ end
    if opcion
        figure(4);
        %~ plot([muestras(:,1), lt(:,1)*max(abs([max(muestras(:,1)),min(muestras(:,1))]))/max(abs([max(lt(:,1)),min(lt(:,1))]))]);
        plot([muestras(:,1), ltprom(:,1)*max(abs([max(muestras(:,1)),min(muestras(:,1))]))/max(abs([max(ltprom(:,1)),min(ltprom(:,1))]))]);
        grid on; hold on;
        figure(5);
        %~ plot([muestras(:,2), lt(:,2)*max(abs([max(muestras(:,2)),min(muestras(:,2))]))/max(abs([max(lt(:,2)),min(lt(:,2))]))]);
        plot([muestras(:,2), ltprom(:,2)*max(abs([max(muestras(:,2)),min(muestras(:,2))]))/max(abs([max(ltprom(:,2)),min(ltprom(:,2))]))]);
        grid on; hold on;
    end
end

%/*------------------------------------------------------------------*/
%% Bloque de Decision
%
%	usage: [beat] = BloqueDecision(signal,Fs,gain)	
% version6 - 

% Solamente se guardan las posiciones de los posibles QRS
%/*------------------------------------------------------------------*/

function [beat] = BloqueDecision(signal,Fs,gain)
	EyeClosing = round(.25*Fs);	% periodo para evitar el sobredetección(ms)
	nsamp=length(signal);		% numero de muestras que tiene la señal    
	indice=1; 					% barre e indica posicion en el estudio
	N = 1; 						% num de complejos detectados
	posMax = 1;posMin=1; 		% posicion del maximo/minimo
	Tm = .100 * gain;			% Definimos el Threshold minimo de 100uV, expresado mV.
	DistMinima = 10;			% Distancia minima entre un Max-Min para ser QRS
	ExpectedRR = Fs * 2.5;		% Intervalo esperado del QRS (2.5seg MAX)
	tiempo = 1;					% Contabiliza el tiempo desde el ultimo umbral valido
	
	% Definimos un Umbral promedio
	T0 = mean(signal(1:8*Fs));	% Promediar los primeros 8 segundos
	% Fijamos nivel de umbral actual
	Ta = 3*T0;
	T1 = 2*T0;

	%/*-Barremos todo el registro hasta end-ventana-------------------*/
	while ((indice+EyeClosing) < nsamp) 
		% Posible QRS 
		if(max(signal(indice:indice+EyeClosing/2)) > T1)
			% Tiempo sin detectar cambio 
			tiempo = 0;
			% devuelve un vector con la pos del Max
			posMax = find(signal(indice:indice+EyeClosing/2) == max(signal(indice:indice+EyeClosing/2)));
			% me quedo con el primer max
			posMax = indice + posMax(end);
			if(posMax > EyeClosing/2)
				posMin = (posMax-EyeClosing/2) + find(signal(posMax-EyeClosing/2:posMax-1) == min(signal(posMax-EyeClosing/2:posMax-1)));
			else
				posMin = 1 + find(signal(1:posMax) == min(signal(1:posMax)));
			end
			% Si es un vector no quedamos con el mas próximo a Q
			posMin = posMin(end);								
			% entonces lo consideramos un QRS
			if(signal(posMax) - signal(posMin) > DistMinima)	
				% posicion QRS
				beat(N) = posMax; 
				% Actualizamos al valor del último pico
				Ta = Ta + (signal(posMax)-Ta)/10;
				T1 = Ta/3;
				% Desplazamiento para evitar sobredetección
				indice = beat(N) + EyeClosing;					
				% Siguiente deteccion
				N = N + 1;										
			% No es un QRS
			end
		% No se detecto. Puede no haber o el Thres está alto!
		else
			tiempo = tiempo + ceil(EyeClosing/4);
			% y si tiempo es > Expec
			if (Ta > Tm) && (tiempo > ExpectedRR) 
				% Bajamos el umbral
				Ta = Ta - Ta*.1;
				T1 = Ta/3; 
			end
		end
		% nos corremos EyeClosing/4
		indice = indice + ceil(EyeClosing/4);
	end
	
	%/*-Tramo que falta desde indice hasta end------------------------*/
	% Ultimo QRS
	if(max(signal(end-EyeClosing/2:end)) > T1)							
		% reubicamos indice
		indice = length(signal)-EyeClosing/2;							
		% Retorna las pos relativas del intervalo
		posMax = find(signal(indice:indice+EyeClosing/2) == max(signal(indice:indice+EyeClosing/2))); 	
		posMax = indice + posMax(end);
		posMin = (posMax-EyeClosing/2) + find(signal(posMax-EyeClosing/2:posMax) == min(signal(posMax-EyeClosing/2:posMax)));
		% Si es un vector no quedamos con el mas próximo a Q
		posMin = posMin(end);
		% entonces lo consideramos un QRS
		if(signal(posMax) - signal(posMin) > DistMinima)
			% posicion QRS
			beat(N) = round((posMax+posMin)/2);
		end
	end
end

%[EOF]
