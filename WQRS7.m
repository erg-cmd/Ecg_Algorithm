%% WQRS7 ---
% Implementacion del algoritmo WQRS descripto en "Robust-OpenSource
% Algorithm to detect the onset and duration complexes".
% Se modifico para funcionar con el ECG-KIT

%  usage: [QRS] = WQRS6(ECG,HEADER,PROGRESO,PAYLOAD)
%  param in
%			ECG: Matriz de nsig * nsamp
%			HEADER: descripcion de la señal
%			PROGRESO: ----barrita
%			PAYLOAD: -----?
%  param out
%			QRS:Cell que contiene las posiciones en muestras de las detecciones 
%			MULTI: segun las convenciones multileed que desconocemos....
%/*------------------------------------------------------------------*/
function [QRS,MULTI]=WQRS7(ECG, HEADER,PROGRESO,PAYLOAD)
	
	% Comentar! Solo debug!
	t=[1:HEADER.nsamp]/HEADER.freq; 
	
    % Ajustamos las muestras a 0
	for k = 1:HEADER.nsig
		%~ %ECG_filtered(:,k) = (ECG_filtered(:,k) - HEADER.adczero(k)) ./ HEADER.gain(k);
		ECG(:,k) = (ECG(:,k) - HEADER.adczero(k));
    end
    
	% Un objeto filtro con 
	H = Filtro_FIR(HEADER.freq);
	%Aplicamos el filtro
	ECG_filtered=filter(H,ECG);	
	[phi,w]=phasedelay(H);
	PhaseDelay = round(phi(5));
	%Corregimos las muestras x el atraso del filtro
	ECG_filtered=circshift(ECG_filtered,[-PhaseDelay 0]);
	ECG_filtered=2*ECG_filtered(1:end-PhaseDelay,:);
	
	% retorna la Longitud de la Curva de cada señal
	LT = CurveLengthTrans(ECG_filtered,HEADER.nsamp-PhaseDelay,HEADER.nsig, HEADER.freq);
	minimo = min(LT);
	
	% Leeds individuales
	QRS = cell(HEADER.nsig,1);
	for n=1:HEADER.nsig
		LT(:,n) = LT(:,n) - minimo(n);
		%~ QRS{n} = BloqueDecision(LT(:,n),HEADER.freq,HEADER.gain(n));
		deteccion = BloqueDecision(LT(:,n),HEADER.freq,HEADER.gain(n));
		QRS{n} = deteccion;		% retorna las detecciones, la posicion es en muestras

		%Solo debug, desactivar esta parte!!!
		figure(n);
		plot(t(1:end-PhaseDelay),LT(:,n));title(['Length Transform ',num2str(n),'Leed']);grid on;hold on;xlabel('Time(sec)');
		plot(t(deteccion),LT(deteccion,n),'c*','Markersize',5);
		for index=1:length(deteccion)
			text(t(deteccion(index)),LT(deteccion(index),n),[num2str(index)]);
		end
	end
	
	% Sumando los leeds
	LTsum = sum(LT,2);
	MULTI = BloqueDecision(LTsum,HEADER.freq,HEADER.gain(1));
	
	% Solo debug, desactivar esta parte!
	%~ figure(3);
	%~ plot(t,LTsum);title('Length Transform (sum)');grid on;hold on;xlabel('Time(sec)');
	%~ plot(t(MULTI),LTsum(MULTI),'c*','Markersize',5);
	%~ for index=1:length(MULTI)
		%~ text(t(MULTI(index)),LTsum(MULTI(index)),[num2str(index)]);
	%~ end
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
	
	
	Fpass = 50;              % Passband Frequency
	Fstop = 160;              % Stopband Frequency
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
%% CurveLengthTrans - Multileed
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

function [largo]=CurveLengthTrans(muestras, nsamples,nsig, fsamp)

	ventana = 0.007;				% tiempo de la ventana en segundos >>0.02
	c = 1.25*200*200*0.1;				% Factor de escalamiento en segundos >>1.25
	w = round(ventana*fsamp); 		% ancho de la ventana en muestras
	t = c/fsamp;					% periodo de muestreo o "factor de escalamiento"
	largo = zeros(nsamples,nsig);	% Vector "largo" del tamaño de muestras llenado con 0
	suma = zeros(1,nsig);

	for m=1:(nsamples-w)
		suma = sum(muestras(m:m+w,:),1);
		largo(m,:) = ((t.^2)+(suma.^2)).^(1/2); 
		suma = zeros(1,nsig);
	end
	
	% las muestras que me no llegan a procesarce por el tamaño del 
	% filtro, busco el min de la Lt
	for m=0:w
		largo(end-m,:)=min(largo(1:end-w,:)); 
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
	Ts = 1/Fs;
	EyeClosing = round(.25/Ts);	% periodo para evitar el sobredetección(ms)
	SegmQRS = round(.01/Ts);	% Segun bibliografia consultada 130ms
	nsamp=length(signal);		% numero de muestras que tiene la señal    
	indice=1; 					% barre e indica posicion en el estudio
	N = 1; 						% num de complejos detectados
	posMax = 1;posMin=1; 		% posicion del maximo/minimo
	Tm = .100 * gain;			% Definimos el Threshold minimo de 100uV, expresado mV.
	DistMinima = 50;			% Consideracion >>>>>>>>>>
	ExpectedRR = Fs * 2.5;		% Intervalo esperado del QRS (2.5seg MAX)
	tiempo = 1;					% Contabiliza el tiempo desde el ultimo umbral valido
	
	% Ventana de busqueda
	ventana = ceil(ExpectedRR/2);
	% Definimos un Umbral promedio
	T0 = mean(signal(1:8*Fs));	% Promediar los primeros 8 segundos
	% Fijamos nivel de umbral actual
	Ta = 2 * T0;
	T1 = T0;

	%/*-1er Complejo--------------------------------------------------*/
	% Dentro de 1 IntervaloRR debe estar 1 complejo
	posMax = find(signal(indice:indice+ExpectedRR) == max(signal(indice:indice+ExpectedRR)));			
	
	% Si el QRS esta pegado al inicio del record, error en la busqueda 
	% si excedes el inicio del record
	if (posMax<SegmQRS)
		posMin = indice + find(signal(indice:posMax) == min(signal(indice:posMax)));
	else
		posMin = (posMax-SegmQRS) + find(signal(posMax-SegmQRS:posMax) == min(signal(posMax-SegmQRS:posMax)));
	end
	posMin = posMin(end);
	% Posicion QRS
	beat(N) = round ((posMax+posMin)/2);
	% Desplazamiento para evitar sobredeteccion
	indice = beat(N) + EyeClosing;
	% Siguiente deteccion
	N = 2;
	
	%/*-Barremos todo el registro hasta end-ventana-------------------*/
	while ((indice+ventana) < nsamp) 
		% Posible QRS 
		if(max(signal(indice:indice+ventana)) > T1)
			% Tiempo sin detectar cambio 
			tiempo = 0;
			% devuelve un vector con la pos del Max
			posMax = find(signal(indice:indice+ventana) > max(signal(indice:indice+ventana))*0.7);
			% me quedo con el primer max
			posMax = indice + posMax(1);
			posMin = (posMax-SegmQRS) + find(signal(posMax-SegmQRS:posMax) == min(signal(posMax-SegmQRS:posMax)));
			% Si es un vector no quedamos con el mas próximo a Q
			posMin = posMin(end);								
			% entonces lo consideramos un QRS
			if(signal(posMax) - signal(posMin) > DistMinima)	
				% posicion QRS
				beat(N) = round((posMax+posMin)/2); 			
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
			tiempo = tiempo + ceil(ExpectedRR/400);
			% y si tiempo es > Expec
			if (Ta > Tm) && (tiempo > ExpectedRR/2) 
				% Bajamos el umbral
				Ta = Ta - Ta*.1;
				T1 = Ta/3; 
			end
		end
		% nos corremos ExpectedRR/3
		indice = indice + ceil(ExpectedRR/400);
	end
	
	%/*-Tramo que falta desde indice hasta end------------------------*/
	% Ultimo QRS
	if(max(signal(end-ventana:end)) > Ta)													
		% reubicamos indice
		indice = length(signal)-ventana;														
		% Retorna las pos relativas del intervalo
		posMax = find(signal(indice:indice+ventana) > max(signal(indice:indice+ventana))*0.7); 	
		posMax = indice + posMax(end);
		posMin = (posMax-SegmQRS) + find(signal(posMax-SegmQRS:posMax) == min(signal(posMax-SegmQRS:posMax)));
		% Si es un vector no quedamos con el mas próximo a Q
		posMin = posMin(end);
		% entonces lo consideramos un QRS
		if(signal(posMax) - signal(posMin) > DistMinima)
			% posicion QRS
			beat(N) = round((posMax+posMin)/2);
			% Actualizamos al valor del último pico
			Ta = signal(posMax)/2;
			% Desplazamiento para evitar sobredetección
			indice = beat(N) + EyeClosing;
			% Siguiente deteccion
			N = N + 1;
		end
	end
end

%[EOF]
