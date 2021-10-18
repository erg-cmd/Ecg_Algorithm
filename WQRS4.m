%/*------------------------------------------------------------------*/
%% WQRS4 - 
% Implementacion del algoritmo WQRS descripto en "Robust-OpenSource
% Algorithm to detect the onset and duration complexes".
% Se modifico para funcionar con el ECG-KIT

%  usage: [QRS] = WQRS4(ECG,HEADER,PROGRESO,PAYLOAD)
%  param in
%			ECG: Matriz de nsig * nsamp
%			HEADER: descripcion de la señal
%			PROGRESO: ----barrita
%			PAYLOAD: -----?
%  param out
%			QRS:Cell que contiene las posiciones en muestras de las detecciones 
%			MULTI: segun las convenciones multileed que desconocemos....
%/*------------------------------------------------------------------*/
function [QRS,MULTI]=WQRS4(ECG, HEADER,PROGRESO,PAYLOAD)
	
	%t=[1:HEADER.nsamp]/HEADER.freq; % Comentar! Solo debug!
	H = Filtro_FIR(HEADER.freq);
	ECG_filtered=filter(H,ECG);	%Aplicamos el filtro

	%[phi,w]=phasedelay(Filtro_FIR(HEADER.freq));			
	[phi,w]=phasedelay(H);									
	PhaseDelay = round(phi(5));								
	ECG_filtered=circshift(ECG_filtered,[-PhaseDelay 0]);	%para la ver 2013
	
	%ajuste de muestras por delay del filtro
	medio = mean(ECG_filtered(1:HEADER.nsamp-PhaseDelay,:));
	for k=1:HEADER.nsig
		ECG_filtered(HEADER.nsamp-PhaseDelay:HEADER.nsamp,k) = medio(k);
	end
	
	%for k=HEADER.nsamp-PhaseDelay:HEADER.nsamp
	%	ECG_filtered(k,:) = mean(ECG_filtered(1:HEADER.nsamp-PhaseDelay,:));
	%end
	
	for k = 1:HEADER.nsig									%ver el efecto antes y despues de LT. Igual es solo para ver la señal.
		ECG_filtered(:,k) = (ECG_filtered(:,k) - HEADER.adczero(k)) ./ HEADER.gain(k);	% x esc de mV
    end
    
	LT = CurveLengthTrans(ECG_filtered,HEADER.nsamp,HEADER.nsig, HEADER.freq);	% retorna la Longitud de la Curva de cada señal
	minimo = min(LT);
	%LTsum = sum(LT,2); 															% sumamos las señales por fila
	%LTsum = LTsum - minimo;														% bajamos el nivel de continua que viene de????
	%MULTI = BloqueDecision(LTsum,HEADER.freq);								% retorna las detecciones, la posicion es en muestras
	
	QRS = cell(HEADER.nsig,1);
	for n=1:HEADER.nsig
		LT(:,n) = LT(:,n) - minimo(n);
		QRS{n} = BloqueDecision(LT(:,n),HEADER.freq);					% retorna las detecciones, la posicion es en muestras
	end
	
	LTsum = sum(LT,2); 															% sumamos las señales por fila
	MULTI = BloqueDecision(LTsum,HEADER.freq);								% retorna las detecciones, la posicion es en muestras
	
	% en caso de funcionar, desactivar esta parte!
	%figure(2);
	%plot(t,LTsum);title('Length Transform (sum)');grid on;hold on;xlabel('Time(sec)');
	%plot(t(MULTI),LTsum(MULTI),'c*','Markersize',5);
	%for index=1:length(MULTI)		%barrido de los puntos
	%	text(t(MULTI(index)),LTsum(MULTI(index)),[num2str(index)]);
	%end
end


%/*-------------------------------------------------------------------*/
%%FIR_Equiri Returns a discrete-time filter object.
% Equiripple Lowpass filter designed using the FIRPM function.

%/*-------------------------------------------------------------------*/
function Hd = Filtro_FIR(Fs)
	% All frequency values are in Hz.
	Fpass = 14;              % Passband Frequency
	Fstop = 70;              % Stopband Frequency
	Dpass = 0.057501127785;  % Passband Ripple
	Dstop = 0.00001;          % Stopband Attenuation
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
	
	ventana = 0.02;						% tiempo de la ventana en segundos
	c = 0.9;							% Factor de escalamiento en segundos
	w = round(ventana*fsamp); 			% ancho de la ventana en muestras
	t = c*fsamp;						%periodo de muestreo o "factor de escalamiento"
	largo = zeros(nsamples,nsig);		%Vector "largo" del tamaño de muestras llenado con 0
	suma = zeros(1,nsig);

    for m=1:(nsamples-w)
        suma = sum(muestras(m:m+w,:),1);
        largo(m,:) = sqrt(power(t,2)+power(suma,2)); 
        suma = zeros(1,nsig);
    end
    
    % las muestras que me no llegan a procesarce por el tamaño del 
    %filtro, busco el min de la Lt
    for m=0:w
        largo(end-m,:)=min(largo(1:end-w,:)); 
    end
end

%/*------------------------------------------------------------------*/
%% Bloque de Decision
%
%	usage: [beat] = BloqueDecision(signal,Fs)	
% version5 - Dividida en 2 partes, determinación del nivel de umbral 
% para la posicion del jQRS (adaptativo) y estrategias de búsquedas
% dada la 1ra parte.

% Solamente se guardan las posiciones de los posibles QRS
%/*------------------------------------------------------------------*/

function [beat] = BloqueDecision(signal,Fs)
	Ts = 1/Fs;
	periodo = 5; 									% Intervalo a muestrear
	EyeClosing = round(.4/Ts); 						% periodo para evitar el sobredetección(ms)
	SegmQRS = round(.1/Ts);							% Segun bibliografia consultada 100ms
	nsamp=length(signal);							% numero de muestras que tiene la señal    
	indice=1; 										% barre e indica posicion en el estudio
	N = 1; 											% num de complejos detectados
	posMax = 1;posMin=1; 							% posicion del maximo/minimo
	Criterio_Perdido = 1.5;							% relacion entre deteccion actual y posterior en tiempo

%/*-------------------------------------------------------------------*/
	% Definimos el intervalo ExpectedRR para definir una ventana de busqueda
	UmbralProm = mean(signal(indice+(Fs*2):indice+(Fs*10))); 	%Definimos un Umbral mínimo
	ExpectedRR = IntervaloRR(indice, signal, Fs, UmbralProm);	%
	ventana = ceil(2 * ExpectedRR);								%
	Umbral = max(signal(indice:indice+ventana))/2; 				% Fijamos nivel de umbral 
	

    %% 1ra Parte - Deteccion por umbral adaptativo
    %/*-1er Complejo--------------------------------------------------*/
	posMax = find(signal(indice:indice+ExpectedRR) == max(signal(indice:indice+ExpectedRR)));			% Dentro de 1 IntervaloRR debe estar 1 complejo
	posMin = (indice*(posMax<SegmQRS)+(posMax-SegmQRS)*~(posMax<SegmQRS)) + find(signal(indice*(posMax<SegmQRS)+(posMax-SegmQRS)*~(posMax<SegmQRS):posMax) == min(signal(indice*(posMax<SegmQRS)+(posMax-SegmQRS)*~(posMax<SegmQRS):posMax)));
	posMin = posMin(end);
	beat(N) = round ((posMax+posMin)/2);		% Posicion QRS
	indice = beat(N) + EyeClosing;				% Desplazamiento para evitar sobredeteccion
	N = 2;										% Siguiente deteccion
	
	%/*-Barremos todo el registro hasta end-ventana-------------------*/
	while ((indice+ventana) <= nsamp) 
		if(max(signal(indice:indice+ventana)) > Umbral)												% Superiores a Umbral 
			posMax = find(signal(indice:indice+ventana) > max(signal(indice:indice+ventana))*0.4); 	% devuelve un vector
			relativo = diff(posMax);																% Si hay diferencia de 1 hay 1 pico, >1 hay 2 picos
			if(max(relativo) > 1)																	% Hay 2 picos
				relativo = vertcat(1, relativo);													% concatena 1 en horizontal (diff te saca 1 elemento)
				posicion = find(relativo == max(relativo));											% encontras sus posiciones dentro de posMax
				posicion = posicion(1);																% Si hay mas de 2 picos, te quedas con el primero
				%!!! VER SI REALMENTE FUNCIONA !!!
				relativo = ((ExpectedRR-EyeClosing)-posMax(posicion-1)) < (posMax(posicion)-(ExpectedRR-EyeClosing));% Comparo distancias relativas 1:Abajo de ExpectedRR/0:Arriba de ExpectedRR
				posMax = indice + find(signal(indice+posMax(posicion)*~relativo:indice+posMax((posicion-1)*relativo+end*~relativo)) == max(signal(indice+posMax(posicion)*~relativo:indice+posMax((posicion-1)*relativo+end*~relativo))));
			else																					% busco dentro del intervalo más próximo a ExpectedRR el maximo
				posMax = indice + find(signal(indice:posMax+indice) == max(signal(indice:posMax+indice)));
			end
            posMin = (posMax-SegmQRS) + find(signal(posMax-SegmQRS:posMax) == min(signal(posMax-SegmQRS:posMax)));
			posMin = posMin(end);							% Si es un vector no quedamos con el mas próximo a Q
            beat(N) = round((posMax+posMin)/2); 			% posicion QRS
			Umbral = signal(posMax)/2;						% Actualizamos al valor del último pico
			indice = beat(N) + EyeClosing;					% Desplazamiento para evitar sobredetección
			N = N + 1;										% Siguiente deteccion
			
			if(~rem(N,periodo))								% Cada cierto "periodo" se corrige ExpectedRR
				[ExpectedRR,periodo] = IntervaloRR(indice,signal,Fs, UmbralProm); % 
				ventana = ceil(1.5 * ExpectedRR);
			end
		else	%% no hay o me lo morfé
			Umbral = Umbral/3 * (max(signal(indice:indice+ventana))>UmbralProm);	% Superiores a Umbral 
			indice = indice + ceil(ExpectedRR/2);
		end
	end
	
	%/*-Tramo que falta desde indice hasta end------------------------*/
	if(max(signal(end-ventana:end)) > Umbral)													% Superiores a Umbral 
		indice = length(signal)-ventana;														% reubicamos indice
		posMax = find(signal(indice:indice+ventana) > max(signal(indice:indice+ventana))*0.4); 	% Devuelve un vector
		relativo = diff(posMax);																% Si hay diferencia de 1 hay 1 pico, >1 hay 2 picos
		if(max(relativo) > 1)																	% Hay 2 picos
			relativo = vertcat( 1, relativo);
			posicion = find(relativo == max(relativo));											% encontras sus posiciones dentro de posMax
			posicion = posicion(1);																% Si hay mas de 2 picos, te quedas con el primero
			relativo = ((ExpectedRR-EyeClosing)-posMax(posicion-1)) < (posMax(posicion)-(ExpectedRR-EyeClosing));% Comparo distancias relativas 1:Abajo de ExpectedRR/0:Arriba de ExpectedRR
			posMax = indice + find(signal(indice+posMax(posicion)*~relativo:indice+posMax((posicion-1)*relativo+end*~relativo)) == max(signal(indice+posMax(posicion)*~relativo:indice+posMax((posicion-1)*relativo+end*~relativo))));
		else																					% busco dentro del intervalo más próximo a ExpectedRR el maximo
			posMax = indice + find(signal(indice:posMax+indice) == max(signal(indice:posMax+indice)));
		end
		posMin = (posMax-SegmQRS) + find(signal(posMax-SegmQRS:posMax) == min(signal(posMax-SegmQRS:posMax)));
		posMin = posMin(end);							% Si es un vector no quedamos con el mas próximo a Q
		beat(N) = round((posMax+posMin)/2); 			% posicion del complejo QRS
		Umbral = signal(posMax)/2;						% Actualizamos al valor del último pico
		indice = beat(N) + EyeClosing;			% Desplazamiento para evitar sobredetección
		N = N + 1;
    end 
    
    %%2da parte: 
    %/*-Busqueda Local donde los complejos que no esten registrados --*/
    %/*-entre dos maximos. Nos ubicamos en la deteccion anterior donde*/
    %/*-donde pensamos que podría haber haber complejos---------------*/
	relativos = diff(beat);							% hacemos posiciones relativas entre maximos
	relativos2 = horzcat(beat(1),relativos(1:end-1)); % agregamos la posicion de 1ro
	perdidos = relativos./relativos2;						% proporcion entre una posicion y su anterior 
	perdidos = find(perdidos > Criterio_Perdido);			% posiciones de deteccion anterior al probable complejo
	
	for indice=1 :length(perdidos)
		Umbral = signal(beat(perdidos(indice)))*0.5;	% Tomamos el maximo de la deteccion 
		if (max(signal(beat(perdidos(indice)) + EyeClosing:beat(perdidos(indice)+1) - EyeClosing)) > Umbral);
		%if (max(signal(beat(perdidos(indice)) + EyeClosing:beat(perdidos(indice)+1))) > UmbralProm);
			posMax = (beat(perdidos(indice)) + EyeClosing) + find(signal(beat(perdidos(indice)) + EyeClosing:beat(perdidos(indice)+1) - EyeClosing) == max(signal(beat(perdidos(indice)) + EyeClosing:beat(perdidos(indice)+1) - EyeClosing)));
			posMax = posMax(end);							% si es un vector nos quedamos con el mas proximo a S
            posMin = (posMax - SegmQRS) + find(signal(posMax-SegmQRS:posMax) == min(signal(posMax-SegmQRS:posMax)));
			posMin = posMin(end);							% si es un vector nos quedamos con el mas proximo a Q
			beat(N) = round ((posMax+posMin)/2);	% posicion complejo QRS

			%/*2do miss-----------------------------------------------*/
			if ((beat(perdidos(indice)+2)-beat(perdidos(indice)+1))/(beat(perdidos(indice)+1)-beat(N)) > Criterio_Perdido)
				Umbral = signal(beat(perdidos(indice)))*0.5;	% Tomamos el maximo de la deteccion 
				if (max(signal(beat(perdidos(indice)+1) + EyeClosing:beat(perdidos(indice)+2)-EyeClosing)) > Umbral);
					posMax = (beat(perdidos(indice)+1) + EyeClosing) + find(signal(beat(perdidos(indice)+1) + EyeClosing:beat(perdidos(indice)+2)-EyeClosing) == max(signal(beat(perdidos(indice)+1) + EyeClosing:beat(perdidos(indice)+2)-EyeClosing)));
					posMax = posMax(end);							% si es un vector nos quedamos con el mas proximo a S
					posMin = (posMax - SegmQRS) + find(signal(posMax-SegmQRS:posMax) == min(signal(posMax-SegmQRS:posMax)));
					posMin = posMin(end);							% si es un vector nos quedamos con el mas proximo a Q
					beat(N+1) = round ((posMax+posMin)/2);	% posicion complejo QRS
					N = N + 1; 
				end
			end
			N = N + 1;
        end
	end
end

%/*------------------------------------------------------------------*/
%% IntervaloRR - Devuelve el intervalo ExpectedRR en un periodo de tiempo
%
%	usage:[ExpectedRR] = IntervaloRR(Posicion, Signal,FrecSamp)
%
%	Realiza la deteccion de maximos en periodo de 10 segundos(definible)
%	localiza los complejos y devuelve el promedio de los intervalos ExpectedRR
%	detectados. Para determinar los 10 segundos es necesaria la Fsamp
%
%/*------------------------------------------------------------------*/

function [ExpectedRR,periodo]=IntervaloRR(posicion, signal, FrecSamp, threshold)

	RRmin = round(0.66*FrecSamp);
	RRmax = round(1.25*FrecSamp);
	RRtipico = round(0.8*FrecSamp);
    periodo = 5;

	condicion = (posicion+FrecSamp*periodo) < length(signal);
	intervalo = signal(posicion:((posicion+FrecSamp*periodo)*condicion) + (length(signal)*~condicion)); % seleccionamos el  intervalo de interés
	maximos = find(intervalo > max(intervalo).*0.5);													% hallamos la posicion de los maximos
	ExpectedRR = diff(maximos);																			% distancia relativa entre maximos
    ExpectedRR = vertcat(1,ExpectedRR);																	% si detecta 1 solo maximo
	valores = find(ExpectedRR > max(ExpectedRR)*0.5);													% separamos los max contiguos
	ExpectedRR = ceil(mean(ExpectedRR(valores)));														% promedio de max y redondeo hacia arriba

	if(ExpectedRR > RRmax || ExpectedRR < RRmin) 														% verificamos que esta dentro del margen
		periodo = 2;
		ExpectedRR = RRtipico;
	end
	% Alternativa seria en vez de mirar para atras, mantener el valor que 
    % tenia anteriormente - Otra es dar el valor de Umbral como referencia
end

%[EOF]



