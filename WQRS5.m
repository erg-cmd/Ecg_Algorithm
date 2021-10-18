%/*------------------------------------------------------------------*/
%% WQRS5 - 
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
function [QRS,MULTI]=WQRS5(ECG, HEADER,PROGRESO,PAYLOAD)
	
	t=[1:HEADER.nsamp]/HEADER.freq; % Comentar! Solo debug!
	H = Filtro_FIR(HEADER.freq);
	ECG_filtered=filter(H,ECG);	%Aplicamos el filtro

	[phi,w]=phasedelay(H);									
	PhaseDelay = round(phi(5));								
	ECG_filtered=circshift(ECG_filtered,[-PhaseDelay 0]);	%para la ver 2013
	
    %ajuste de muestras por delay del filtro
	medio = mean(ECG_filtered(1:HEADER.nsamp-PhaseDelay,:));
	for k=1:HEADER.nsig
		ECG_filtered(HEADER.nsamp-PhaseDelay:HEADER.nsamp,k) = medio(k);
	end
    
    for k = 1:HEADER.nsig									%ver el efecto antes y despues de LT. Igual es solo para ver la señal.
		%ECG_filtered(:,k) = (ECG_filtered(:,k) - HEADER.adczero(k)) ./ HEADER.gain(k);	% x esc de mV
		ECG_filtered(:,k) = (ECG_filtered(:,k) - HEADER.adczero(k));	% x esc de mV
    end
    
	LT = CurveLengthTrans(ECG_filtered,HEADER.nsamp,HEADER.nsig, HEADER.freq);	% retorna la Longitud de la Curva de cada señal
	minimo = min(LT);
	
	QRS = cell(HEADER.nsig,1);
	for n=1:HEADER.nsig
		LT(:,n) = LT(:,n) - minimo(n);
		deteccion = BloqueDecision(LT(:,n),HEADER.freq,HEADER.gain(n));
		QRS{n} = deteccion;		% retorna las detecciones, la posicion es en muestras
		
		%en caso de funcionar, desactivar esta parte!
		figure(n);
		plot(t,LT(:,n));title(['Length Transform ',num2str(n),'mmm']);grid on;hold on;xlabel('Time(sec)');
		plot(t(deteccion),LT(deteccion,n),'c*','Markersize',5);
		for index=1:length(deteccion)		%barrido de los puntos
			text(t(deteccion(index)),LT(deteccion(index),n),[num2str(index)]);
		end
	end
	
	
	LTsum = sum(LT,2); 														% sumamos las señales por fila
	MULTI = BloqueDecision(LTsum,HEADER.freq,HEADER.gain(1));				% retorna las detecciones, la posicion es en muestras
	
	%en caso de funcionar, desactivar esta parte!
	figure(3);
	plot(t,LTsum);title('Length Transform (sum)');grid on;hold on;xlabel('Time(sec)');
	plot(t(MULTI),LTsum(MULTI),'c*','Markersize',5);
	for index=1:length(MULTI)		%barrido de los puntos
		text(t(MULTI(index)),LTsum(MULTI(index)),[num2str(index)]);
	end
end


%/*-------------------------------------------------------------------*/
%%FIR_Equiri Returns a discrete-time filter object.
% Equiripple Lowpass filter designed using the FIRPM function.

%/*-------------------------------------------------------------------*/
function Hd = Filtro_FIR(Fs)
	% All frequency values are in Hz.
    %/*-Notch a 50Hz
    %Fpass = 15;              % Passband Frequency
    %Fstop = 50;              % Stopband Frequency
    %Dpass = 0.057501127785;  % Passband Ripple
    %Dstop = 0.0001;          % Stopband Attenuation
    %dens  = 20;              % Density Factor
    
    Fpass = 14;              % Passband Frequency
    Fstop = 70;              % Stopband Frequency
	Dpass = 0.057501127785;  % Passband Ripple% Calculate the order from the parameters using FIRPMORD.
	Dstop = 0.00001;          % Stopband Attenuation[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
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
	
	ventana = 0.013;						% tiempo de la ventana en segundos >>0.02
	c = 1.25;							% Factor de escalamiento en segundos >>0.9
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
%	usage: [beat] = BloqueDecision(signal,Fs,gain)	
% version5 - Dividida en 2 partes, determinación del nivel de umbral 
% para la posicion del jQRS (adaptativo) y estrategias de búsquedas
% dada la 1ra parte.

% Solamente se guardan las posiciones de los posibles QRS
%/*------------------------------------------------------------------*/

function [beat] = BloqueDecision(signal,Fs,gain)
	Ts = 1/Fs;
	periodo = 5; 						% Intervalo a muestrear
	EyeClosing = round(.3/Ts);			% periodo para evitar el sobredetección(ms)
	SegmQRS = round(.1/Ts);				% Segun bibliografia consultada 100ms
	nsamp=length(signal);				% numero de muestras que tiene la señal    
	indice=1; 							% barre e indica posicion en el estudio
	N = 1; 								% num de complejos detectados
	posMax = 1;posMin=1; 				% posicion del maximo/minimo
	Criterio_Perdido = 1.5;				% relacion entre deteccion actual y posterior en tiempo
	Tm = .100 * gain;						% Definimos el Threshold minimo de 100uV, expresado mV.
	DistMinima = 20;
	%/*---------------------------------------------------------------*/
	% Definimos el intervalo ExpectedRR para definir una ventana de busqueda
	UmbralProm = mean(signal(indice+(Fs*2):indice+(Fs*10))); 	%Definimos un Umbral promedio
	ExpectedRR = IntervaloRR(indice, signal, Fs, UmbralProm);	%posicion esperada QRS
	ventana = ceil(ExpectedRR);									%ventana de busqueda
	Umbral = max(signal(indice:indice+ventana))/2; 				% Fijamos nivel de umbral actual

	% 1ra Parte - Deteccion por umbral adaptativo
	%/*-1er Complejo--------------------------------------------------*/
	posMax = find(signal(indice:indice+ExpectedRR) == max(signal(indice:indice+ExpectedRR)));			% Dentro de 1 IntervaloRR debe estar 1 complejo
	posMin = (indice*(posMax<SegmQRS)+(posMax-SegmQRS)*~(posMax<SegmQRS)) + find(signal(indice*(posMax<SegmQRS)+(posMax-SegmQRS)*~(posMax<SegmQRS):posMax) == min(signal(indice*(posMax<SegmQRS)+(posMax-SegmQRS)*~(posMax<SegmQRS):posMax)));
	posMin = posMin(end);
	beat(N) = round ((posMax+posMin)/2);		% Posicion QRS
	indice = beat(N) + EyeClosing;				% Desplazamiento para evitar sobredeteccion
	N = 2;										% Siguiente deteccion
	
	%/*-Barremos todo el registro hasta end-ventana-------------------*/
	while ((indice+ventana) < nsamp) 
		if(max(signal(indice:indice+ventana)) > Umbral)												% Posible QRS 
			posMax = find(signal(indice:indice+ventana) > max(signal(indice:indice+ventana))*0.7); 	% devuelve un vector con la pos del Max
			posMax = indice + posMax(1); 						% me quedo con el primer max
			posMin = (posMax-SegmQRS) + find(signal(posMax-SegmQRS:posMax) == min(signal(posMax-SegmQRS:posMax)));
			posMin = posMin(end);								% Si es un vector no quedamos con el mas próximo a Q
			if(signal(posMax) - signal(posMin) > DistMinima)	% entonces lo consideramos un QRS
				beat(N) = round((posMax+posMin)/2); 			% posicion QRS
				Umbral = Umbral + (signal(posMax)-Umbral)/10;	% Actualizamos al valor del último pico
				indice = beat(N) + EyeClosing;					% Desplazamiento para evitar sobredetección
				N = N + 1;										% Siguiente deteccion
			else												% No es un QRS
				indice = indice + ceil(ExpectedRR/2);			% Nos corremos hasta posMax y buscamos de nuevo 
			end
			
			if(~rem(N,periodo))									% Cada cierto "periodo" se corrige ExpectedRR
				[ExpectedRR,periodo] = IntervaloRR(indice,signal,Fs, UmbralProm); % 
				ventana = ceil(ExpectedRR);
			end
		else	%no se detecto. Puede no haber o el Thres está alto!
			if (Umbral > Tm) % y si tiempo es > ExpectedRR
				Umbral = Umbral/3;							% Bajamos el umbral
			end
			indice = indice + ceil(ExpectedRR/5);			% nos corremos ExpectedRR/3
		end
	end
	
	%/*-Tramo que falta desde indice hasta end------------------------*/
	if(max(signal(end-ventana:end)) > Umbral)													% Ultimo QRS
		indice = length(signal)-ventana;														% reubicamos indice
		posMax = find(signal(indice:indice+ventana) > max(signal(indice:indice+ventana))*0.7); 	% Retorna las pos relativas del intervalo
		posMax = indice + posMax(end);
		posMin = (posMax-SegmQRS) + find(signal(posMax-SegmQRS:posMax) == min(signal(posMax-SegmQRS:posMax)));
		posMin = posMin(end);							% Si es un vector no quedamos con el mas próximo a Q
		if(signal(posMax) - signal(posMin) > DistMinima)% entonces lo consideramos un QRS
        
			beat(N) = round((posMax+posMin)/2); 			% posicion QRS
			Umbral = signal(posMax)/2;						% Actualizamos al valor del último pico
			indice = beat(N) + EyeClosing;					% Desplazamiento para evitar sobredetección
			N = N + 1;										% Siguiente deteccion
		end
    end 
    
    %%2da parte: 
    %/*-Busqueda Local donde los complejos que no esten registrados --*/
    %/*-entre dos maximos. Nos ubicamos en la deteccion anterior donde*/
    %/*-donde pensamos que podría haber haber complejos---------------*/
	relativos = diff(beat);									% hacemos posiciones relativas entre maximos
	relativos2 = horzcat(beat(1),relativos(1:end-1)); 		% agregamos la posicion de 1ro
	perdidos = relativos./relativos2;						% proporcion entre una posicion y su anterior 
	perdidos = find(perdidos > Criterio_Perdido);			% posiciones de deteccion anterior al probable complejo
	
	for indice=1 :length(perdidos)
		posMax = (beat(perdidos(indice)) + EyeClosing) + find(signal(beat(perdidos(indice)) + EyeClosing:beat(perdidos(indice)+1) - EyeClosing) == max(signal(beat(perdidos(indice)) + EyeClosing:beat(perdidos(indice)+1) - EyeClosing)));
		if (~isempty(posMax))
			posMax = posMax(end);% si es un vector nos quedamos con el mas proximo a S
			posMin = (posMax - SegmQRS) + find(signal(posMax-SegmQRS:posMax) == min(signal(posMax-SegmQRS:posMax)));
			posMin = posMin(end);										% si es un vector nos quedamos con el mas proximo a Q
			if(signal(posMax) - signal(posMin) > DistMinima)% entonces lo consideramos un QRS
			%if(signal(posMax)-signal(posMin) > ConsOcular)
                beat(N) = round ((posMax+posMin)/2);						% posicion complejo QRS

				%/*2do miss-----------------------------------------------*/
				if ((beat(perdidos(indice)+2)-beat(perdidos(indice)+1))/(beat(perdidos(indice)+1)-beat(N)) > Criterio_Perdido)
					posMax = (beat(perdidos(indice)+1) + EyeClosing) + find(signal(beat(perdidos(indice)+1) + EyeClosing:beat(perdidos(indice)+2)-EyeClosing) == max(signal(beat(perdidos(indice)+1) + EyeClosing:beat(perdidos(indice)+2)-EyeClosing)));
					if (~isempty(posMax))
						posMax = posMax(end);% si es un vector nos quedamos con el mas proximo a S
						posMin = (posMax - SegmQRS) + find(signal(posMax-SegmQRS:posMax) == min(signal(posMax-SegmQRS:posMax)));
						posMin = posMin(end);							% si es un vector nos quedamos con el mas proximo a Q
						if(signal(posMax) - signal(posMin) > DistMinima)% entonces lo consideramos un QRS
						%if(signal(posMax)-signal(posMin) > ConsOcular)
                            beat(N+1) = round ((posMax+posMin)/2);	% posicion complejo QRS
							N = N + 1; 
						end
					end
				end
				N = N + 1;
			end
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
