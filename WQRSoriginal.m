%% WQRS_original 
%
% Modificacion - > 1ro interprete que trabaja con las unidades en mV
% pero trabaja directo con numero de cuentas

function [QRS,MULTI] = WQRSoriginal (ECG, HEADER, progreso, payload)
	
	% Comentar! Solo debug!
	t=[1:HEADER.nsamp]/HEADER.freq; 
    
    %~ */LPF- notch en freq de linea------------------------------------*/
	%~ num=[1 0 0 0 0 -2 0 0 0 0 1];						% num del filtro
	%~ den=[1 -2 1];										% den del filtro
	%~ ECG_filtered = filter(0.04*num,den,ecg); 			% Filtramos la senal y ajustamos para 0db
	%~ % corregimos el atraso de muestras
    %~ PhaseDelay = 11;
    %~ medio = mean(ECG_filtered(1:HEADER.nsamp-PhaseDelay,:));
	%~ for k=1:HEADER.nsig
		%~ ECG_filtered(HEADER.nsamp-PhaseDelay:HEADER.nsamp,k) = medio(k);
	%~ end
	%*/LPF de mi version----------------------------------------------*/
	% Un objeto filtro con 
	H = Filtro_FIR(HEADER.freq);
	%Aplicamos el filtro
	ECG_filtered=filter(H,ECG);	
	[phi,w]=phasedelay(H);
	PhaseDelay = round(phi(5));
	%Corregimos las muestras x el atrasao del filtro
	ECG_filtered=circshift(ECG_filtered,[-PhaseDelay 0]);
	
    %ajuste de muestras por delay del filtro
	medio = mean(ECG_filtered(1:HEADER.nsamp-PhaseDelay,:));
	for k=1:HEADER.nsig
		ECG_filtered(HEADER.nsamp-PhaseDelay:HEADER.nsamp,k) = medio(k);
	end
    
    %Agregado-1ra Variacion------------------------------------------------
    % necesario para el exito de la lengthtransform--DUDOSO-----------*/
    %~ for k = 1:HEADER.nsig
		%ECG_filtered(:,k) = (ECG_filtered(:,k) - HEADER.adczero(k)) ./ HEADER.gain(k);	% x esc de mV
		%~ ECG_filtered(:,k) = (ECG_filtered(:,k) - HEADER.adczero(k));
	%~ end
    %agregado-1ra Variacion-----------------------------------------------------*/
    
	%/*-Length Transform----------------------------------------------*/
	LT = lengthtransform (ECG_filtered,HEADER);
	minimo = min(LT); % impactada reduciendo el costo computacional
	
	QRS = cell(HEADER.nsig,1);
	for n=1:HEADER.nsig
		LT(:,n) = LT(:,n) - minimo(n);
		deteccion = DecisionBlock_original(LT(:,n), HEADER);% retorna las detecciones, la posicion es en muestras
		QRS{n} = deteccion;
		
		%Solo debug, desactivar esta parte!!!
		figure(n);
		plot(t,LT(:,n));title(['Length Transform ',num2str(n),' Leed']);grid on;hold on;xlabel('Time(sec)');
		plot(t(deteccion),LT(deteccion,n),'c*','Markersize',5);
		for index=1:length(deteccion)
			text(t(deteccion(index)),LT(deteccion(index),n),[num2str(index)]);
		end
	end
	
	LTsum = sum(LT,2); 										% sumamos las señales por fila
	MULTI = DecisionBlock_original(LTsum,HEADER);
	
	% Solo para debug - Comentar todo esto!!!!
	figure(3);
	plot(t,LTsum);title('Length Transform (SUM)');grid on;hold on;xlabel('Time(sec)');
	plot(t(MULTI),LTsum(MULTI),'c*','Markersize',5);
	for index=1:length(MULTI)		%barrido de los puntos
		text(t(MULTI(index)),LTsum(MULTI(index)),[num2str(index)]);
	end
end

%-----------------------------------------------------------------------
%% lengthtransform
% le das el tiempo y calcula la lt entre [ t-LTwindow+1 : t ]
% 
%
%-----------------------------------------------------------------------%
function lt = lengthtransform (ecg, header)
	MaxQRSw = 0.13;    								%/* maximum QRS width (130ms) */
	sps = header.freq;								%freq de muestreo
	LTwindow = round(sps * MaxQRSw);				%/* length transform window size */
	ltfs = 1.25*power(header.gain(1)*100000,2)/sps; 	% lt factor de escalamiento (ya lo pone ^2)
	dy = zeros(1,header.nsig);						% Creas un vector vacio del tamano de muestras
	lt = zeros(header.nsamp,header.nsig);			%Vector "lt" del tamaño de muestras llenado con 0
	%/*-------------------------------------------------------------*/
	for m = 1: header.nsamp-LTwindow+1       
		dy = sum(ecg(m:m+LTwindow-1,:),1);
		lt(m+LTwindow-1,:) = sqrt(ltfs+power(dy,2));  % hacemos la length transform
		dy = zeros(1,header.nsig);			% vacias el vector
	end
	
	%/*--------------------------------------------------------------*/
	media = mean(lt(LTwindow:end,:)); % 1er Arreglo dudoso
	for m=1:LTwindow-1
		lt(m,:) = media;
	end
end

%-----------------------------------------------------------------------
%% DecisionBlock
% Evaluamos la posicion de los complejos--------------------------------
%
%
%-----------------------------------------------------------------------
function [deteccion]=DecisionBlock_original(signal,header)
    Fs = header.freq;        % 
	EyeClosing = Fs * 0.25;		% Seteamos tiempo a 250ms
	Tm = .100 * header.gain(1);	% Definimos el Threshold minimo de 100uV, expresado mV.
	ExpectedPeriod = Fs * 2.5;	% Si no se encuentra QRS dentro de este tiempo se ajusta Threshold
	T0 = mean(signal(1:8*Fs));	% Promediar los primeros 8 segundos
	Ta = 3 * T0;				% Threshold Alto y Bajo
	t = 1;N = 1;
	Constante = 10;
	deteccion(N) = 1;
	tiempo = 1;
	% Pregunta por un periodo de aprendizaje q resumo asi
	T1 = 2 * T0;	% 2 veces el Threshold
	
	% Loop principal
	while t < length(signal)-EyeClosing
		if signal(t) > T1 		% posible QRS
			tiempo  = 0;		% me indica el tiempo que paso desde la ultima deteccion
			
			Vmax = max(signal(t+1:t+EyeClosing/2));
			Vmin = min(signal(t-EyeClosing/2:t-1));
			%if(Vmax(end) > (0.5*Vmax(end) + Vmin(end)))
			if(Vmax(end) > Vmin(end)+Constante)
				posMax = find(signal(t+1:t+EyeClosing/2) == Vmax(end));
				posMin = find(signal(t-EyeClosing/2:t-1) == Vmin(end));
				deteccion(N) = t + (posMax(end)-posMin(end)+(EyeClosing/2)+1)/2;
				N = N + 1;
				Ta = Ta + (signal(t+posMax(end)) - Ta)/10 ; 	% ajuste de Threshold
				T1 = Ta / 3;				% ajuste Thres
				t = t + EyeClosing;
			end
		else
			tiempo = tiempo + 5; 
			%Pasado el tiempo de QRS esperado se baja el Threshold
			if ((tiempo > ExpectedPeriod) && (Ta > Tm))
				Ta = Ta - Ta*.1;
				T1 = Ta/3;
			end
		end
		t = t + 5;
    end
    deteccion = uint32(deteccion);
end

%/*-------------------------------------------------------------------*/
%%FIR_Equiri Returns a discrete-time filter object.
% Equiripple Lowpass filter designed using the FIRPM function.

%/*-------------------------------------------------------------------*/
function Hd = Filtro_FIR(Fs)
	% All frequency values are in Hz.
	Fpass = 14;              % Passband Frequency
	Fstop = 70;              % Stopband Frequency
	Dpass = 0.057501127785;  % Passband Ripple% 
	Dstop = 0.00001;          % Stopband Attenuation
	dens  = 20;              % Density Factor
	
	% Calculate the order from the parameters using FIRPMORD.
	[N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
	
	% Calculate the coefficients using the FIRPM function.
	b  = firpm(N, Fo, Ao, W, {dens});
	Hd = dfilt.dffir(b);
end


% [EOF]
