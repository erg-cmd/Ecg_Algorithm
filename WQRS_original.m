%% WQRS_original 
%
%
%

function [QRS,MULTI] = WQRSoriginal (ecg, HEADER, progreso, payload)
	
    %t=[1:HEADER.nsamp]/HEADER.freq; % Comentar! Solo debug!
    
    %*/LPF- notch en freq de linea------------------------------------*/
	num=[1 0 0 0 0 -2 0 0 0 0 1];						% num del filtro
	den=[1 -2 1];										% den del filtro
	ECG_filtered = filter(0.04*num,den,ecg); 			% Filtramos la senal y ajustamos para 0db
	
	
	% corregimos el atraso de muestras
    PhaseDelay = 11;
    medio = mean(ECG_filtered(1:HEADER.nsamp-PhaseDelay,:));
	for k=1:HEADER.nsig
		ECG_filtered(HEADER.nsamp-PhaseDelay:HEADER.nsamp,k) = medio(k);
	end
    
    %Agregado-1ra Variacion------------------------------------------------
    % necesario para el exito de la lengthtransform----------------------*/
    for k = 1:HEADER.nsig									%ver el efecto antes y despues de LT. Igual es solo para ver la señal.
		ECG_filtered(:,k) = (ECG_filtered(:,k) - HEADER.adczero(k)) ./ HEADER.gain(k);	% x esc de mV
	end
    %agregado-1ra Variacion-----------------------------------------------------*/
    
	%/*-Length Transform----------------------------------------------*/
	LT = lengthtransform (ECG_filtered,HEADER);
	minimo = min(LT);
	
	QRS = cell(HEADER.nsig,1);
	for n=1:HEADER.nsig
		LT(:,n) = LT(:,n) - minimo(n);
		QRS{n} = DecisionBlock_original(LT(:,n), HEADER.freq);		% retorna las detecciones, la posicion es en muestras
	end
	
	LTsum = sum(LT,2); 										% sumamos las señales por fila
	MULTI = DecisionBlock_original(LTsum,HEADER.freq);
	
	% Solo para debug - Comentar todo esto!!!!
	%figure(2);
	%plot(t,LTsum);title('Length Transform (sum)');grid on;hold on;xlabel('Time(sec)');
	%plot(t(MULTI),LTsum(MULTI),'c*','Markersize',5);
	%for index=1:length(MULTI)		%barrido de los puntos
	%	text(t(MULTI(index)),LTsum(MULTI(index)),[num2str(index)]);
	%end
end

%-----------------------------------------------------------------------
%% lengthtransform
% le das el tiempo y calcula la lt entre [ t-LTwindow+1 : t ]
% 
%
%-----------------------------------------------------------------------%
function lt = lengthtransform (ecg, header)
	MaxQRSw = 0.13;    						%/* maximum QRS width (130ms) */
	sps = header.freq;						% freq de muestreo
	LTwindow = round(sps * MaxQRSw);				%/* length transform window size */
	ltfs = 1.25*header.gain(1)*header.gain(1)*sps; 				% lt factor de escalamiento (ya lo pone ^2)
	dy = zeros(1,header.nsig);			% Creas un vector vacio del tamano de muestras
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
function [deteccion]=DecisionBlock_original(signal,Fs)
	EyeClosing = Fs * 0.25;			% Seteamos tiempo a 250ms
	%Tm = 100;						% Valor de Thres minimo
    Tm = mean(signal)/5;            % Valor de Thres minimo
	ExpectedPeriod = Fs * 2.5;		% Si no se encuentra QRS dentro de este tiempo se ajusta Threshold
	T0 = mean(signal(1:8*Fs));		% Promediar los primeros 8 segundos
	Ta = 3 * T0;					% Threshold Alto y Bajo
	t = 1;
	N = 1;
	deteccion(N) = 1;
	tiempo = 1;
	% Pregunta por un periodo de aprendizaje q resumo asi
	T1 = 2 * T0;	% 2 veces el Threshold
	
	% Loop principal
	while t<length(signal)-EyeClosing
		if signal(t) > T1 % posible QRS
			tiempo  = 0;		% me indica el tiempo que paso desde la ultima deteccion
			
			Vmax = max(signal(t+1:t+EyeClosing/2));
			Vmin = min(signal(t-EyeClosing/2:t-1));
			if(Vmax(end) > (0.5*Vmax(end) + Vmin(end)))
				posMax = find(signal(t+1:t+EyeClosing/2) == Vmax(end));
				posMin = find(signal(t-EyeClosing/2:t-1) == Vmin(end));
				deteccion(N) = t + (posMax(end)-posMin(end)+(EyeClosing/2)+1)/2;
				N = N + 1;
				Ta = Ta + (signal(t+posMax(end)) - Ta)/10 ; 	% ajuste de Threshold
				T1 = Ta / 3;				% ajuste Thres
				
				t = t + EyeClosing;
			end
			
		else
			tiempo = tiempo + 1; 
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

% [EOF]
