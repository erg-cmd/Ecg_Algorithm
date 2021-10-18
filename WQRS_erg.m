%% WQRS_erg
% Implementacion modesta del algoritmo WQRS descripto en "Robust-OpenSource
% Algorithm to detect the onset and duration complexes"
% para prueba trabajamos con  el estudio 200m 2^11 DC:1024 200 cuentasx mV


function [latidos]=WQRS_erg()
	load 100m;
	% hay que leer headers! para ver las variables del entorno fisico
	Fs =360;	%Hz
	Ts = 1/Fs;
	y = val(2,1:100000);
	t=[1:length(y)]/Fs; %x Esc de Tiempo
	y1=filter(FIR_Equiri,y);
	[phi,w]=phasedelay(FIR_Equiri);
	PhaseDelay = round(phi(5));
	y1=circshift(y1,-PhaseDelay,2);
	y1(end-PhaseDelay:end) = mean (y1);
	y1 = (y1-1024)/200;	% x esc de mV
	y2=LongCurva(10,y1,1/0.08);
	%latidos = DecisionBlock(y2);
    y2 = y2 - min(y2);
    latidos = DecisionBlock3(y2);
	figure(2);
	plot(t,y2);grid on; hold on ;xlabel('Time(sec)');
	for index=1:length(latidos)		%barrido de los pnutos
	plot(t(latidos(index).posicion),y2(latidos(index).posicion),'c*');
    text(t(latidos(index).posicion),y2(latidos(index).posicion),[num2str(index)]);
    end
    
   %figure(2);
    %plot(t,y1);grid on; hold on ;xlabel('Time(sec)');
	%for index=1:length(latidos)		%barrido de los pnutos
	%plot(t(latidos(index).posicion),y1(latidos(index).posicion),'c*');
   %text(t(latidos(index).posicion),y1(latidos(index).posicion),[num2str(index)]);
    %end
    
end
% [EOF]
