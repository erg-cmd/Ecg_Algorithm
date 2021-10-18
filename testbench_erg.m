%%	testbench_erg - 
%	
%
%
%

function [QRS,B,MtxConf1,TP1_Ann,TP1_Detec,FN1_idx,FP1_idx,MtxConf2,TP2_Ann,TP2_Detec,FN2_idx,FP2_idx] = testbench_erg( num )
    archivo = ['/home/elias/Documentos/enMatlab/record/',num2str(num)];
    if exist(cat(2,archivo,'.dat'))
        ECGw = ECGwrapper( 'recording_name', archivo);
        ECG = ECGw.read_signal(1, ECGw.ECG_header.nsamp);
        [QRS, B] = WQRS10( ECG, ECGw.ECG_header,0,0);
        [MtxConf1,TP1_Ann,TP1_Detec,FN1_idx,FP1_idx] = bxb(ECGw.ECG_annotations,QRS{1,1},ECGw.ECG_header);
        [MtxConf2,TP2_Ann,TP2_Detec,FN2_idx,FP2_idx] = bxb(ECGw.ECG_annotations,QRS{2,1},ECGw.ECG_header);
        %plot_ecg_strip(ECGw);
    else
        sprintf('\tNo existe el archivo\n')
    end
end