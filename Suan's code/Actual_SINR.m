function [SINR_dB] = Actual_SINR(user,beam,beam_index,xlsAlgle,User_AoA,AP_Pn,...
    BeamGain_pat,Stream_Number,AP_lambda,User_Distance)
Tx_Watt = 10^(18/10);
Rx_gain = interp1(xlsAlgle,BeamGain_pat(beam(beam_index),:),User_AoA(user));
Rx_gain = (AP_lambda/(4*pi*User_Distance(user)))^2 * Rx_gain * Tx_Watt;
interference_gain = zeros(1,Stream_Number);
for a=1:Stream_Number%interference
    if a~=beam_index
        interference_gain(1,a) = interp1(xlsAlgle,BeamGain_pat(beam(a),:),User_AoA(user));
        interference_gain(1,a) = (AP_lambda/(4*pi*User_Distance(user)))^2 * interference_gain(1,a) * Tx_Watt;
    end
end
SINR = Rx_gain/(sum(interference_gain)+AP_Pn);
SINR_dB = 10*log10(SINR);





