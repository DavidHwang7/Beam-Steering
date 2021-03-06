function [SINR_dB] = Approximate_SINR(user,beam,beam_index,xlsAlgle,User_AoA,AP_Pn,...
    BeamGain_pat,AP_lambda,User_Distance)

Rx_gain = interp1(xlsAlgle,BeamGain_pat(beam(beam_index),:),User_AoA(user));
Rx_gain = (AP_lambda/(4*pi*User_Distance(user)))^2 * Rx_gain;
interference_gain = zeros(size(BeamGain_pat,1),1);
for interference_beam = 1:size(BeamGain_pat,1)
   if  interference_beam ~= beam(beam_index)
       interference_gain(interference_beam) = interp1(xlsAlgle,BeamGain_pat(interference_beam,:),User_AoA(user));
       %Path Loss
       interference_gain(interference_beam) = (AP_lambda/(4*pi*User_Distance(user)))^2 * interference_gain(interference_beam);
   end
end
if find(beam == find(interference_gain==max(interference_gain)))
    max_interference = max(interference_gain)+AP_Pn;
else
    max_interference = AP_Pn;
end
SINR = Rx_gain/max_interference;
SINR_dB = 10*log10(SINR);

