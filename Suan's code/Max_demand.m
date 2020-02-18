function  [Design_DemandSum_Sort,Design_UserSet_Sort_Demand,Design_UserSet_Sort_User,utilization,throughput] = ...
    Max_demand(Report_Method,AP_Beam_CB,AP_lambda,AP_Pn,RB_number,RB_symbol,User_AoA,User_Distance,...
    Design_DemandSum_Sort,Design_UserSet_Sort_Demand,Design_UserSet_Sort_User,...
    BeamGain_pat,Stream_Number,xlsAlgle,xls_minSNR,xls_Efficiency)
%{
 Max demand+Report Method
%}
Design_transmit_bit_per_slot = 0;
Design_usedRB_per_slot = 0;%used RB
demand = sum(Design_UserSet_Sort_Demand,2);
beam = Design_DemandSum_Sort(1:Stream_Number);%取依照sum demand sort過後的前面4個beam
Design_RB_scheduler = zeros(size(AP_Beam_CB,2),RB_number);%beam x RB
Design_remain_RB = repmat(RB_number,Stream_Number,1);%4*1個100
for beam_index=1:Stream_Number %4個beam去找最多demand的user去填滿RB
    if (demand(beam(beam_index))==0)%這個beam的user已經被服務完

    else
        j=1;%users are sorted,找較多demand的user
        while Design_remain_RB(beam_index,1)>0
            if Design_UserSet_Sort_Demand(beam(beam_index),j)==0 %這個beam沒有適合的user了

                break;
            else
                user = Design_UserSet_Sort_User(beam(beam_index),j);
                if strcmp(Report_Method,'Actual_SINR')
                    SINR_dB = Actual_SINR(user,beam,beam_index,xlsAlgle,User_AoA,AP_Pn,...
                        BeamGain_pat,Stream_Number,AP_lambda,User_Distance);
                elseif strcmp(Report_Method,'Approximate_SINR')
                    SINR_dB = Approximate_SINR(user,beam,beam_index,xlsAlgle,User_AoA,AP_Pn,...
                        BeamGain_pat,AP_lambda,User_Distance);
                elseif strcmp(Report_Method,'Approximate_SINR_two')
                        SINR_dB = Approximate_SINR_two(user,beam,beam_index,xlsAlgle,User_AoA,AP_Pn,...
                        BeamGain_pat,AP_lambda,User_Distance);
                elseif strcmp(Report_Method,'SIR')
                    SINR_dB = SIR(beam,beam_index,xlsAlgle,AP_Beam_CB,AP_Pn,...
                        BeamGain_pat,Stream_Number);
                end
                if SINR_dB < min(xls_minSNR)
                    rate = 0;%transmit bit per RB
                elseif SINR_dB > max(xls_minSNR)
                    rate = fix(RB_symbol*max(xls_Efficiency));
                else
                    rate = fix(RB_symbol*xls_Efficiency( fix(interp1(xls_minSNR,find(xls_minSNR),SINR_dB))  ));
                end 
                if rate==0%不要此user
                else
                    require_RB = ceil(Design_UserSet_Sort_Demand(beam(beam_index),j)/(rate));%rate=transmit bit per RB
                    if require_RB > Design_remain_RB(beam_index,1)
                        require_RB = Design_remain_RB(beam_index,1);
                    end
                    Design_usedRB_per_slot = Design_usedRB_per_slot + require_RB; 
                    Design_RB_scheduler(beam(beam_index),RB_number-Design_remain_RB(beam_index,1)+1:RB_number-Design_remain_RB(beam_index,1)+require_RB) =  user;
                    %{
                    if (Design_UserSet_Sort_Demand(beam(i),j)-require_RB*rate)<0
                        transmit_bits = Design_UserSet_Sort_Demand(beam(i),j);
                    else
                        transmit_bits = require_RB*rate;
                    end
                    Design_transmit_bit_per_slot = Design_transmit_bit_per_slot + transmit_bits;
                    Design_UserSet_Sort_Demand(beam(i),j) = Design_UserSet_Sort_Demand(beam(i),j)-transmit_bits;
                    %}                    
                    Design_remain_RB(beam_index,1) = Design_remain_RB(beam_index,1)-require_RB;
                end
                j=j+1;

            end

        end
    end
end
% real transmission
for beam_index=1:Stream_Number
    RB_index = 1;
    user = Design_RB_scheduler(beam(beam_index),RB_index);
    while user~=0
        RB_num = size(find(Design_RB_scheduler(beam(beam_index),:)==user),2);
        SINR_dB = Actual_SINR(user,beam,beam_index,xlsAlgle,User_AoA,AP_Pn,...
            BeamGain_pat,Stream_Number,AP_lambda,User_Distance);
       
        if SINR_dB < min(xls_minSNR)
            rate = 0;%transmit bit per RB
        elseif SINR_dB > max(xls_minSNR)
            rate = fix(RB_symbol*max(xls_Efficiency));
        else
            rate = fix(RB_symbol*xls_Efficiency( fix(interp1(xls_minSNR,find(xls_minSNR),SINR_dB))  ));
        end 
        user_index = find(Design_UserSet_Sort_User(beam(beam_index),:)==user);
        if (Design_UserSet_Sort_Demand(beam(beam_index),user_index)-RB_num*rate)<0
            transmit_bits = Design_UserSet_Sort_Demand(beam(beam_index),user_index);
        else
            transmit_bits = RB_num*rate;
        end
        Design_transmit_bit_per_slot = Design_transmit_bit_per_slot + transmit_bits;
        Design_UserSet_Sort_Demand(beam(beam_index),user_index) = Design_UserSet_Sort_Demand(beam(beam_index),user_index)-transmit_bits;
        RB_index = RB_index + RB_num;
        if RB_index>size(Design_RB_scheduler,2)
            break;
        end
        user = Design_RB_scheduler(beam(beam_index),RB_index);
    end
end


%---re-sort beams according to their sum of user's demands---
AP1_Beam_DemandSum = sum(Design_UserSet_Sort_Demand,2);
[M,Design_DemandSum_Sort] = sort(AP1_Beam_DemandSum,'descend');
if sum(M)==0%no user has demand

else
    for beam_index=1:size(AP_Beam_CB,2)
        vec = find(Design_UserSet_Sort_User(beam_index,:) ~=0);
        [Design_UserSet_Sort_Demand(beam_index,vec),AP1_Beam_UserSet_Sort_Index(beam_index,vec)] = sort(Design_UserSet_Sort_Demand(beam_index,vec),'descend');
        Design_UserSet_Sort_User(beam_index,vec) = Design_UserSet_Sort_User(beam_index,AP1_Beam_UserSet_Sort_Index(beam_index,vec));
    end
end
    
  
utilization = Design_usedRB_per_slot/(RB_number*Stream_Number);
throughput=  Design_transmit_bit_per_slot;