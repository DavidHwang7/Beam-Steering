function  [Design_DemandSum_Sort,Design_UserSet_Sort_Demand,Design_UserSet_Sort_User,utilization,throughput] = ...
    Topk_Demand_Max_Throughput(Report_Method,TopK,AP_Beam_CB,AP_lambda,AP_Pn,RB_number,RB_symbol,User_AoA,User_Distance,...
    Design_DemandSum_Sort,Design_UserSet_Sort_Demand,Design_UserSet_Sort_User,...
    BeamGain_pat,Stream_Number,xlsAlgle,xls_minSNR,xls_Efficiency)
%{
    Top k: Choose the beams whose demands are top k and then force search 4
    beams which tansmit more bits within a slot.
%}

Design_transmit_bit_per_slot = 0;
Design_usedRB_per_slot = 0;%used RB
demand = sum(Design_UserSet_Sort_Demand,2);

%% find the combination of beam
combination_beam = nchoosek(Design_DemandSum_Sort(1:TopK),Stream_Number);
tmp_requireRB = zeros(size(combination_beam,1),1);
tmp_tansmit_bit = zeros(size(combination_beam,1),1);
tmp_RB_scheduler = zeros(size(BeamGain_pat,1),RB_number,size(combination_beam,1));%beam x RB
for combination_index=1:size(combination_beam,1)
    tmp_beam = combination_beam(combination_index,:);
    Design_remain_RB = repmat(RB_number,Stream_Number,1);%4*1個100
    tmp_AP1_Beam_UserSet_Sort_Demand = Design_UserSet_Sort_Demand;
    for beam_index=1:Stream_Number %4個beam去找最多demand的user去填滿RB
        if (demand(tmp_beam(beam_index))==0)%這個beam的user已經被服務完

        else
            j=1;%users are sorted,找較多demand的user
            while Design_remain_RB(beam_index,1)>0 && j<=size(tmp_AP1_Beam_UserSet_Sort_Demand,2)
                if tmp_AP1_Beam_UserSet_Sort_Demand(tmp_beam(beam_index),j)==0 %這個beam沒有適合的user了%!!!!
                    break;
                else
                    user = Design_UserSet_Sort_User(tmp_beam(beam_index),j);
                   
                    if strcmp(Report_Method,'Actual_SINR')
                        SINR_dB = Actual_SINR(user,tmp_beam,beam_index,xlsAlgle,User_AoA,AP_Pn,...
                            BeamGain_pat,Stream_Number,AP_lambda,User_Distance);
                    elseif strcmp(Report_Method,'Approximate_SINR')
                        SINR_dB = Approximate_SINR(user,tmp_beam,beam_index,xlsAlgle,User_AoA,AP_Pn,...
                            BeamGain_pat,Stream_Number,AP_lambda,User_Distance);
                    elseif strcmp(Report_Method,'Approximate_SINR_two')
                        SINR_dB = Approximate_SINR_two(user,tmp_beam,beam_index,xlsAlgle,User_AoA,AP_Pn,...
                            BeamGain_pat,Stream_Number,AP_lambda,User_Distance);
                    elseif strcmp(Report_Method,'SIR')
                        SINR_dB = SIR(tmp_beam,beam_index,xlsAlgle,AP_Beam_CB,AP_Pn,...
                            BeamGain_pat,Stream_Number);
                    end
                    
                    if SINR_dB < min(xls_minSNR)
                        rate = 0;%transmit bit per RB
                    elseif SINR_dB > max(xls_minSNR)
                        rate = fix(RB_symbol*max(xls_Efficiency));
                    else
                        rate = fix(RB_symbol*xls_Efficiency( fix(interp1(xls_minSNR,find(xls_minSNR),SINR_dB))  ));
                    end 
                    %rate
                    %fix(RB_symbol*interpolation(SINR_dB,xls_minSNR,xls_Efficiency))

                    if rate==0%不要此user
                    else
                        require_RB = ceil(tmp_AP1_Beam_UserSet_Sort_Demand(tmp_beam(beam_index),j)/(rate));%rate=transmit bit per RB
                        if require_RB > Design_remain_RB(beam_index,1)
                            require_RB = Design_remain_RB(beam_index,1);
                        end
                        tmp_requireRB(combination_index) = tmp_requireRB(combination_index) + require_RB; 
                        tmp_RB_scheduler(tmp_beam(beam_index),RB_number-Design_remain_RB(beam_index,1)+1:RB_number-Design_remain_RB(beam_index,1)+require_RB,combination_index) =  user;
                    
                        if (tmp_AP1_Beam_UserSet_Sort_Demand(tmp_beam(beam_index),j)-require_RB*rate)<0
                            transmit_bits = tmp_AP1_Beam_UserSet_Sort_Demand(tmp_beam(beam_index),j);
                        else
                            transmit_bits = require_RB*rate;
                        end
                        tmp_tansmit_bit(combination_index) = tmp_tansmit_bit(combination_index) + transmit_bits;
                        tmp_AP1_Beam_UserSet_Sort_Demand(tmp_beam(beam_index),j) = tmp_AP1_Beam_UserSet_Sort_Demand(tmp_beam(beam_index),j)-transmit_bits;
                        Design_remain_RB(beam_index,1) = Design_remain_RB(beam_index,1)-require_RB;
                    end
                    j=j+1;

                end

            end
            
        end
    end

end

if size(find(tmp_tansmit_bit==max(tmp_tansmit_bit)),1)>1
    best_cimbination = find(tmp_requireRB==max(tmp_requireRB(find(tmp_tansmit_bit==max(tmp_tansmit_bit)))),1 );
    beam = combination_beam(best_cimbination,:);
else
    best_cimbination = find(tmp_tansmit_bit==max(tmp_tansmit_bit),1);
    beam = combination_beam(best_cimbination,:);
end

%% Really to transmit
Design_RB_scheduler = tmp_RB_scheduler(:,:,best_cimbination);
for i=1:Stream_Number
    RB_index = 1;
    user = Design_RB_scheduler(beam(i),RB_index);
    while user~=0 
        RB_num = size(find(Design_RB_scheduler(beam(i),:)==user),2);
        SINR_dB = Actual_SINR(user,beam,i,xlsAlgle,User_AoA,AP_Pn,...
            BeamGain_pat,Stream_Number,AP_lambda,User_Distance);
        if SINR_dB < min(xls_minSNR)
            rate = 0;%transmit bit per RB
        elseif SINR_dB > max(xls_minSNR)
            rate = fix(RB_symbol*max(xls_Efficiency));
        else
            rate = fix(RB_symbol*xls_Efficiency( fix(interp1(xls_minSNR,find(xls_minSNR),SINR_dB))  ));
        end 
        user_index = find(Design_UserSet_Sort_User==user);
        if (Design_UserSet_Sort_Demand(user_index)-RB_num*rate)<0
            transmit_bits = Design_UserSet_Sort_Demand(user_index);
        else
            transmit_bits = RB_num*rate;
        end
        Design_usedRB_per_slot = Design_usedRB_per_slot + RB_num; 
        Design_transmit_bit_per_slot = Design_transmit_bit_per_slot + transmit_bits;
        Design_UserSet_Sort_Demand(user_index) = Design_UserSet_Sort_Demand(user_index)-transmit_bits;
        RB_index = RB_index + RB_num;
        if RB_index>size(Design_RB_scheduler,2)
            break;
        end
        user = Design_RB_scheduler(beam(i),RB_index);
    end
end



%% Resort beam
Design_DemandSum = sum(Design_UserSet_Sort_Demand,2);
[M,Design_DemandSum_Sort] = sort(Design_DemandSum,'descend');
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