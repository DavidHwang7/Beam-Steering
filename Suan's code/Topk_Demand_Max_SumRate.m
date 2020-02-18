function [Design_DemandSum_Sort,Design_UserSet_Sort_Demand,Design_UserSet_Sort_User,utilization,throughput] = ...
    Topk_Demand_Max_SumRate(Report_Method,TopK,AP_Beam_CB,AP_lambda,AP_Pn,RB_number,RB_symbol,User_AoA,User_Distance,...
    Design_DemandSum_Sort,Design_UserSet_Sort_Demand,Design_UserSet_Sort_User,...
    BeamGain_pat,Stream_Number,xlsAlgle,xls_minSNR,xls_Efficiency)

%{
    建立SIR table，如果b1 b2 b3 b4的gain分別為10 20 30 40，b1的SIR就是10/(20+30+40)，
    再算capacity rate，table排序之後，Ck取4個beam的組合中選排序最高的作為服務的beam
%}
Bandwidth = 180000;% 180KHz per RB
if strcmp(Report_Method,'Actual_SINR')
    combi = nchoosek(1:size(BeamGain_pat,1),Stream_Number);%C7取4
    Table_rate = zeros(size(combi,1),size(BeamGain_pat,1)); 
    Table_SINR = zeros(size(combi,1),size(BeamGain_pat,1)); 
    for i=1:size(combi,1)
        Rx_gain = zeros(Stream_Number,1);
        interference = zeros(Stream_Number,Stream_Number);
        for index = 1:Stream_Number
            Rx_gain(index) = interp1(xlsAlgle,BeamGain_pat(combi(i,index),:),AP_Beam_CB(combi(i,index)));
            for index2 = 1:Stream_Number
                if index~=index2
                    interference(index,index2) = interp1(xlsAlgle,BeamGain_pat(combi(i,index2),:),AP_Beam_CB(combi(i,index)));
                end               
            end
            Table_SINR(i,combi(i,index)) = Rx_gain(index)/sum(interference(index,:));
            Table_SINR(i,combi(i,index)) = 10*log10( Table_SINR(i,combi(i,index)));
            Table_rate(i,combi(i,index)) = Bandwidth*log2(1+Table_SINR(i,combi(i,index)));
        end
    end
end

if strcmp(Report_Method,'Approximate_SINR')
    combi = nchoosek(1:size(BeamGain_pat,1),Stream_Number);%C7取4
    Table_rate = zeros(size(combi,1),size(BeamGain_pat,1));
    Table_SINR = zeros(size(combi,1),size(BeamGain_pat,1));
    for i=1:size(combi,1)
        test_beam = find(Table_rate(i,:));
        for beam_index=1:Stream_Number
            Rx_gain = interp1(xlsAlgle,BeamGain_pat(test_beam(beam_index),:),AP_Beam_CB(test_beam(beam_index)));
            interference_gain = zeros(size(BeamGain_pat,1),1);
            for interference_beam = 1:size(BeamGain_pat,1)
               if  interference_beam ~= test_beam(beam_index)
                   interference_gain(interference_beam) = interp1(xlsAlgle,BeamGain_pat(interference_beam,:),AP_Beam_CB(test_beam(beam_index)));
               end
            end
            if find(test_beam== find(interference_gain==max(interference_gain)))
                max_interference = max(interference_gain)+AP_Pn;
            else
                max_interference = AP_Pn;
            end
            Table_SINR(i,test_beam(beam_index)) = 10*log10(Rx_gain/max_interference);
            Table_rate(i,test_beam(beam_index)) = log2(1+Rx_gain/max_interference);
        end

    end

end

% sort combination according to sum-rate calculated by SIR
Table_sumrate = sum(Table_rate,2);
[M,Pola_Table_sumrate_sort] = sort(Table_sumrate,'descend');

demand = sum(Design_UserSet_Sort_Demand,2);
finding=0;
combination_beam = nchoosek(Design_DemandSum_Sort(1:TopK),Stream_Number);
beam = zeros(1,Stream_Number);
for k=1:size(Table_rate,1) %1->sum-rate最高
    tmp_beam = find(Table_rate(Pola_Table_sumrate_sort(k),:));
    for j=1:size(combination_beam,1)%find best beam from top k combination
        if size(find(demand(tmp_beam)),1)>0 && size(find(sort(combination_beam(j,:))==tmp_beam),2)==Stream_Number
            beam = tmp_beam;
            finding=1;
            break;
        end
    end
    if finding
        break;
    end

end

Design_transmit_bit_per_slot = 0;
Design_usedRB_per_slot = 0;%used RB
Design_RB_scheduler = zeros(size(BeamGain_pat,1),RB_number);%beam x RB
Design_remain_RB = repmat(RB_number,Stream_Number,1);%4*1個100
for i=1:Stream_Number %4個beam去找最多demand的user去填滿RB
    if (demand(beam(i))==0)%這個beam的user已經被服務完

    else
        j=1;%users are sorted,找較多demand的user
        while Design_remain_RB(i,1)>0
            if Design_UserSet_Sort_Demand(beam(i),j)==0 %這個beam沒有適合的user了

                break;
            else
                user = Design_UserSet_Sort_User(beam(i),j);
                if strcmp(Report_Method,'Actual_SINR')
                    SINR_dB = Actual_SINR(user,beam,i,xlsAlgle,User_AoA,AP_Pn,...
                        BeamGain_pat,Stream_Number,AP_lambda,User_Distance);
                elseif strcmp(Report_Method,'Approximate_SINR')
                    SINR_dB = Approximate_SINR(user,beam,i,xlsAlgle,User_AoA,AP_Pn,...
                        BeamGain_pat,Stream_Number,AP_lambda,User_Distance);
                elseif strcmp(Report_Method,'SIR')
                    SINR_dB = SIR(beam,i,xlsAlgle,AP_Beam_CB,AP_Pn,...
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
                    require_RB = ceil(Design_UserSet_Sort_Demand(beam(i),j)/(rate));%rate=transmit bit per RB
                    if require_RB > Design_remain_RB(i,1)
                        require_RB = Design_remain_RB(i,1);
                    end
                    Design_usedRB_per_slot = Design_usedRB_per_slot + require_RB; 
                    Design_RB_scheduler(beam(i),RB_number-Design_remain_RB(i,1)+1:RB_number-Design_remain_RB(i,1)+require_RB) =  user;
                   %{
                    if (Design_UserSet_Sort_Demand(beam(i),j)-require_RB*rate)<0
                        transmit_bits = Design_UserSet_Sort_Demand(beam(i),j);
                    else
                        transmit_bits = require_RB*rate;
                    end
                    Design_transmit_bit_per_slot = Design_transmit_bit_per_slot + transmit_bits;
                    Design_UserSet_Sort_Demand(beam(i),j) = Design_UserSet_Sort_Demand(beam(i),j)-transmit_bits;
                    %}
                    Design_remain_RB(i,1) = Design_remain_RB(i,1)-require_RB;
                end
                j=j+1;

            end

        end
    end
end
% real transmission
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
        user_index = find(Design_UserSet_Sort_User(beam(i),:)==user);
        if (Design_UserSet_Sort_Demand(beam(i),user_index)-RB_num*rate)<0
            transmit_bits = Design_UserSet_Sort_Demand(beam(i),user_index);
        else
            transmit_bits = RB_num*rate;
        end
        Design_transmit_bit_per_slot = Design_transmit_bit_per_slot + transmit_bits;
        Design_UserSet_Sort_Demand(beam(i),user_index) = Design_UserSet_Sort_Demand(beam(i),user_index)-transmit_bits;
        RB_index = RB_index + RB_num;
        if RB_index>size(Design_RB_scheduler,2)
            break;
        end
        user = Design_RB_scheduler(beam(i),RB_index);
    end
end
%---re-sort beams according to their sum of user's demands---

AP1_Beam_DemandSum = sum(Design_UserSet_Sort_Demand,2);
[M,Design_DemandSum_Sort] = sort(AP1_Beam_DemandSum,'descend');
if sum(M)==0%no user has demand

else
    for i=1:size(AP_Beam_CB,2)
        vec = find(Design_UserSet_Sort_User(i,:) ~=0);
        [Design_UserSet_Sort_Demand(i,vec),AP1_Beam_UserSet_Sort_Index(i,vec)] = sort(Design_UserSet_Sort_Demand(i,vec),'descend');
        Design_UserSet_Sort_User(i,vec) = Design_UserSet_Sort_User(i,AP1_Beam_UserSet_Sort_Index(i,vec));
    end
end
    
   

utilization = Design_usedRB_per_slot/(RB_number*Stream_Number);
throughput=  Design_transmit_bit_per_slot;
    
    
    