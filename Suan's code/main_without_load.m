function [] = main_without_load
%{
    use simulated beam pattern
%}
clear;
close all;
Bandwidth = 180000;% 180KHz per RB

%%  access excel of CQI for data rate and min SNR
xlsFile2 = './file/CQI_index.xlsx';
xls_Efficiency = xlsread(xlsFile2,1,'D2:D16').';%LTE
xls_minSNR = xlsread(xlsFile2,1,'E2:E16').';%LTE
%% Simulated pattern 
XLabel_Cell = {{'4','8','16','32','64'},{'30','15','10','5'},{'20','30','40','50'},{'2','4','6','8'}};
XLabel_Name = {'Number of antennas','Gap of beam','Number of User','Number of streams'};
num_case = 4;
case_size = [5 4 4 4];
case_antenna = [4 8 16 32 64];
case_beam = [30 15 10 5];
case_user = [20 30 40 50];
case_stream = [2 4 6 8];
num_of_design = 3;% Max D./SINR  , Topk/SINR  , Max rate <-----------------------------------
Plot_index = 0;
for num_case_index = 1:num_case
Antenna_Number = 4;
Beam_Gap = 30;
User_Number = 40;
Stream_Number = 4;
RB_number = 100;%number of RB per slot=100
RB_symbol = 84;%number of symbol per RB
AP_Pn_dB = -90;
AP_Pn = 10^(AP_Pn_dB/10);
AP_lambda = (3*10^8)/(20*10^6);%光速/20MHz
Packet_size = 3000;

%% setting for user:number,angle,distance,demand
design_throughput = zeros(case_size(num_case_index),num_of_design);
design_utilize = zeros(case_size(num_case_index),num_of_design);
for case_instance = 1:case_size(num_case_index)
    if num_case_index==1
        Antenna_Number =  case_antenna(case_instance);
    elseif num_case_index==2
        Antenna_Number = 8;
        Beam_Gap = case_beam(case_instance);
    elseif num_case_index==3
        User_Number = case_user(case_instance);
    else
        Stream_Number = case_stream(case_instance);
        Beam_Gap = 15;
    end
    BeamGain_pat = Gain_simulation(Antenna_Number,Beam_Gap);
    num_beam = (180/Beam_Gap)+1;
    AP_Beam_CB = 0:Beam_Gap:180;
    xlsAlgle = 0:5:180;
    User_Demand = zeros(1,User_Number);
    User_Demand(1,:) = Packet_size;
    
    test_size = 1; % calculate many times and average
    for test_index=1:test_size
    % Reset UE setting
    User_AoA = 0 + rand(1,User_Number)*180; 
    User_Distance = rand(1,User_Number)*150; % ~150m
    %% optimal beam for user
    AP1_Beam_UserSet_size = ceil(User_Number);
    AP1_Beam_UserSet = zeros(size(AP_Beam_CB,2),AP1_Beam_UserSet_size);
    AP1_Beam_UserSet_Demand = zeros(size(AP1_Beam_UserSet));
    [M,User_Beam]=min(abs(repmat(User_AoA,size(AP_Beam_CB,2),1)-repmat(AP_Beam_CB.',1,User_Number)),[],1);%!!!
    for i=1:size(AP_Beam_CB,2)
        AP1_Beam_UserSet(i,:) = [find(User_Beam==i) zeros(1,AP1_Beam_UserSet_size-size(find(User_Beam==i),2))];%!! if there are beams having zero UE will crash
        vec = find(AP1_Beam_UserSet(i,:) ~=0);
        AP1_Beam_UserSet_Demand(i,vec) = User_Demand(AP1_Beam_UserSet(i,vec));
    end
    %%  sort beams according to their sum of user's demands
    AP1_Beam_DemandSum = zeros(size(AP_Beam_CB,2),1);
    for i=1:size(AP_Beam_CB,2)
        copy_row_AP1_Beam_UserSet = AP1_Beam_UserSet(i,1:find(AP1_Beam_UserSet(i,:)==0,1)-1).';
        AP1_Beam_DemandSum(i) = sum(User_Demand(copy_row_AP1_Beam_UserSet),2);%bug:if beam has no user
    end
    [M,AP1_Beam_DemandSum_Sort] = sort(AP1_Beam_DemandSum,'descend');
    %%  sort users according to their demand for each user set of beams   
    AP1_Beam_UserSet_Sort_Demand = zeros(size(AP1_Beam_UserSet));
    AP1_Beam_UserSet_Sort_Index = zeros(size(AP1_Beam_UserSet));
    AP1_Beam_UserSet_Sort_User = zeros(size(AP1_Beam_UserSet));
    for i=1:size(AP_Beam_CB,2)
        vec = find(AP1_Beam_UserSet_Demand(i,:) ~=0);
        [AP1_Beam_UserSet_Sort_Demand(i,vec),AP1_Beam_UserSet_Sort_Index(i,vec)] = sort(AP1_Beam_UserSet_Demand(i,vec),'descend');
        AP1_Beam_UserSet_Sort_User(i,vec) = AP1_Beam_UserSet(i,AP1_Beam_UserSet_Sort_Index(i,vec));
    end
    %%
    index_of_design = 1;

    %% Max demand + Actual_SINR + Secondary UE
    %{
    index_of_design = index_of_design + 1;
    Design_UserSet_Sort_Demand = AP1_Beam_UserSet_Sort_Demand;
    Design_UserSet_Sort_User = AP1_Beam_UserSet_Sort_User;
    Design_DemandSum_Sort = AP1_Beam_DemandSum_Sort;
    [Design_DemandSum_Sort,Design_UserSet_Sort_Demand,Design_UserSet_Sort_User,design_utilize(case_index,index_of_design),design_throughput(case_index,index_of_design)] = ...
            Max_demand_secondary('Actual_SINR',AP_Beam_CB,AP_lambda,AP_Pn,RB_number,RB_symbol,User_AoA,User_Distance,...
            Design_DemandSum_Sort,Design_UserSet_Sort_Demand,Design_UserSet_Sort_User,Pola_Table,...
            xlsBeamGain_pos_pos,xlsBeamGain_pos_neg,xlsBeamGain_neg_neg,xlsBeamGain_neg_pos,xlsAlgle,xls_minSNR,xls_Efficiency);
    %}
    
    %% Max demand+Actual_SINR
    %{
      demand最多的四個beam，選最好的極化組合，各挑一個deamnd最高的user分配RB
    %}
    %index_of_design = index_of_design + 1;
    Design_UserSet_Sort_Demand = AP1_Beam_UserSet_Sort_Demand;
    Design_UserSet_Sort_User = AP1_Beam_UserSet_Sort_User;
    Design_DemandSum_Sort = AP1_Beam_DemandSum_Sort;
    [Design_DemandSum_Sort,Design_UserSet_Sort_Demand,Design_UserSet_Sort_User,utilize,throughput] = ...
            Max_demand('Actual_SINR',AP_Beam_CB,AP_lambda,AP_Pn,RB_number,RB_symbol,User_AoA,User_Distance,...
            Design_DemandSum_Sort,Design_UserSet_Sort_Demand,Design_UserSet_Sort_User ,...
            BeamGain_pat,Stream_Number,xlsAlgle,xls_minSNR,xls_Efficiency);
    design_utilize(case_instance,index_of_design) = design_utilize(case_instance,index_of_design) + utilize;
    design_throughput(case_instance,index_of_design) = design_throughput(case_instance,index_of_design) + throughput;
    %% Max Demand+SINR'
    %{
    index_of_design = index_of_design + 1; 
    Design_UserSet_Sort_Demand = AP1_Beam_UserSet_Sort_Demand;
    Design_UserSet_Sort_User = AP1_Beam_UserSet_Sort_User;
    Design_DemandSum_Sort = AP1_Beam_DemandSum_Sort;
    [Design_DemandSum_Sort,Design_UserSet_Sort_Demand,Design_UserSet_Sort_User,utilize,throughput] = ...
            Max_demand('Approximate_SINR',AP_Beam_CB,AP_lambda,AP_Pn,RB_number,RB_symbol,User_AoA,User_Distance,...
            Design_DemandSum_Sort,Design_UserSet_Sort_Demand,Design_UserSet_Sort_User ,...
            BeamGain_pat,Stream_Number,xlsAlgle,xls_minSNR,xls_Efficiency);
    design_utilize(case_index,index_of_design) = design_utilize(case_index,index_of_design) + utilize;
    design_throughput(case_index,index_of_design) = design_throughput(case_index,index_of_design) + throughput;
    %}

    %% Top k demand, max throughput + Actual_SINR + Secondary
    %{
    TopK = 5;
    index_of_design = index_of_design + 1;
    Design_UserSet_Sort_Demand = AP1_Beam_UserSet_Sort_Demand;
    Design_UserSet_Sort_User = AP1_Beam_UserSet_Sort_User;
    Design_DemandSum_Sort = AP1_Beam_DemandSum_Sort;
    [Design_DemandSum_Sort,Design_UserSet_Sort_Demand,Design_UserSet_Sort_User,design_utilize(case_index,index_of_design),design_throughput(case_index,index_of_design)] = ...
            Topk_Demand_Max_Throughput_secondary('Actual_SINR',TopK,AP_Beam_CB,AP_lambda,AP_Pn,RB_number,RB_symbol,User_AoA,User_Distance,...
        Design_DemandSum_Sort,Design_UserSet_Sort_Demand,Design_UserSet_Sort_User ,...
        BeamGain_pat,Stream_Number,xlsAlgle,xls_minSNR,xls_Efficiency);
    %}

    %% Top k demand, max throughput + Actual_SINR
    if num_case_index==4
        TopK = 9;
    else
        TopK = 6;
    end
    
    index_of_design = index_of_design + 1;
    Design_UserSet_Sort_Demand = AP1_Beam_UserSet_Sort_Demand;
    Design_UserSet_Sort_User = AP1_Beam_UserSet_Sort_User;
    Design_DemandSum_Sort = AP1_Beam_DemandSum_Sort;
    [Design_DemandSum_Sort,Design_UserSet_Sort_Demand,Design_UserSet_Sort_User,utilize,throughput] = ...
            Topk_Demand_Max_Throughput('Actual_SINR',TopK,AP_Beam_CB,AP_lambda,AP_Pn,RB_number,RB_symbol,User_AoA,User_Distance,...
        Design_DemandSum_Sort,Design_UserSet_Sort_Demand,Design_UserSet_Sort_User ,...
        BeamGain_pat,Stream_Number,xlsAlgle,xls_minSNR,xls_Efficiency);
    design_utilize(case_instance,index_of_design) = design_utilize(case_instance,index_of_design) + utilize;
    design_throughput(case_instance,index_of_design) = design_throughput(case_instance,index_of_design) + throughput;

    %% Top k demand, max throughput + SINR'
    %{
    TopK = 5;
    index_of_design = index_of_design + 1;
    Design_UserSet_Sort_Demand = AP1_Beam_UserSet_Sort_Demand;
    Design_UserSet_Sort_User = AP1_Beam_UserSet_Sort_User;
    Design_DemandSum_Sort = AP1_Beam_DemandSum_Sort;
    [Design_DemandSum_Sort,Design_UserSet_Sort_Demand,Design_UserSet_Sort_User,utilize,throughput] = ...
            Topk_Demand_Max_Throughput('Approximate_SINR',TopK,AP_Beam_CB,AP_lambda,AP_Pn,RB_number,RB_symbol,User_AoA,User_Distance,...
        Design_DemandSum_Sort,Design_UserSet_Sort_Demand,Design_UserSet_Sort_User ,...
        BeamGain_pat,Stream_Number,xlsAlgle,xls_minSNR,xls_Efficiency);
    design_utilize(case_index,index_of_design) = design_utilize(case_index,index_of_design) + utilize;
    design_throughput(case_index,index_of_design) = design_throughput(case_index,index_of_design) + throughput;
    %}
    %% Top k demand, max throughput + SIR
    %% Max Sumrate + Actual_SINR

    TopK = size(BeamGain_pat,1);
    index_of_design = index_of_design + 1;
    Design_UserSet_Sort_Demand = AP1_Beam_UserSet_Sort_Demand;
    Design_UserSet_Sort_User = AP1_Beam_UserSet_Sort_User;
    Design_DemandSum_Sort = AP1_Beam_DemandSum_Sort;
    [Design_DemandSum_Sort,Design_UserSet_Sort_Demand,Design_UserSet_Sort_User,utilize,throughput] = ...
        Topk_Demand_Max_SumRate('Actual_SINR',TopK,AP_Beam_CB,AP_lambda,AP_Pn,RB_number,RB_symbol,User_AoA,User_Distance,...
    Design_DemandSum_Sort,Design_UserSet_Sort_Demand,Design_UserSet_Sort_User ,...
    BeamGain_pat,Stream_Number,xlsAlgle,xls_minSNR,xls_Efficiency);
    design_utilize(case_instance,index_of_design) = design_utilize(case_instance,index_of_design) + utilize;
    design_throughput(case_instance,index_of_design) = design_throughput(case_instance,index_of_design) + throughput;
  %}  

    end
    design_utilize(case_instance,:) = design_utilize(case_instance,:)./test_size;
    design_throughput(case_instance,:) = design_throughput(case_instance,:)./test_size;
end

%% Save and plot figure
design_throughput = design_throughput ./ 0.0005;
design_utilize = design_utilize .* 100;
x_index = 1:case_size(num_case_index);
Plot_index = Plot_index + 1;
figure(Plot_index);
plot(x_index,design_throughput(:,1),'-xb',x_index,design_throughput(:,2),'-xm',x_index,design_throughput(:,3),'-xk');
set(gca,'XTick',x_index);
set(gca,'xticklabel',XLabel_Cell{num_case_index});
xlabel(XLabel_Name(num_case_index));
ylabel('Effective Throughput (bps)');
legend('Max D./SINR','Top5 D./SINR','Max Sumrate');
Plot_index = Plot_index + 1;
figure(Plot_index);
plot(x_index,design_utilize(:,1),'-xb',x_index,design_utilize(:,2),'-xm',x_index,design_utilize(:,3),'-xk');
set(gca,'XTick',x_index);
set(gca,'xticklabel',XLabel_Cell{num_case_index});
xlabel(XLabel_Name(num_case_index));
ylabel('Utilization (%)');
legend('Max D./SINR','Top5 D./SINR','Max Sumrate');
Figure_Legend_name = {'Max D./SINR','Top5 D./SINR','Max Sumrate'};
Figure_XLabel_Name = XLabel_Name(num_case_index);
Figure_XLabel_Cell = XLabel_Cell{num_case_index};
%{
filepath = './data/Throughput_';
filepath = join([filepath,int2str(num_case_index)]);
save(filepath,'design_throughput','Figure_Legend_name','Figure_XLabel_Name','Figure_XLabel_Cell');
filepath = './data/Utilization_';
filepath = join([filepath,int2str(num_case_index)]);
save(filepath,'design_utilize','Figure_Legend_name','Figure_XLabel_Name','Figure_XLabel_Cell');
%}
end


