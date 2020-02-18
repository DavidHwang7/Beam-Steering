close all;
addpath ./ewa_function; 
d=0.5; % the number of wavelengths between antenna array
%ph0=60; % the initialized degree of main lobe direction
N=16; % the number of antenna in antenna array for an mainlobe
num_sample=400;
Xtx1=250; % transimtter1 coordinate in x-axis (1 meter for 1 unit)
Ytx1=500; % transimtter1 coordinate in y-axis
Xtx2=750; % transimtter2 coordinate in x-axis (1 meter for 1 unit)
Ytx2=500; % transimtter2 coordinate in y-axis
bandwidth = 20; % MHz, for calculating data rate


%Ptx=(10^1.8)/1000; % transmitter's transmit power (watt)
%Ptx_dBm = 10*log10(Ptx*1000); % transmit power (dBm) = 18 dBm
Wavelength=10e-3; % millimeter wavelength (meter)
Noise_dBm = -90; % Noise power (dBm) = -90 dBm
phi = (0 : num_sample) * pi / num_sample;     % equally-spaced over [0,pi]       
psi = 2 * pi * d * cos(phi);


%% Experiment 1: 固定50 users , 固定4個beam, codebook分別間隔 90、60、45、30、15、10、5(度)
num_Rx=50;
num_beam=4;
num_time_blocks=100;
%codebook_interval = [90 60 45 30 15 10 5];
codebook_interval = [60 45 30 15 10 5];
[~, num_codebook] = size(codebook_interval);
Xrx = randi([0,1000],1,50);
Yrx = randi([0,1000],1,50);
%[~, num_Rx1] = size(Xrx1);

% 每個user set都run 100~1000次
dist1=sqrt((Xrx-Xtx1).^2+(Yrx-Ytx1).^2); % the distance between transmitter1 and receivers
dist2=sqrt((Xrx-Xtx2).^2+(Yrx-Ytx2).^2); % the distance between transmitter2 and receivers
[~ ,prefer_list1] = sort(dist1);
[~ ,prefer_list2] = sort(dist2);


side1=Xrx-Xtx1;
side2=Xrx-Xtx2;

RxAngle1 = zeros(1, num_Rx);
RxAngle2 = zeros(1, num_Rx);
for i=1:num_Rx
    RxAngle1(i)=atan((Yrx(i)-Ytx1)/(Xrx(i)-Xtx1))+pi/2;
    RxAngle2(i)=atan((Yrx(i)-Ytx2)/(Xrx(1)-Xtx2))+pi/2;

end
Prx_dBm = zeros(num_codebook, num_Rx, num_time_blocks);
Prx_dBm2 = zeros(num_codebook, num_Rx, num_time_blocks);
SINR = zeros(num_codebook, num_Rx, num_time_blocks);
SINR2 = zeros(num_codebook, num_Rx, num_time_blocks);
%Rate = zeros(num_codebook, num_Rx, num_time_blocks);
Rate2 = zeros(num_codebook, num_Rx, num_time_blocks);
Rate_avg1=zeros(1,num_codebook);

for n=1:num_codebook
    degree_interval = codebook_interval(n);
    current_codebook_size = 180/degree_interval+1;
    a_current = zeros(current_codebook_size, N);
    A_current=zeros(current_codebook_size, num_sample+1);
    g_current=zeros(current_codebook_size, num_sample+1);
    for i=1:current_codebook_size
        ph0=(i-1)*degree_interval;
        a_current(i,:)= uniform(d, ph0, N);
        A_current(i,:) = dtft(a_current(i,:), -psi);    % array factor, note dsft(a,psi)=dtft(a,-psi)
        g_current(i,:) = abs(A_current(i,:)).^2;  % power gain
    end
    Ptx_total=(10^1.8)/1000; % transmitter's transmit power (watt)
    Ptx_beam = Ptx_total;
    Ptx_dBm = 10*log10(Ptx_beam*1000); % transmit power (dBm) = 18 dBm
    
    %Ptx_beam=Ptx_total/num_beam;
    %Ptx_dBm = 10*log10(Ptx_beam*1000); % transmit power (dBm) = 18 dBm
    degree_index_rx1=zeros(1,num_Rx);
    degree_index_rx2=zeros(1,num_Rx);
    maxBeam_index_rx1=zeros(1,num_Rx);
    maxBeam_index_rx2=zeros(1,num_Rx);
    %beam_chosen = zeros(1,4);
    
    for i=1:num_Rx
        [~ , degree_index_rx1(i)]=min(abs(phi(:)-abs(RxAngle1(i))));
        [~ , degree_index_rx2(i)]=min(abs(phi(:)-abs(RxAngle2(i))));
        [~ , maxBeam_index_rx1(i)]=max(g_current(:, degree_index_rx1(i)));
        [~ , maxBeam_index_rx2(i)]=max(g_current(:, degree_index_rx2(i)));
    end
    
    user_service_time=zeros(1,num_Rx);
    users_served_run=1;
    prefer_index1=1;
    prefer_index2=1;
    
    %%% 跑 num_time_blocks 個 time blocks，每個time block共服務8個user 
    for t=1:num_time_blocks
    %% Tx1,Tx2在一time block中各選4個不同user，並分配不同的beam
        user_chosen1=zeros(1,num_beam);
        user_chosen2=zeros(1,num_beam);
        beam_chosen1=zeros(1,num_beam);
        beam_chosen2=zeros(1,num_beam);
        
        i1=1;
        i2=1;

        while i1<=num_beam || i2<=num_beam
            
           min_served_time = min(user_service_time);
           if min_served_time==users_served_run
               users_served_run = users_served_run+1;
               prefer_index1=1;
               prefer_index2=1;
           end
           
           if i1<=num_beam && isempty(find(user_chosen2==prefer_list1(prefer_index1), 1)) && user_service_time(prefer_list1(prefer_index1))< users_served_run
                user_chosen1(i1)=prefer_list1(prefer_index1);
                beam_temp = maxBeam_index_rx1(user_chosen1(i1));
                if side1(user_chosen1(i1))<0 && beam_temp~=1 && beam_temp~=current_codebook_size
                    beam_temp = -beam_temp;
                end
                
                jj=2;
                g_temp=g_current(:, degree_index_rx1(user_chosen1(i1)));
                [~, beam_sort]=sort(g_temp,'descend');
                while isempty(find(beam_chosen1==beam_temp,1))==false
                    beam_temp = beam_sort(jj);
                    if side1(user_chosen1(i1))<0 && beam_temp~=1 && beam_temp~=current_codebook_size
                        beam_temp = -beam_temp;
                    end
                    jj =jj+1;
                end
                beam_chosen1(i1)=beam_temp;
                prefer_index1 = prefer_index1+1;
                user_service_time(user_chosen1(i1)) = user_service_time(user_chosen1(i1))+1;
                i1=i1+1;
           elseif i1<=num_beam 
               prefer_index1 = prefer_index1+1;
           end
           
           if i2<=num_beam && isempty(find(user_chosen1==prefer_list2(prefer_index2), 1)) && user_service_time(prefer_list2(prefer_index2))< users_served_run
                user_chosen2(i2)=prefer_list2(prefer_index2);
                beam_temp = maxBeam_index_rx2(user_chosen2(i2));
                if side2(user_chosen2(i2))<0 && beam_temp~=1 && beam_temp~=current_codebook_size
                    beam_temp = -beam_temp;
                end
                
                jj=2;
                g_temp=g_current(:, degree_index_rx2(user_chosen2(i2)));
                [~, beam_sort]=sort(g_temp,'descend');
                while isempty(find(beam_chosen2==beam_temp,1))==false
                    beam_temp=beam_sort(jj);
                    if side2(user_chosen2(i2))<0 && beam_temp~=1 && beam_temp~=current_codebook_size
                        beam_temp = -beam_temp;
                    end
                    jj=jj+1;
                end
                beam_chosen2(i2)=beam_temp;
                prefer_index2 = prefer_index2+1;
                user_service_time(user_chosen2(i2)) = user_service_time(user_chosen2(i2))+1;
                i2 = i2+1;
           elseif i2<=num_beam 
               prefer_index2 = prefer_index2+1;
           end
        end
    %% 計算receive power和interference
        i1=1;
        i2=1;
        while i1<=num_beam || i2<=num_beam
            %%%%% Tx1 %%%%%
            if beam_chosen1(i1)<0
                beam_current=-beam_chosen1(i1);
            else
                beam_current=beam_chosen1(i1);
            end
            TXAntennaGain=10*log10(g_current(beam_current,degree_index_rx1(user_chosen1(i1)))); %dB
            RXAntennaGain=1;%dB
            M1 = Wavelength / (4 * pi * dist1(user_chosen1(i1)));
            M2 = Wavelength / (4 * pi * dist2(user_chosen1(i1)));
            Prx_dBm(n,user_chosen1(i1),t) = Ptx_dBm + TXAntennaGain + RXAntennaGain + (20*log10(M1)); %Friis equation
            I_total1=0; % interference
            j1=1;
            j2=1;
            while j1<=num_beam || j2<=num_beam
                if j1 ~= i1 || side1(user_chosen1(i1))*side1(user_chosen1(j1))>=0 %%判斷是否為同一個user＆判斷是否在同一側，同側的話side相乘會>=0
                    if beam_chosen1(j1)<0
                        beam_current=-beam_chosen1(j1);
                    else
                        beam_current=beam_chosen1(j1);
                    end
                    TXAntennaGain=10*log10(g_current(beam_current,degree_index_rx1(user_chosen1(i1)))); %dB
                    I_dBm = Ptx_dBm + TXAntennaGain + RXAntennaGain + (20*log10(M1));
                    I_total1 = I_total1 + 10^(I_dBm/10);
                end
                if (side1(user_chosen1(i1))<=0 && side2(user_chosen2(j2))<=0) || (side1(user_chosen1(i1))>=0 && side2(user_chosen2(j2))<=0)
                    if beam_chosen2(j2)<0
                        beam_current=-beam_chosen2(j2);
                    else
                        beam_current=beam_chosen2(j2);
                    end
                    TXAntennaGain=10*log10(g_current(beam_current,degree_index_rx2(user_chosen1(i1)))); %dB
                    I_dBm = Ptx_dBm + TXAntennaGain + RXAntennaGain + (20*log10(M2));
                    I_total1 = I_total1 + 10^(I_dBm/10);
                end
                j1 = j1+1;
                j2 = j2+1;
            end
            
            Prx1 = 10^(Prx_dBm(n,user_chosen1(i1),t)/10);
            SINR(n,user_chosen1(i1),t) = Prx1 / (I_total1 + 10^(Noise_dBm/10));
            i1 = i1+1;
            
            %%%%% Tx2 %%%%%
            if beam_chosen2(i2)<0
                beam_current=-beam_chosen2(i2);
            else
                beam_current=beam_chosen2(i2);
            end
            TXAntennaGain=10*log10(g_current(beam_current,degree_index_rx2(user_chosen2(i2)))); %dB
            RXAntennaGain=1;%dB
            M1 = Wavelength / (4 * pi * dist1(user_chosen2(i2)));
            M2 = Wavelength / (4 * pi * dist2(user_chosen2(i2)));
            Prx_dBm(n,user_chosen2(i2),t) = Ptx_dBm + TXAntennaGain + RXAntennaGain + (20*log10(M2)); %Friis equation
            I_total2=0; % interference
            j1=1;
            j2=1;
            while j1<=num_beam || j2<=num_beam
                if j2 ~= i2 || side2(user_chosen2(i2))*side2(user_chosen2(j2))>=0 %%判斷是否為同一個user＆判斷是否在同一側，同側的話side相乘會>=0
                    if beam_chosen2(j2)<0
                        beam_current=-beam_chosen2(j2);
                    else
                        beam_current=beam_chosen2(j2);
                    end
                    TXAntennaGain=10*log10(g_current(beam_current,degree_index_rx2(user_chosen2(i2)))); %dB
                    I_dBm = Ptx_dBm + TXAntennaGain + RXAntennaGain + (20*log10(M2));
                    I_total2 = I_total2 + 10^(I_dBm/10);
                end
                if (side2(user_chosen2(i2))>=0 && side1(user_chosen1(j1))>=0) || (side2(user_chosen2(i2))<=0 && side1(user_chosen1(j1))>=0)
                    if beam_chosen1(j1)<0
                        beam_current=-beam_chosen1(j1);
                    else
                        beam_current=beam_chosen1(j1);
                    end
                    TXAntennaGain=10*log10(g_current(beam_current,degree_index_rx1(user_chosen2(i2)))); %dB
                    I_dBm = Ptx_dBm + TXAntennaGain + RXAntennaGain + (20*log10(M1));
                    I_total2 = I_total2 + 10^(I_dBm/10);
                end
                j1 = j1+1;
                j2 = j2+1;
            end
            
            Prx2 = 10^(Prx_dBm(n,user_chosen2(i2),t)/10);
            SINR(n,user_chosen2(i2),t) = Prx2 / (I_total2 + 10^(Noise_dBm/10));
            i2 = i2+1;
            
            %Rate(n,user_chosen1(i1),t) = bandwidth * log2(1 + SINR1(n,i));
            %SINR_dB(n,i)=10*log10(SINR(n,i));
            %fprintf('Prx%d(dBm) = %.4f\n',i,Prx_dBm(i));
            %fprintf('Rx%d''s SINR(dB) = %.4f\n\n',i,SINR_dB);
            %disp(Prx_dBm);
            %fprintf('Tx Antenna Gain toward Rx%d = %.4f\n\n',i,TXAntennaGain);
            %disp(TXAntennaGain);
            
            
        end
        
    end
    
    %SINR_avg_dB(n) = 10*log10(mean(SINR(n,:)));
    %Rate_avg1(n)=mean(Rate1(n,:));
    clear a_current A_current g_current;
end

Rate = bandwidth * log2(1 + SINR);
Total_Mb_per_user = sum(Rate,3);
Rate_average_codebook = sum(Total_Mb_per_user,2)/num_time_blocks;
figure;
plot(codebook_interval, Rate_average_codebook, '-*r');
xlabel('codebook interval: 60 45 30 15 10 5');
ylabel('Data rate (Mbps)');    
    
clear Prx_dBm;
    

%{
%% Experiment 2: codebook固定間隔15度, 20~80 users

num_users = [20 30 40 50 60 70 80];
num_users_set = size(num_users,2);
degree_interval = 15;
codebook_size = 180/degree_interval + 1;
a = zeros(codebook_size, N);
A = zeros(codebook_size, num_sample+1);
g = zeros(codebook_size, num_sample+1);
for i=1:codebook_size
    ph0=(i-1)*degree_interval;
    a(i,:)= uniform(d, ph0, N);
    A(i,:) = dtft(a(i,:), -psi);    % array factor, note dsft(a,psi)=dtft(a,-psi)
    g(i,:) = abs(A(i,:)).^2;  % power gain
end

Rate_avg = zeros(1,num_users_set);

for n=1:size(num_users,2)
    num_Rx = num_users(n);
    Prx_dBm = zeros(1, num_Rx);
    SINR = zeros(1, num_Rx);
    Rate = zeros(1, num_Rx);
    
    Xrx = randi([100,500],1,num_Rx);
    rand_index = randperm(num_Rx);
    draw_rand_index = rand_index(1:num_Rx/2);
    Xrx(draw_rand_index) = -1 * Xrx(draw_rand_index);
    Yrx = randi([100,500],1,num_Rx);
    dist=sqrt((Xrx-Xtx).^2+(Yrx-Ytx).^2); % the distance between transmitter and receivers

    RxAngle=zeros(1,num_Rx);
    RxAngle_degree=zeros(1,num_Rx);
    
    for i=1:num_Rx
        RxAngle(i)=atan((Yrx(i)-Ytx)/(Xrx(i)-Xtx)); 
        if Xrx(i)-Xtx<0 && Yrx(i)-Ytx>=0
            RxAngle(i)=RxAngle(i)+pi;
        end
        if Xrx(i)-Xtx<0 && Yrx(i)-Ytx<0
            RxAngle(i)=RxAngle(i)-pi;
        end
        RxAngle_degree(i)=RxAngle(i)*180/pi;
    end
    
    Ptx_total=(10^1.8)/1000; % transmitter's transmit power (watt)
    Ptx_beam = Ptx_total; %Ptx_beam=Ptx_total/num_beam;
    Ptx_dBm = 10*log10(Ptx_beam*1000); % transmit power (dBm) = 18 dBm
    
    degree_index_rx=zeros(1,num_Rx);
    maxBeam_index_rx=zeros(1,num_Rx);
    
    for i=1:num_Rx
        %figure; dbz(phi, g_normalized(i,:), 30 ,20);
        [~ , degree_index_rx(i)]=min(abs(phi(:)-abs(RxAngle(i))));
        [max_gain , maxBeam_index_rx(i)]=max(g(:, degree_index_rx(i)));
    end
    
    beam_frequency_table = tabulate(maxBeam_index_rx(:));
    [beam_frequency_sort, beam_index_sort]=sort(transpose(beam_frequency_table(:,2)));
    beam_chosen = transpose(beam_frequency_table(beam_index_sort(end:-1:end-(4-1)),1));
    
    for i=1:num_Rx
        flag = find(beam_chosen==maxBeam_index_rx(i));
        if isempty(flag)
            [~, beam_index_from_chosen] = max(g(beam_chosen, degree_index_rx(i)));
            maxBeam_index_rx(i)=beam_chosen(beam_index_from_chosen);
            %disp(beam_index_from_chosen);
            TXAntennaGain=10*log10(g(maxBeam_index_rx(i),degree_index_rx(i))); %dB
        else
            TXAntennaGain=10*log10(g(maxBeam_index_rx(i),degree_index_rx(i))); %dB
        end
        clear flag;
        RXAntennaGain=1;%dB
        M = Wavelength / (4 * pi * dist(i));
        Prx_dBm(i) = Ptx_dBm + TXAntennaGain + RXAntennaGain + (20*log10(M)); %Friis equation
        I_total=0; % interference
        for j=1:num_Rx
            if j == i
                continue;
            else
                flag = find(beam_chosen==maxBeam_index_rx(j));
                if isempty(flag)
                    [~, beam_index_from_chosen] = max(g(beam_chosen, degree_index_rx(j)));
                    maxBeam_index_rx(j)=beam_chosen(beam_index_from_chosen);
                    TXAntennaGain=10*log10(g(maxBeam_index_rx(j),degree_index_rx(i))); %dB
                else
                    TXAntennaGain=10*log10(g(maxBeam_index_rx(j),degree_index_rx(i))); %dB
                end
                clear flag;
                %TXAntennaGain=10*log10(g(maxBeam_index_rx(n,j),degree_index_rx(i),n)); %dB
                I_dBm = Ptx_dBm + TXAntennaGain + RXAntennaGain + (20*log10(M));
                I_total = I_total + 10^(I_dBm/10);
            end 
        end
        Prx = 10^(Prx_dBm(i)/10);
        SINR(i) = Prx / (I_total + 10^(Noise_dBm/10));
        Rate(i) = bandwidth * log2(1 + SINR(i));
        %SINR_dB(n,i)=10*log10(SINR(n,i));
    end
    Rate_avg(n)=mean(Rate);
    
end

figure;
plot(1:num_users_set, Rate_avg);
xlabel('Number of users: 20 30 40 50 60 70 80');
ylabel('Data rate (Mbps)');    

%{
figure; plot(psi/pi, sqrt(g_normalized));
legend('Rx1','Rx2','Rx3','Rx4'); 
%}

%}




