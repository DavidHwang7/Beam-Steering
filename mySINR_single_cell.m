close all;
addpath ./ewa_function; 
d=0.5; % the number of wavelengths between antenna array
%ph0=60; % the initialized degree of main lobe direction
N=16; % the number of antenna in antenna array for an mainlobe
num_sample=400;
Xtx=500; % transimtter coordinate in x-axis (1 meter for 1 unit)
Ytx=500; % transimtter coordinate in y-axis
bandwidth = 20; % MHz, for calculating data rate


%Ptx=(10^1.8)/1000; % transmitter's transmit power (watt)
%Ptx_dBm = 10*log10(Ptx*1000); % transmit power (dBm) = 18 dBm
Wavelength=10e-3; % millimeter wavelength (meter)
Noise_dBm = -90; % Noise power (dBm) = -90 dBm
phi = (0 : num_sample) * pi / num_sample;     % equally-spaced over [0,pi]       
psi = 2 * pi * d * cos(phi);


%% Experiment 1: 固定50 users , 固定4個beam, codebook分別間隔60、45、30、15、10、5(度)
num_Rx=50;
num_beam=4;
num_time_blocks=100; % 每個user set都run 100~1000次
%codebook_interval = [90 60 45 30 15 10 5];
codebook_interval = [60 45 30 15 10 5];
[~, num_codebook] = size(codebook_interval);
Xrx = randi([0,1000],1,50);
Yrx = randi([0,1000],1,50);
%rand_index = randperm(50);
%draw_rand_index = rand_index(1:25);
%Xrx(draw_rand_index) = -1 * Xrx1(draw_rand_index);
%Yrx = randi([100,500],1,50);


Ptx_total=(10^1.8)/1000; % transmitter's transmit power (watt)
Ptx_beam = Ptx_total;
%Ptx_beam=Ptx_total/num_beam;
Ptx_dBm = 10*log10(Ptx_beam*1000); % transmit power (dBm) = 18 dBm


dist=sqrt((Xrx-Xtx).^2+(Yrx-Ytx).^2); % the distance between transmitter and receivers
side=Xrx-Xtx;

RxAngle=zeros(1,num_Rx);
RxAngle_degree=zeros(1,num_Rx);

for i=1:num_Rx
    RxAngle(i)=atan((Yrx(i)-Ytx)/(Xrx(i)-Xtx))+pi/2;
end

Prx_dBm = zeros(num_codebook, num_Rx, num_time_blocks);
SINR1 = zeros(num_codebook, num_Rx, num_time_blocks);

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
    
    degree_index_rx=zeros(1,num_Rx);
    maxBeam_index_rx=zeros(1,num_Rx);
    
    for i=1:num_Rx
        %figure; dbz(phi, g_normalized(i,:), 30 ,20);
        [~ , degree_index_rx(i)]=min(abs(phi(:)-abs(RxAngle(i))));
        [max_gain , maxBeam_index_rx(i)]=max(g_current(:, degree_index_rx(i)));
        %max = max(g(:,degree_index_rx(i),n));
        %fprintf('Max gain beam:%d Rx:%d',n,i); 
        %disp(max_gain);
        %Prx_dBm(i)=0; % receiver power from mainlobe
    end
    %{
    if current_codebook_size > num_beam
        beam_frequency_table = tabulate(maxBeam_index_rx(:));
        [beam_frequency_sort, beam_index_sort]=sort(transpose(beam_frequency_table(:,2)));
        beam_chosen = transpose(beam_frequency_table(beam_index_sort(end:-1:end-(4-1)),1));
    else
        beam_chosen = [1 2 3 4];
    end
    %}
    user_index=1;
    for t=1:num_time_blocks
        user_chosen=zeros(1,num_beam);
        beam_chosen=zeros(1,num_beam);
        
    %%%%% 選User和Beam %%%%%    
        for i=1:num_beam
            user_chosen(i)=user_index;
            beam_temp = maxBeam_index_rx(user_index);
            if side(user_index)<0 && beam_temp~=1 && beam_temp~=current_codebook_size
                beam_temp = -beam_temp;
            end

            jj=2;
            g_temp=g_current(:, degree_index_rx(user_index));
            [~, beam_sort]=sort(g_temp,'descend');
            while isempty(find(beam_chosen==beam_temp,1))==false
                beam_temp = beam_sort(jj);
                if side(user_index)<0 && beam_temp~=1 && beam_temp~=current_codebook_size
                    beam_temp = -beam_temp;
                end
                jj =jj+1;
            end
            beam_chosen(i)=beam_temp;
            
            user_index = user_index+1;
            if user_index > num_Rx
                user_index=1;
            end
        end
        
    %%%%% 算interference %%%%%    
        for i=1:num_beam
            if beam_chosen(i)<0
                beam_current=-beam_chosen(i);
            else
                beam_current=beam_chosen(i);
            end 
            TXAntennaGain=10*log10(g_current(beam_current,degree_index_rx(user_chosen(i)))); %dB
            RXAntennaGain=1;%dB
            M = Wavelength / (4 * pi * dist(user_chosen(i)));
            Prx_dBm(n,user_chosen(i),t) = Ptx_dBm + TXAntennaGain + RXAntennaGain + (20*log10(M)); %Friis equation
            I_total=0; % interference
            j=1;
            for j=1:num_beam
                if j ~= i || side(user_chosen(i))*side(user_chosen(j))>=0 %%判斷是否為同一個user＆判斷是否在同一側，同側的話side相乘會>=0
                    if beam_chosen(j)<0
                        beam_current=-beam_chosen(j);
                    else
                        beam_current=beam_chosen(j);
                    end
                    TXAntennaGain=10*log10(g_current(beam_current,degree_index_rx(user_chosen(i)))); %dB
                    I_dBm = Ptx_dBm + TXAntennaGain + RXAntennaGain + (20*log10(M));
                    I_total = I_total + 10^(I_dBm/10);
                end
            end
            Prx = 10^(Prx_dBm(n,user_chosen(i),t)/10);
            SINR1(n,user_chosen(i),t) = Prx / (I_total + 10^(Noise_dBm/10));
            %Rate1(n,user_chosen(i),t) = bandwidth * log2(1 + SINR1(n,user_chosen(i),t));
            
        end
    end
    %{
    for i=1:num_Rx
        flag = find(beam_chosen==maxBeam_index_rx(i));
        if isempty(flag)
            [~, beam_index_from_chosen] = max(g_current(beam_chosen, degree_index_rx(i)));
            maxBeam_index_rx(i)=beam_chosen(beam_index_from_chosen);
            %disp(beam_index_from_chosen);
            TXAntennaGain=10*log10(g_current(maxBeam_index_rx(i),degree_index_rx(i))); %dB
        else
            TXAntennaGain=10*log10(g_current(maxBeam_index_rx(i),degree_index_rx(i))); %dB
        end
        clear flag;
        RXAntennaGain=1;%dB
        M = Wavelength / (4 * pi * dist1(i));
        Prx_dBm(n,i) = Ptx_dBm + TXAntennaGain + RXAntennaGain + (20*log10(M)); %Friis equation
        I_total=0; % interference
        for j=1:num_Rx
            if j == i
                continue;
            else
                flag = find(beam_chosen==maxBeam_index_rx(j));
                if isempty(flag)
                    [~, beam_index_from_chosen] = max(g_current(beam_chosen, degree_index_rx(j)));
                    maxBeam_index_rx(j)=beam_chosen(beam_index_from_chosen);
                    TXAntennaGain=10*log10(g_current(maxBeam_index_rx(j),degree_index_rx(i))); %dB
                else
                    TXAntennaGain=10*log10(g_current(maxBeam_index_rx(j),degree_index_rx(i))); %dB
                end
                %TXAntennaGain=10*log10(g(maxBeam_index_rx(n,j),degree_index_rx(i),n)); %dB
                I_dBm = Ptx_dBm + TXAntennaGain + RXAntennaGain + (20*log10(M));
                I_total = I_total + 10^(I_dBm/10);
            end 
        end
        Prx = 10^(Prx_dBm(n,i)/10);
        SINR1(n,i) = Prx / (I_total + 10^(Noise_dBm/10));
        Rate1(n,i) = bandwidth * log2(1 + SINR1(n,i));
        %SINR_dB(n,i)=10*log10(SINR(n,i));
        %fprintf('Prx%d(dBm) = %.4f\n',i,Prx_dBm(i));
        %fprintf('Rx%d''s SINR(dB) = %.4f\n\n',i,SINR_dB);
        %disp(Prx_dBm);
        %fprintf('Tx Antenna Gain toward Rx%d = %.4f\n\n',i,TXAntennaGain);
        %disp(TXAntennaGain);
    end
    %SINR_avg_dB(n) = 10*log10(mean(SINR(n,:)));
    Rate_avg1(n)=mean(Rate1(n,:));
    %}
    
    clear a_current A_current g_current;
end

Rate = bandwidth * log2(1 + SINR1);
Total_Mb_per_user = sum(Rate,3);
Rate_average_codebook_single = sum(Total_Mb_per_user,2)/num_time_blocks;
figure;
%plot(1:num_codebook, Rate_average_codebook);
plot(codebook_interval, Rate_average_codebook_single, '-*b');
set(gca,'XTick',[5 10 15 30 45 60]);
%axis([0, 60]);  %确定x?与y?框?大小
xlabel('codebook interval (degree)');
ylabel('Data rate (Mbps)');

%{
x=1:1:5;%x?上的?据，第一?值代表?据?始，第二?值代表?隔，第三?值代表?止
plot(x,a,'-*b',x,b,'-or'); %?性，?色，??
axis([0,6,0,700]);  %确定x?与y?框?大小
set(gca,'XTick',0:1:6); %x?范?1-6，?隔1
set(gca,'YTick',0:100:700); %y?范?0-700，?隔100
%}
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




