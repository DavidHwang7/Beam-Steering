function [g_saved] = Gain_simulation(N,B)
% N is the number of antenna in antenna array for an mainlobe
% B:(30,15,10,5)每隔幾度可以打一個BEAM
addpath ./ewa_function; 
d=0.5; % the number of wavelengths between antenna array
num_sample = 36; % 180/5 每隔5度取一個sample(角度)
num_beam = (180/B)+1;
beam_direction = 0:B:180; %以45度為中心，Beam可打出的角度
phi = (0 : num_sample) * pi / num_sample;     % equally-spaced over [0,pi]       
psi = 2 * pi * d * cos(phi);

a=zeros(1,N); 
A=zeros(1,num_sample+1);
g=zeros(1,num_sample+1);
g_saved=zeros(num_beam,37); % save -45 ~ 45

for i=1:num_beam
    ph0 = beam_direction(i); % the initialized degree of main lobe direction
    a(1,:) = uniform(d, ph0, N); % antenna array weight factor
    A(1,:) = dtft(a(1,:), -psi);  % array factor, note dsft(a,psi)=dtft(a,-psi)
    g(1,:)= abs(A(1,:)).^2;  % power gain
    g_saved(i,:)=g(1:37); % g(2,10:28) = 45 ~ 135
end



%filepath = '.\beam_data4';
%save(filepath,'g_saved_v2');










