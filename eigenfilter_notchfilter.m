%% Eigen filter method and its application: Notch filter
%
% Srilakshmi Alla
%%  Clearing memory
clc;
clear all;
close all;

%%
% N=50; % Given Order of the filter

in=input('in(try to keep it between 5 to 8):'); 

N=in*10; % if "in" goes from 5 to 8 we may be get good response
% When order is increased the error decreases and the response we
% generates equals desired response

L=16*N;   % Length 
fs=24000; % Sampling Frequency
ws=2*pi*fs; 
Ts=1/fs;

wn=(0:fs/(2*L):12000)*(2*pi); % Total Frequency Range

% wn1=(0:fs/(2*L):4000)*(2*pi);     %  Passband
% wn2=(5000:fs/(2*L):7000)*(2*pi);  %  Stopband
% wn3=(7500:fs/(2*L):8500)*(2*pi);  %  Passband
% wn4=(9000:fs/(2*L):12000)*(2*pi); %  Stopband


% Another method for finding the regions

sz=length(wn); % Length of wn

z=fs/(2*sz);

n1=0:4000/z;
n2=round(5000/z:7000/z);
n3=round(7500/z:8500/z);
n4=round(9000/z:12000/z);
 

wn1=wn(n1+1);  % Passband
wn2=wn(n2);    % Stopband
wn3=wn(n3);    % Passband
wn4=wn(n4);    % Stopband

% Both of these methods work

%% Desired Frequency Response

D1=(wn1>=0);
D2=zeros(1,length(wn2));
D3=(wn3/(1000*2*pi))-7.5;
D4=zeros(1,length(wn4));

s=2*pi;
figure(1);
plot(wn1/s,D1,wn2/s,D2,wn3/s,D3,wn4/s,D4,'linewidth',2);
grid on;
xlabel('Frequency f(Hz)');
ylabel('D');
title('Desired Frequency Response');

%% Step1: Q matrix

M=N/2;

% Initialize everything to 0
Q=zeros(M+1,M+1);
Q1=zeros(M+1,M+1);
Q2=zeros(M+1,M+1);
Q3=zeros(M+1,M+1);
Q4=zeros(M+1,M+1);

% Calculation of Q matrix

for m=0:M
    for n=0:M
        
        Q1(m+1,n+1)=1/ws*trapz(wn1,cos(m*wn1*Ts).*cos(n*wn1*Ts));
        Q2(m+1,n+1)=1/ws*trapz(wn2,cos(m*wn2*Ts).*cos(n*wn2*Ts));
        Q3(m+1,n+1)=1/ws*trapz(wn3,cos(m*wn3*Ts).*cos(n*wn3*Ts));
        Q4(m+1,n+1)=1/ws*trapz(wn4,cos(m*wn4*Ts).*cos(n*wn4*Ts));
        
        Q=Q1+Q2+Q3+Q4;
        
    end
end



%% Step1 : p matrix

% Initialize everything to 0
p=zeros(M+1,1);
p1=zeros(M+1,1);
p2=zeros(M+1,1);
p3=zeros(M+1,1);
p4=zeros(M+1,1);


% Calculation of p matrix
for m=0:M
    
    p1(m+1,:)=1/ws*trapz(wn1,cos(m*wn1*Ts).*D1);
    p2(m+1,:)=1/ws*trapz(wn2,cos(m*wn2*Ts).*D2);
    p3(m+1,:)=1/ws*trapz(wn3,cos(m*wn3*Ts).*D3);
    p4(m+1,:)=1/ws*trapz(wn4,cos(m*wn4*Ts).*D4);
    
    
    p=p1+p2+p3+p4;
    
    
end


%% Step 1:d matrix

% Initialize everything to 0

d=zeros(1,1);
d1=zeros(1,1);
d2=zeros(1,1);
d3=zeros(1,1);
d4=zeros(1,1);

% Calculation of d matrix

d1=1/ws*trapz(wn1,D1.*D1);
d2=1/ws*trapz(wn2,D2.*D2);
d3=1/ws*trapz(wn3,D3.*D3);
d4=1/ws*trapz(wn4,D4.*D4);


d=d1+d2+d3+d4;

%% Step 2:Qt

% Calculation of Qt matrix

Qt=[Q p;p' d];

%% Step3: Calculation of minimum eigenvector a0_hat of Qt

Qtt=real(Qt);
[V,Di]=eig(Qt); % eigen values and eigen vectors corresponding to those eigen values
[Dis,Ix]=sort((diag(Di)),'ascend'); % Sorting the diagonal elements to pick up the minimum eigen value

V1=V(:,Ix);
a0_hat=V1(:,1); % eigen vector corresponding to minimum eigen value
%% Verification whether Qt is positive definite or not

if Dis>0
    disp('Qt is positive definite');
else
    disp('Qt is not postive definite');
end

%% Step4: Normalization of a0_hat to obtain a0

V3=a0_hat/-(a0_hat(length(a0_hat),1)); % Normalizing the eigen vector by dividing it with the negative of last element 

a0=V3(1:length(V3)-1); % Removing the last element to get desired a0

%% Step5: Finding Impulse Response and frequency response

% Finding h(n) from a0


for k0=2:M+1
    
  hn0(M+1)=a0(1,1);
  hn0(M-k0+2)=0.5*a0(k0,1);
    
end

% Finding remaining impulse response using symmetry property

for n11=1:N+1
    hn0(N-n11+2)=hn0(n11);
end

% Impulse Response

figure(2);
stem(hn0/max(hn0));
title('Impulse Response');

% Frequency Response

WW=0:ws/2/L:ws/2;
[VVx,w0]=freqz(hn0,1,WW/2/pi,fs);

figure(3);

plot(w0,abs(VVx),'k','linewidth',2);
grid on;
hold on;
plot(wn1/s,D1,wn2/s,D2,wn3/s,D3,wn4/s,D4,'linewidth',2);
hold off;
title('Frequency response of designed filter and desired frequency response')
%% Application of Eigenfilter:Notch Filter

%% Step1: Finding Qt matrix

% Qt is taken from part A
%% Step2: Finding E_hat

% Initialize E to 0
E=zeros(2,M+1);

% Notches
w0=4500*2*pi;
w1=8000*2*pi;

% Calculation of E

n=0:M;
E(1,:)=(1/ws)*cos(w0*n*Ts);
E(2,:)=(1/ws)*cos(w1*n*Ts);

% Append 0 to E to make it E_hat

E_hat(1,:)=[E(1,:) 0];
E_hat(2,:)=[E(2,:) 0];

%% Step3: Finding orthonormal basis of the null space of E_hat

B=null(E_hat);

%% Step4:Calculation of minimum eigenvector

W1=B'*Qt*B;

[Vy,Diy]=eig(W1);
[Disy,Ixy]=sort((diag(Diy)),'ascend');
V11=Vy(:,Ixy);

W=V11(:,1); % minimum eigen vector

%% Step5:Finding a0_hat(a0y)

a0y=B*W;

%% Step6:Normalize a0_hat to obtain a0

V33=a0y/-a0y(length(a0y),1);

a01=V33(1:length(V33)-1);

%% Step7: Finding Impulse Response from a0

% Impulse Response

for k=2:M+1
          
  hn(M+1)=a01(1,1);
  hn(M-k+2)=0.5*a01(k,1);
    
end

for n11=1:N+1
    hn(N-n11+2)=hn(n11);
end


figure(4);
stem(hn);
title('Impulse Response');


WW=0:ws/2/L:ws/2;



figure(25);

VV12=freqz(hn,1,WW/2/pi,fs);
plot(WW/2/pi,abs(VV12),'g','linewidth',2);
grid on;
hold on;
plot(wn1/s,D1,wn2/s,D2,wn3/s,D3,wn4/s,D4,'linewidth',3);
plot([4500 4500],[0 1],'-.k*','linewidth',2);
plot([8000 8000],[0 1],'-.k*','linewidth',2);
hold off;
title('Frequency Response of Designed Notch filter and desired frequency response with notches')
%% Summary
%  Always Qt should be positive definite.
%  Got a better response if order is 70 or 80
