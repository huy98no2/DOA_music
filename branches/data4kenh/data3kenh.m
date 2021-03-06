%Author: Huy Vuong
%Project: Sonar 3 channel + Distance
%Lab: WICOMLAB
clear;
clc;
[y, fs]=audioread('CV_TN_80_165_200kHz_2.wav');
%frame=length(y)/50;
frame = 84000;
[Bf1 Af1]=butter(5,[5/12 1/2]);%40kHz----->50kHz
[Bf2 Af2]=butter(5,[55/96 65/96]);%55kHz----->65kHz
[Bf3 Af3]=butter(5,[75/96 85/96]);%75kHz----->85kHz
y1=filter(Bf1,Af1,y);% Loc y lay y1=45k
y2=filter(Bf2,Af2,y);%Loc y lay y2=60k
y3=filter(Bf3,Af3,y);%Loc y lay y2=80k
y1=y1/(max(y1));
y2=y2/(max(y2));
y3=y3/(max(y3));
y1 = y1(1:576000);
y2 = y2(1:576000);
y3 = y3(1:576000);
%===============
figure(1)
subplot(3,1,1);
plot(y1);
subplot(3,1,2);
plot(y2);
subplot(3,1,3);
plot(y3);
%===============
dt = 1/fs;                   % seconds per sample
StopTime = 3;                % seconds
t = (0:dt:StopTime-dt)';     % seconds
%%Sine wave:
Fc1 = 33000;                     % 45k-33k = 12k
Fc2 = 48000;                     % 60k-48k = 12k
Fc3 = 68000;                     % 80k-68k = 12k
x1 = cos(2*pi*Fc1*t);%Song sine de ha tan 1
x2 = cos(2*pi*Fc2*t);%Song sine de ha tan 2
x3 = cos(2*pi*Fc3*t);%Song sine de ha tan 3
y_nhan1 = y1.*x1;
y_nhan2 = y2.*x2;
y_nhan3 = y3.*x3;
[b,a]=butter(5,[11000,13000]/(fs/2),'bandpass');
y_12k_1 = filter(b,a,y_nhan1);  %filtered signal
y_12k_2 = filter(b,a,y_nhan2);  %filtered signal
y_12k_3 = filter(b,a,y_nhan3);  %filtered signal
figure(2)
subplot(3,1,1);
plot(y_12k_1);
title('tin hieu 12k kenh 1');
subplot(3,1,2);
plot(y_12k_2);
title('tin hieu 12k kenh 2');
subplot(3,1,3);
plot(y_12k_3);
title('tin hieu 12k kenh 3');
data_3kenh_12k = [y_12k_1 y_12k_2 y_12k_3];
figure(4)
subplot(2,1,1);
plot(data_3kenh_12k);

%==================
%Tim goc toi
%==================
M=3;
f = 12000;
v = 1500;
lamda=v/f;
gain = 0;
k=2*pi/lamda;
d=lamda/2;

yy = hilbert(data_3kenh_12k);
Nb = length(data_3kenh_12k); %so mau
UU = yy';
dt = 1/fs;
t = 0:dt:(length(data_3kenh_12k)*dt)-dt;
Ruu=UU*UU'/Nb;
[eigVector,eigValue]=eig(Ruu);
eigValueMax=max(max(eigValue));
signals = 2;
eigVectorNoise=eigVector(:,1:M-signals);
i=1; 
for theta=1:.1:180
    i=i+1; 
    A0_tmp=10^(gain/10)*exp(j*k*(0:M-1)*d*(cos(theta*pi/180))); 
    A0=A0_tmp.';
    P(i)=10*log((A0'*A0)/(A0'*eigVectorNoise*eigVectorNoise'*A0))/10;
    %P(i)=((A0'*A0)/(A0'*eigVectorNoise*eigVectorNoise'*A0))
    %[maxP,index] = max(real(P));
end

for index=1:length(P)
    if (index == 1) && (P(1) > P(2))
        disp(index/10);
    elseif (index == length(P)) && (P(length(P)) > P(length(P)-1))
        disp(index/10);
    elseif (index~= length(P))&& (index ~= 1) && (P(index) > P(index-1)) && (P(index) > P(index+1))
        disp(index/10);
    end
end 

subplot(2,1,2);
theta=1:.1:180;
plot(theta,real(P(1,1:1791)),'k','linewidth',2); 
xlabel('DOA(do)');
ylabel('Pho khong gian MUSIC(dB)'); 
hold on;
%================== 
%Plot Enveloped Signal
%==================
y_12k_1 = abs(y_12k_1);
  y1_enveloped=envelope(y_12k_1,150,'rms'); %lay mau 150 mau/lan
  y1_enveloped=y1_enveloped/(max(y1_enveloped));
  y1_enveloped = y1_enveloped(1:100000);
  figure(3)
  subplot(3,1,1);
  plot(y1_enveloped);
  xlabel('Time(s)');
  ylabel('Amplitude(V)');
  title('ENVELOPE SIGNAL CHANNEL 1');
  
  
y_12k_3 = abs(y_12k_3);
  y3_enveloped=envelope(y_12k_3,250,'rms'); %lay mau 150 mau/lan
  y3_enveloped=y3_enveloped/(max(y3_enveloped));
  y3_enveloped = y3_enveloped(1:100000);
  figure(3)
  subplot(3,1,3);
  plot(y3_enveloped);
  xlabel('Time(s)');
  ylabel('Amplitude(V)');
  title('ENVELOPE SIGNAL CHANNEL 3');
  
  y_12k_2 = abs(y_12k_2);
  y2_enveloped=envelope(y_12k_2,250,'rms'); %lay mau 150 mau/lan
  y2_enveloped=y2_enveloped/(max(y2_enveloped));
  y2_enveloped = y2_enveloped(1:100000);
  figure(3)
  subplot(3,1,2);
  plot(y2_enveloped);
  xlabel('Time(s)');
  ylabel('Amplitude(V)');
  title('ENVELOPE SIGNAL CHANNEL 2');
  
%================== 
%Xu li kenh 1
%==================
 i = 1;
 stt = 0;
 D1=[];
 khoangcach1=0;
while i<=(frame)
    switch stt
        case 0
            if(y1_enveloped(i)>= 0.4) 
                index1 = i;
                i = i+400;
                stt = 1;
            else
                i = i+1;
                stt = 0;
            end
        case 1
            if(y1_enveloped(i)>= max(y1_enveloped(i:i+16000)))
                index2 = i;
                i = i+400;
                time = index2-index1;
                khoangcach1=time*1500/2;
                disp('Distance 1(m)');
                khoangcach1=khoangcach1*1/fs;
                disp(khoangcach1);
                stt = 0;
            else
                i = i+1;
                stt = 1;
            end
        otherwise
            i = i+1;
            stt = 0;
    end
   D1=[D1 khoangcach1];
end

%================== 
%Xu li kenh 2
%==================
 ii = 1;
 stt2 = 0;
 D2=[];
 khoangcach2=0;
while ii<=(frame)
    switch stt2
        case 0
            if(y2_enveloped(ii)>= 0.4) 
                index1 = ii;
                ii = ii+400;
                stt2 = 1;
            else
                ii = ii+1;
                stt2 = 0;
            end
        case 1
            if(y2_enveloped(ii)>= max(y2_enveloped(ii:ii+1600)))
                index2 = ii;
                ii = ii+400;
                time = index2-index1;
                khoangcach2=time*1500/2;
                disp('Distance 2(m)');
                khoangcach2=khoangcach2*1/fs;
                disp(khoangcach2);
                stt2 = 0;
            else
                ii = ii+1;
                stt2 = 1;
            end
        otherwise
            ii = ii+1;
            stt2 = 0;
    end
   D2=[D2 khoangcach2];
end

%================== 
%Xu li kenh 3
%==================
 iii = 1;
 stt3 = 0;
 D3=[];
 khoangcach3=0;
while iii<=(frame)
    switch stt3
        case 0
            if(y3_enveloped(iii)>= 0.5) 
                index1 = iii;
                iii = iii+400;
                stt3 = 1;
            else
                iii = iii+1;
                stt3 = 0;
            end
        case 1
            if(y3_enveloped(iii)>= max(y3_enveloped(iii:iii+1600)))
                index2 = iii;
                iii = iii+400;
                time = index2-index1;
                khoangcach3=time*1500/2;
                disp('Distance 3(m)');
                khoangcach3=khoangcach3*1/fs;
                disp(khoangcach3);
                stt3 = 0;
            else
                iii = iii+1;
                stt3 = 1;
            end  
        otherwise
            iii = iii+1;
            stt3 = 0;
    end
   D3=[D3 khoangcach3];
end


%=============
%Ve anh
%=============

cc=[size(D1,2) , size(D2,2),size(D3,2)];
   b_max=max(cc);
   D1=[D1 zeros(1,b_max-size(D1,2))];
   D2=[D2 zeros(1,b_max-size(D2,2))];
   D3=[D3 zeros(1,b_max-size(D3,2))];
   DD1=[]; DD_1=[];
   DD2=[]; DD_2=[];
   DD3=[]; DD_3=[];
   j=1;
   while j<(length(D1)-10)
    DD1=sum(D1(j:(j+9)));
    DD1=DD1/length(D1(j:j+9));
    DD_1=[DD_1 DD1];
    
    DD2=sum(D2(j:j+9));
    DD2=DD2/length(D2(j:j+9));
    DD_2=[DD_2 DD2];
   
    DD3=sum(D3(j:j+9));
    DD3=DD3/length(D3(j:j+9));
    DD_3=[DD_3 DD3];
    j=j+10;
   end
   D=[D1;D2;D3];
D=-D;
figure
size_D=size(D);
[X,Y]=meshgrid(1:size_D(2),1:size_D(1));
ax1=subplot(1,1,1);
mesh(X,Y,D);
figure
contourf(X,Y,D);
colorbar;     
% plot3(X,Y,D);
ylabel('Legth(0.3m)');
zlabel('Depth(m)');
xlabel('Time(0.01s)')
%axis([0  100000 0 4 -10 0]);
title(' The bottom surface Cong vien thong nhat 2D Image');
colormap jet; 

DD1=0;
for i=1:size(D,1)
    for j=1:size(D,2)
    DD1=DD1+D(i,j);
    end
end 
DD1=DD1/(size(D,1)*size(D,2));
disp('Do sau trung binh');
disp(DD1);