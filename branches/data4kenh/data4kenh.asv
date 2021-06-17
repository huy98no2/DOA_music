%Author: Huy Vuong
%Project: Sonar 4 channel + Distance +MUSIC
%Lab: WICOMLAB
clear;
clc;
[y, fs]=audioread('congvietTN_1704_14h30.wav');
%frame=length(y)/50;
frame = 84000;
 [Bf1 Af1]=butter(5,[5/48, 5/24]);%10kHz ---->20kHz
 [Bf2 Af2]=butter(5,[5/24, 5/16]);%20kHz----->30kHz
 [Bf3 Af3]=butter(5,[5/16, 5/12]);%30kHz----->40kHz
 [Bf4 Af4]=butter(5,[5/12, 1/2]); %40kHz----->50kHz
y1=filter(Bf1,Af1,y);%Loc y lay y1=15k
y2=filter(Bf2,Af2,y);%Loc y lay y2=25k
y3=filter(Bf3,Af3,y);%Loc y lay y2=35k
y4=filter(Bf4,Af4,y);%Loc y lay y2=45k
%Chuan hoa tin hieu cac kenh
y1=y1/(max(y1));
y2=y2/(max(y2));
y3=y3/(max(y3));
y4=y4/(max(y4));
y1 = y1(1:576000);
y2 = y2(1:576000);
y3 = y3(1:576000);
y4 = y4(1:576000);
%===============
figure(1)
subplot(4,1,1);
plot(y1);
subplot(4,1,2);
plot(y2);
subplot(4,1,3);
plot(y3);
subplot(4,1,4);
plot(y4);
% subplot(4,2,1);
% plot(psd(spectrum.periodogram,y1,'Fs',fs,'NFFT',length(y1)));
% subplot(4,2,2);
% plot(psd(spectrum.periodogram,y2,'Fs',fs,'NFFT',length(y2)));
% subplot(4,2,3);
% plot(psd(spectrum.periodogram,y3,'Fs',fs,'NFFT',length(y3)));
% subplot(4,2,4);
% plot(psd(spectrum.periodogram,y4,'Fs',fs,'NFFT',length(y4)));
%===============
dt = 1/fs;                   % seconds per sample
StopTime = 3;                % seconds
t = (0:dt:StopTime-dt)';     % seconds
%%Sine wave:
Fc1 = 3000;                      % 15k- 3k = 12k
Fc2 = 13000;                     % 25k-13k = 12k
Fc3 = 23000;                     % 35k-23k = 12k
Fc4 = 33000;                     % 45k-33k = 12k
x1 = cos(2*pi*Fc1*t);%Song sine de ha tan 1
x2 = cos(2*pi*Fc2*t);%Song sine de ha tan 2
x3 = cos(2*pi*Fc3*t);%Song sine de ha tan 3
x4 = cos(2*pi*Fc4*t);%Song sine de ha tan 4
y_nhan1 = y1.*x1;
y_nhan2 = y2.*x2;
y_nhan3 = y3.*x3;
y_nhan4 = y4.*x4;
[b,a]=butter(5,[11000,13000]/(fs/2),'bandpass');
y_12k_1 = filter(b,a,y_nhan1);  %filtered signal
y_12k_2 = filter(b,a,y_nhan2);  %filtered signal
y_12k_3 = filter(b,a,y_nhan3);  %filtered signal
y_12k_4 = filter(b,a,y_nhan4);  %filtered signal
figure(2);
subplot(4,1,1);
plot(y_12k_1);
title('tin hieu 12k kenh 1');
subplot(4,1,2);
plot(y_12k_2);
title('tin hieu 12k kenh 2');
subplot(4,1,3);
plot(y_12k_3);
title('tin hieu 12k kenh 3');
subplot(4,1,4);
plot(y_12k_4);
title('tin hieu 12k kenh 4');

% subplot(4,2,1);
% plot(psd(spectrum.periodogram,y_12k_1,'Fs',fs,'NFFT',length(y_12k_1)));
% subplot(4,2,2);
% plot(psd(spectrum.periodogram,y_12k_2,'Fs',fs,'NFFT',length(y_12k_2)));
% subplot(4,2,3);
% plot(psd(spectrum.periodogram,y_12k_3,'Fs',fs,'NFFT',length(y_12k_3)));
% subplot(4,2,4);
% plot(psd(spectrum.periodogram,y_12k_4,'Fs',fs,'NFFT',length(y_12k_4)));

data_4kenh_12k = [y_12k_1 y_12k_2 y_12k_3 y_12k_4];
figure(4);
subplot(2,1,1);
plot(data_4kenh_12k);

%==================
%Tim goc toi
%==================
M=4;
f = 12000;
v = 1500;
lamda=v/f;
gain = 0;
k=2*pi/lamda;
d=lamda/2;
signals = 3;
%==============
%Goc tai vi tri 1
%==============
yy1 = hilbert(data_4kenh_12k(1000:14400,:));
Nb1 = length(data_4kenh_12k(1000:14400,:)); %so mau
UU1 = yy1';
dt = 1/fs;
t = 0:dt:(length(data_4kenh_12k(1000:14400,:))*dt)-dt;
Ruu1=UU1*UU1'/Nb1;
[eigVector1,eigValue1]=eig(Ruu1);
eigValueMax1=max(max(eigValue1));
eigVectorNoise1=eigVector1(:,1:M-signals);
i=1; 
for theta1=1:.1:180
    i=i+1; 
    A0_tmp1=10^(gain/10)*exp(j*k*(0:M-1)*d*(cos(theta1*pi/180))); 
    A01=A0_tmp1.';
    P1(i)=10*log((A01'*A01)/(A01'*eigVectorNoise1*eigVectorNoise1'*A01))/10;
    %P(i)=((A0'*A0)/(A0'*eigVectorNoise*eigVectorNoise'*A0))
    %[maxP,index] = max(real(P));
end

for index1=1:length(P1)
%     if (index == 1) && (P(1) > P(2))
%         disp('goc toi vi tri 1');
%         disp(index/10-90);
%     elseif (index == length(P)) && (P(length(P)) > P(length(P)-1))
%         disp('goc toi');
%         disp(index/10-90);
    if (index1~= length(P1))&& (index1 ~= 1) && (P1(index1) > P1(index1-1)) && (P1(index1) > P1(index1+1))
        %disp('goc toi vi tri 1');
        %disp(index1/10-90);
    end
end 

%==============
%Goc tai vi tri 2
%==============
yy2 = hilbert(data_4kenh_12k(23000:33500,:));
Nb2 = length(data_4kenh_12k(23000:33500,:)); %so mau
UU2 = yy2';
%dt = 1/fs;
t2 = 0:dt:(length(data_4kenh_12k(23000:33500,:))*dt)-dt;
Ruu2=UU2*UU2'/Nb2;
[eigVector2,eigValue2]=eig(Ruu2);
eigValueMax2=max(max(eigValue2));
eigVectorNoise2=eigVector2(:,1:M-signals);
i2=1; 
for theta2=1:.1:180
    i2=i2+1; 
    A0_tmp2=10^(gain/10)*exp(j*k*(0:M-1)*d*(cos(theta2*pi/180))); 
    A02=A0_tmp2.';
    P2(i2)=10*log((A02'*A02)/(A02'*eigVectorNoise2*eigVectorNoise2'*A02))/10;
    %P(i)=((A0'*A0)/(A0'*eigVectorNoise*eigVectorNoise'*A0))
    %[maxP,index] = max(real(P));
end

for index2=1:length(P2)
%     if (index == 1) && (P(1) > P(2))
%         disp('goc toi vi tri 1');
%         disp(index/10-90);
%     elseif (index == length(P)) && (P(length(P)) > P(length(P)-1))
%         disp('goc toi');
%         disp(index/10-90);
    if (index2~= length(P2))&& (index2 ~= 1) && (P2(index2) > P2(index2-1)) && (P2(index2) > P2(index2+1))
        %disp('goc toi vi tri 2');
        %disp(index2/10-90);
    end
end 
%==============
%Goc tai vi tri 3
%==============
yy3 = hilbert(data_4kenh_12k(481800:493800,:));
Nb3 = length(data_4kenh_12k(481800:493800,:)); %so mau
UU3 = yy3';
%dt = 1/fs;
t3 = 0:dt:(length(data_4kenh_12k(481800:493800,:))*dt)-dt;
Ruu3=UU3*UU3'/Nb3;
[eigVector3,eigValue3]=eig(Ruu3);
eigValueMax3=max(max(eigValue3));
eigVectorNoise3=eigVector3(:,1:M-signals);
i3=1; 
for theta3=1:.1:180
    i3=i3+1; 
    A0_tmp3=10^(gain/10)*exp(j*k*(0:M-1)*d*(cos(theta3*pi/180))); 
    A03=A0_tmp3.';
    P3(i3)=10*log((A03'*A03)/(A03'*eigVectorNoise3*eigVectorNoise3'*A03))/10;
    %P(i)=((A0'*A0)/(A0'*eigVectorNoise*eigVectorNoise'*A0))
    %[maxP,index] = max(real(P));
end

for index3=1:length(P3)
%     if (index == 1) && (P(1) > P(2))
%         disp('goc toi vi tri 1');
%         disp(index/10-90);
%     elseif (index == length(P)) && (P(length(P)) > P(length(P)-1))
%         disp('goc toi');
%         disp(index/10-90);
    if (index3~= length(P3))&& (index3 ~= 1) && (P3(index3) > P3(index3-1)) && (P3(index3) > P3(index3+1))
        %disp('goc toi vi tri 3');
        %disp(index3/10-90);
    end
end 
%==============
%Goc tai vi tri 3
%==============
yy4 = hilbert(data_4kenh_12k(500000:513000,:));
Nb4 = length(data_4kenh_12k(500000:513000,:)); %so mau
UU4 = yy4';
%dt = 1/fs;
t4 = 0:dt:(length(data_4kenh_12k(500000:513000,:))*dt)-dt;
Ruu4=UU4*UU4'/Nb4;
[eigVector4,eigValue4]=eig(Ruu4);
eigValueMax4=max(max(eigValue4));
eigVectorNoise4=eigVector4(:,1:M-signals);
i4=1; 
for theta4=1:.1:180
    i4=i4+1; 
    A0_tmp4=10^(gain/10)*exp(j*k*(0:M-1)*d*(cos(theta4*pi/180))); 
    A04=A0_tmp4.';
    P4(i4)=10*log((A04'*A04)/(A04'*eigVectorNoise4*eigVectorNoise4'*A04))/10;
    %P(i)=((A0'*A0)/(A0'*eigVectorNoise*eigVectorNoise'*A0))
    %[maxP,index] = max(real(P));
end

for index4=1:length(P4)
%     if (index == 1) && (P(1) > P(2))
%         disp('goc toi vi tri 1');
%         disp(index/10-90);
%     elseif (index == length(P)) && (P(length(P)) > P(length(P)-1))
%         disp('goc toi');
%         disp(index/10-90);
    if (index4~= length(P4))&& (index4 ~= 1) && (P4(index4) > P4(index4-1)) && (P4(index4) > P4(index4+1))
        %disp('goc toi vi tri 4');
        %disp(index4/10-90);
    end
end
%=========
%Plot pho khong gian DOA
%=========
% subplot(2,1,2);
% theta=-89:.1:90;
% plot(theta,real(P(1,1:1791)),'k','linewidth',2); 
% xlabel('DOA(do)');
% ylabel('Pho khong gian MUSIC(dB)'); 
% hold on;
%================== 
%Plot Enveloped Signal
%==================
y_12k_1 = abs(y_12k_1);
  y1_enveloped=envelope(y_12k_1,150,'rms'); %lay mau 150 mau/lan
  y1_enveloped=y1_enveloped/(max(y1_enveloped));
  y1_enveloped = y1_enveloped(1:100000);
  figure(3)
  subplot(4,1,1);
  plot(y1_enveloped);
  xlabel('Time(s)');
  ylabel('Amplitude(V)');
  title('ENVELOPE SIGNAL CHANNEL 1');
  
y_12k_2 = abs(y_12k_2);
  y2_enveloped=envelope(y_12k_2,250,'rms'); %lay mau 150 mau/lan
  y2_enveloped=y2_enveloped/(max(y2_enveloped));
  y2_enveloped = y2_enveloped(1:100000);
  figure(3)
  subplot(4,1,2);
  plot(y2_enveloped);
  xlabel('Time(s)');
  ylabel('Amplitude(V)');
  title('ENVELOPE SIGNAL CHANNEL 2');
  
y_12k_3 = abs(y_12k_3);
  y3_enveloped=envelope(y_12k_3,250,'rms'); %lay mau 150 mau/lan
  y3_enveloped=y3_enveloped/(max(y3_enveloped));
  y3_enveloped = y3_enveloped(1:100000);
  figure(3)
  subplot(4,1,3);
  plot(y3_enveloped);
  xlabel('Time(s)');
  ylabel('Amplitude(V)');
  title('ENVELOPE SIGNAL CHANNEL 3');
  
  y_12k_4 = abs(y_12k_4);
  y4_enveloped=envelope(y_12k_4,150,'rms'); %lay mau 150 mau/lan
  y4_enveloped=y4_enveloped/(max(y4_enveloped));
  y4_enveloped = y4_enveloped(1:100000);
  figure(3)
  subplot(4,1,4);
  plot(y4_enveloped);
  xlabel('Time(s)');
  ylabel('Amplitude(V)');
  title('ENVELOPE SIGNAL CHANNEL 4');
  
  
  
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
            if(y1_enveloped(i)>= max(y1_enveloped(i:i+1600)))
                index2 = i;
                i = i+400;
                time1 = index2-index1;
                khoangcach1=time1*1500/2;
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
%  ii = 1;
%  stt2 = 0;
%  D2=[];
%  khoangcach2=0;
% while ii<=(frame)
%     switch stt2
%         case 0
%             if(y2_enveloped(ii)>= 0.5) 
%                 index1 = ii;
%                 ii = ii+400;
%                 stt2 = 1;
%             else
%                 ii = ii+1;
%                 stt2 = 0;
%             end
%         case 1
%             if(y2_enveloped(ii)>= max(y2_enveloped(ii:ii+1600)))
%                 index2 = ii;
%                 ii = ii+400;
%                 time2 = index2-index1;
%                 khoangcach2=time2*1500/2;
%                 disp('Distance 2(m)');
%                 khoangcach2=khoangcach2*1/fs;
%                 disp(khoangcach2);
%                 stt2 = 0;
%             else
%                 ii = ii+1;
%                 stt2 = 1;
%             end
%         otherwise
%             ii = ii+1;
%             stt2 = 0;
%     end
%    D2=[D2 khoangcach2];
% end

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
                time3 = index2-index1;
                khoangcach3=time3*1500/2;
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

%================== 
%Xu li kenh 4
%==================
 iiii = 1;
 stt4 = 0;
 D4=[];
 khoangcach4=0;
while iiii<=(frame)
    switch stt4
        case 0
            if(y4_enveloped(iiii)>= 0.5) 
                index1 = iiii;
                iiii = iiii+400;
                stt4 = 1;
            else
                iiii = iiii+1;
                stt4 = 0;
            end
        case 1
            if(y3_enveloped(iiii)>= max(y3_enveloped(iiii:iiii+1600)))
                index2 = iiii;
                iiii = iiii+400;
                time4 = index2-index1;
                khoangcach4=time4*1500/2;
                disp('Distance 4(m)');
                khoangcach4=khoangcach4*1/fs;
                disp(khoangcach4);
                stt4 = 0;
            else
                iiii = iiii+1;
                stt4 = 1;
            end  
        otherwise
            iiii = iiii+1;
            stt4 = 0;
    end
   D4=[D4 khoangcach4];
end

%=============
%Ve anh
%=============

% cc=[size(D1,2) , size(D2,2),size(D3,2),size(D4,2)];
%    b_max=max(cc);
%    D1=[D1 zeros(1,b_max-size(D1,2))];
%    D2=[D2 zeros(1,b_max-size(D2,2))];
%    D3=[D3 zeros(1,b_max-size(D3,2))];
%    D4=[D4 zeros(1,b_max-size(D4,2))];
%    DD1=[]; DD_1=[];
%    DD2=[]; DD_2=[];
%    DD3=[]; DD_3=[];
%    DD4=[]; DD_4=[];
%    j=1;
%    while j<(length(D1)-10)
%     DD1=sum(D1(j:(j+9)));
%     DD1=DD1/length(D1(j:j+9));
%     DD_1=[DD_1 DD1];
%     
%     DD2=sum(D2(j:j+9));
%     DD2=DD2/length(D2(j:j+9));
%     DD_2=[DD_2 DD2];
%    
%     DD3=sum(D3(j:j+9));
%     DD3=DD3/length(D3(j:j+9));
%     DD_3=[DD_3 DD3];
%     
%     DD4=sum(D4(j:j+9));
%     DD4=DD4/length(D4(j:j+9));
%     DD_4=[DD_4 DD4];
%     j=j+10;
%    end
%    D=[D1;D2;D3;D4];
% D=-D;
% figure
% size_D=size(D);
% [X,Y]=meshgrid(1:size_D(2),1:size_D(1)); %Tao luoi toa do 2D
% ax1=subplot(1,1,1);
% mesh(X,Y,D);%Ve anh 3D
% figure
% contourf(X,Y,D); % ve anh 2-D
% colorbar;     
% % plot3(X,Y,D);
% ylabel('Legth');
% zlabel('Depth(m)');
% xlabel('Time')
% %axis([0  100000 0 4 -10 0]);
% title(' The bottom surface Cong vien thong nhat 2D Image');
% colormap jet; %theme mau
% 
% DD1=0;
% for i=1:size(D,1)
%     for j=1:size(D,2)
%     DD1=DD1+D(i,j);
%     end
% end 
% DD1=DD1/(size(D,1)*size(D,2));
% disp('Do sau trung binh');
% disp(DD1);