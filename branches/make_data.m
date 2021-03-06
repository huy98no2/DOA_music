% Make data
%huyvq0

mode = 0; % 0: nhieu = 0 ; 1: nhieu khac 0
[U,fs] = audioread('D:/THUAT_TOAN_MUSIC/file_wav/28-11/test_12k_2811_3.wav');
U1 = U(:,1);
U2 = U(:,2);
%random = -6e-4 + (9e-4)*rand(76704,1);
switch mode
    case 1 
    
        nhieu_1_part1 = 0.05*U1(30000:33840); 
        nhieu_2_part1 = 0.05*U2(30000:33840);
        nhieu_1_part2 = 0.05*U1(30000:102864);
        nhieu_2_part2 = 0.05*U2(30000:102864);
    case 0 
        nhieu_1_part1 = zeros(3840, 1); 
        nhieu_2_part1 = zeros(3840, 1);
        nhieu_1_part2 = zeros(72864,1);
        nhieu_2_part2 = zeros(72864,1);
end
reflex1 = 0.3*U1(30000:30095)-0.001;
reflex2 = 0.3*U2(30000:30095)-0.001;
data1 = [ zeros(10000,1); U1(30000:30095); nhieu_1_part1;reflex1;nhieu_1_part2; U1(30096:30191);nhieu_1_part1;reflex1;nhieu_1_part2;U1(30286:30381);
          nhieu_1_part1;reflex1;nhieu_1_part2;U1(30000:30095); nhieu_1_part1;reflex1;nhieu_1_part2; U1(30096:30191);nhieu_1_part1;reflex1;nhieu_1_part2;
          U1(30286:30381);nhieu_1_part1;reflex1;nhieu_1_part2 ];
data2 = [ zeros(10000,1); U2(30000:30095); nhieu_2_part1;reflex2;nhieu_2_part2; U2(30096:30191);nhieu_2_part1;reflex2;nhieu_2_part2;U2(30286:30381);
          nhieu_2_part1;reflex2;nhieu_2_part2;U2(30000:30095); nhieu_2_part1;reflex2;nhieu_2_part2; U2(30096:30191);nhieu_2_part1;reflex2;nhieu_2_part2;
          U2(30286:30381);nhieu_2_part1;reflex2;nhieu_2_part2 ];
data = [data1 data2];
subplot(3,1,1);
plot(data);
audiowrite('data_2kenh_radar_12k.wav',data,192000);
%Loc 12k
[b,a]=butter(5,[11000,13000]/(fs/2),'bandpass');
filtsig=filter(b,a,data);  %filtered signal
maxkenh1 = max(filtsig(:,1));
maxkenh2 = max(filtsig(:,2));
heso = maxkenh2/maxkenh1; %Tim do lech bien do giua 2 kenh
y1 = heso*filtsig(:,1);
y2 = filtsig(:,2);
y = [ y1 y2 ];%Ghep 2 kenh 
subplot(3,1,2);
plot(y); 



% %=====================
% %test tinh khoang cach
% %=====================
% y1 = abs(y1);
%   y1=envelope(y1,150,'rms'); %lay mau 150 mau/lan
%   y1=y1/(max(y1));
%   subplot(3,1,3);
%   plot(y1(1:100000));
%   xlabel('Time(s)');
%   ylabel('Amplitude(V)');
%   title('ENVELOPE SIGNAL');
% 
%  i = 1;
%  stt = 0;
%  D1=[];
% while i<=(length(y1)/100)
%     switch stt
%         case 0
%             if(y1(i)>= 0.03) %=== default = 0.3
%                 index1 = i;
%                 i = i+800;
%                 stt = 1;
%             else
%                 i = i+1;
%                 stt = 0;
%             end
%         case 1
%             if(y1(i)>= 0.002)%== dafault = 0.02
%                 index2 = i;
%                 i = i+200;
%                 time = index2-index1;
%                 distance=time*1500/2;
% %                 disp('Time');
% %                 disp(time*1/fs);
%                 disp('Distance(m)');
%                 distance=distance*1/fs;
%                 disp(distance);
%                 stt = 0;
%             else
%                 i = i+1;
%                 stt = 1;
%             end
%         otherwise
%             i = i+1;
%             stt = 0;
%     end
%    %D1=[D1 distance];
% end
