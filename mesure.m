clear;
close all;
%%
%read raw data


fs = 60;
cutoff_freq=[3 15];

filename_line = 'C:\Users\qotng\Desktop\GPG\HCI_³í¹®ÁØºñ\¸ÅÆ²·¦ ÁØºñ\¸ÅÆ²·¦ÄÚµå\test.xlsx';

number=1;

%Line raw data
line_time = xlsread(filename_line, number, 'B2:CDX2');
line_xp= xlsread(filename_line, number, 'B3:CDX3'); %300°³ 
line_yp = xlsread(filename_line, number, 'B4:CDX4');

line_time = line_time(1,1:end)/1000;
line_xp = line_xp(1,1:end);

if(length(line_xp) > length(line_yp))
    line_xp = line_xp(1,1:end-1);
    line_time = line_time(1,1:end-1);
end 


%%
%Spiral raw data

% 1.1 make data form
%normal data

[~, txt_base_x] = xlsread('C:\Users\qotng\Desktop\GPG\¸ÅÆ²·¦\½ºÆÄÀÌ·²\normal_data.xlsx', 'B2:B11');
[~, txt_base_y] = xlsread('C:\Users\qotng\Desktop\GPG\¸ÅÆ²·¦\½ºÆÄÀÌ·²\normal_data.xlsx', 'C2:C11');
[~, txt_pos_x] = xlsread('C:\Users\qotng\Desktop\GPG\¸ÅÆ²·¦\½ºÆÄÀÌ·²\normal_data.xlsx', 'D2:D11');
[~, txt_pos_y] = xlsread('C:\Users\qotng\Desktop\GPG\¸ÅÆ²·¦\½ºÆÄÀÌ·²\normal_data.xlsx', 'E2:E11');

x_base_string = string(txt_base_x);
y_base_string = string(txt_base_y);
x_pos_string = string(txt_pos_x);    
y_pos_string = string(txt_pos_y);
    
%sorting data(N)
   
for idx = 1:length(txt_base_x);
    x_base_split = strsplit(x_base_string(idx,:), ':');
    y_base_split = strsplit(y_base_string(idx,:), ':');
    x_pos_split = strsplit(x_pos_string(idx,:), ':');
    y_pos_split = strsplit(y_pos_string(idx,:), ':');
    x_base_num = str2double(x_base_split);
    y_base_num = str2double(y_base_split);
    x_pos_num = str2double(x_pos_split);
    y_pos_num = str2double(y_pos_split);
    if idx == 1
        task00 = [x_base_num; y_base_num; x_pos_num; y_pos_num];
    else if idx == 2
        task01 = [x_base_num; y_base_num; x_pos_num; y_pos_num];
    else if idx == 3
        task02 = [x_base_num; y_base_num; x_pos_num; y_pos_num];
    else if idx == 4
        task03 = [x_base_num; y_base_num; x_pos_num; y_pos_num];
    else if idx == 5
        task04 = [x_base_num; y_base_num; x_pos_num; y_pos_num];
    else if idx == 6
        task05 = [x_base_num; y_base_num; x_pos_num; y_pos_num];
    else if idx == 7
        task06 = [x_base_num; y_base_num; x_pos_num; y_pos_num];
    else if idx == 8
        task07 = [x_base_num; y_base_num; x_pos_num; y_pos_num];
    else if idx == 9
        task08 = [x_base_num; y_base_num; x_pos_num; y_pos_num];
    else if idx == 10
        task09 = [x_base_num; y_base_num; x_pos_num; y_pos_num];
            end
            end
            end
            end
            end
            end
            end
        end
        end
    end
    end

    %abnormal data
    [~, txt_base_x] = xlsread('C:\Users\qotng\Desktop\GPG\¸ÅÆ²·¦\½ºÆÄÀÌ·²\abnormal_data.xlsx', 'B2:B11');
    [~, txt_base_y] = xlsread('C:\Users\qotng\Desktop\GPG\¸ÅÆ²·¦\½ºÆÄÀÌ·²\abnormal_data.xlsx', 'C2:C11');
    [~, txt_pos_x] = xlsread('C:\Users\qotng\Desktop\GPG\¸ÅÆ²·¦\½ºÆÄÀÌ·²\abnormal_data.xlsx', 'D2:D11');
    [~, txt_pos_y] = xlsread('C:\Users\qotng\Desktop\GPG\¸ÅÆ²·¦\½ºÆÄÀÌ·²\abnormal_data.xlsx', 'E2:E11');
    
    x_base_string = string(txt_base_x);
    y_base_string = string(txt_base_y);
    x_pos_string = string(txt_pos_x);
    y_pos_string = string(txt_pos_y);
    
    %sorting data(Ab)
    for idx = 1:length(txt_base_x);
x_base_split = strsplit(x_base_string(idx,:), ':');
y_base_split = strsplit(y_base_string(idx,:), ':');
x_pos_split = strsplit(x_pos_string(idx,:), ':');
y_pos_split = strsplit(y_pos_string(idx,:), ':');
x_base_num = str2double(x_base_split);
y_base_num = str2double(y_base_split);
x_pos_num = str2double(x_pos_split);
y_pos_num = str2double(y_pos_split);
    if idx == 1
        task10 = [x_base_num; y_base_num; x_pos_num; y_pos_num];
    else if idx == 2
        task11 = [x_base_num; y_base_num; x_pos_num; y_pos_num];
    else if idx == 3
        task12 = [x_base_num; y_base_num; x_pos_num; y_pos_num];
        else if idx == 4
        task13 = [x_base_num; y_base_num; x_pos_num; y_pos_num];
        else if idx == 5
        task14 = [x_base_num; y_base_num; x_pos_num; y_pos_num];
        else if idx == 6
        task15 = [x_base_num; y_base_num; x_pos_num; y_pos_num];
        else if idx == 7
        task16 = [x_base_num; y_base_num; x_pos_num; y_pos_num];
        else if idx == 8
        task17 = [x_base_num; y_base_num; x_pos_num; y_pos_num];
        else if idx == 9
        task18 = [x_base_num; y_base_num; x_pos_num; y_pos_num];
        else if idx == 10
        task19 = [x_base_num; y_base_num; x_pos_num; y_pos_num];
            end
            end
            end
            end
            end
            end
            end
        end
        end
    end
    end   
%%
 %1.2. select data
 x_base = task18(1,:);
 y_base = task18(2,:);
 x_position = task18(3,:);
 y_position = task18(4,:);
 


%%
%%Filtering
%lowpass(xp_line,3,fs);
%bandpass(xp_line,cutoff_freq,fs);

%Line
line_low_x=lowpass(line_xp,1,fs)-480;
line_low_y=lowpass(line_yp,1,fs)-(1748/2);
[line_band_x,~] = bandpass(line_xp,cutoff_freq,fs);
[line_band_y,~] = bandpass(line_yp,cutoff_freq,fs);

%Spiral
 %low_spiral=lowpass(x_position,3,fs)-480;
 [spiral_band_x,~] = bandpass(x_position,cutoff_freq,fs);
 [spiral_band_y,~] = bandpass(y_position,cutoff_freq,fs);

%Line PCA
line_xy = vertcat(line_band_x,line_band_y);
new_l=myPCA(line_xy);

Csig = line_xy * line_xy'; % Covariance
[V,D] = eig(Csig); % eigenvector
W = V(:,end); % weight 
projection = W'*line_xy;  % projection

line_low_xy = vertcat(line_low_x,line_low_y);
new_l_low=myPCA(line_low_xy); 

%Spiral PCA
spiral_xy = vertcat(spiral_band_x,spiral_band_y);
new = myPCA(spiral_xy);   

%%
%%TM

%Line
[line_envelope,~]= envelope(new_l);
A_line_TM=rms(line_envelope)/88.18;

%Spiral
 [spiral_envelope,~]= envelope(new);
 B_spiral_TM=rms(spiral_envelope)/88.18;

%%
%Line Frequency
[line_fft_result,line_f]=tremor_fft(new_l,fs);



for i=1:length(line_f)
    if(line_f(1,i) >= 10)
       tmp_first = line_f(1,i);
       tmp_index_first = i;
      
       break;
    end 
   
end

for i=1:length(line_f)
    if(line_f(1,i) <= 25)
       tmp_last = line_f(1,i);
       tmp_index_last = i;
    end
end
 

f_m = mean(line_fft_result(1,tmp_index_first:tmp_index_last));
n_st=sqrt(var(line_fft_result(1,tmp_index_first:tmp_index_last)));
 
for i=1:length(line_fft_result)
   f_mm(1,i) = f_m;  
   m_var1(1,i) = f_m+n_st;
   m_var2(1,i) = f_m+2*n_st;
   m_var3(1,i) = f_m+3*n_st;
   m_var4(1,i) = f_m+40*n_st;
end


figure();
plot(line_f,line_fft_result);
hold on
plot(line_f,f_mm);
hold on
plot(line_f,m_var1);
hold on
plot(line_f,m_var2);
hold on
plot(line_f,m_var3);
hold on
plot(line_f,m_var4);

hold off
[line_max_fft_result,line_index_fft_result]=max(line_fft_result);


% [line_max_fft_result,line_index_fft_result]=peak_function(new_l);

%%peak°ªÀº Hz
A_line_Freqency=line_f(line_index_fft_result);
%%
%Spiral Frequency
[spiral_fft_result,spiral_f]=tremor_fft(new,fs);

[max_spiral_fft_result,index_spiral_fft_result]=max(spiral_fft_result);

%%peak°ªÀº Hz
B_spiral_Freqency=spiral_f(index_spiral_fft_result);

%%
%%Line error distance
line_err_distance=rms(line_low_x);
A_line_err_cm = line_err_distance/88.18;

%%
%%Spiral error distance
[B_spiral_err_cm,low_spiral_x,low_spiral_y ]= err_spiral_cm(x_base,y_base,x_position,y_position);
%%
%%time
A_line_time = line_time(1,end); 
list_time = line_time;

%%
%%velocity

%Line

real_dis_line=0;
for i=1:length(line_xp)
   base_l(1,i) = 0;  
end
for i= 1:length(line_xp)
   line_error_d11(1,i) =sqrt((line_low_x(1,i)-base_l(1,i)).^2) ;
end

%ED

line_error_d=0;
for i= 1:length(line_xp)
   line_error_d = line_error_d +sqrt((line_low_x(1,i)-base_l(1,i)).^2) ;
end
A_line_Error_distance = line_error_d/length(line_xp)/88.18;

% ½ÇÁ¦·Î ±×¸° lineÀÇ ±æÀÌ

for i= 2:length(line_xp)
   real_dis_line =real_dis_line+sqrt((line_low_x(1,i)-line_low_x(1,i-1)).^2) ;
end
A_line_draw_length = real_dis_line/88.18;

A_line_velocity = A_line_draw_length/A_line_time;

%Spiral
r_dis_spiral=0;
for i= 2:length(x_position)
   r_dis_spiral =r_dis_spiral+sqrt((low_spiral_x(1,i)-low_spiral_x(1,i-1)).^2 +(low_spiral_y(1,i)-low_spiral_y(1,i-1)).^2) ;
end

B_spiral_real_dis = r_dis_spiral/88.18;

%B_spiral_velocity = B_spiral_real_dis/time

%%
%%plot
% 
% %Line
figure();
% plot(line_base_y,line_base_x);
% hold on
plot(line_yp,line_xp,'Linewidth',2);
%  plot(list_time,line_xp);
axis equal;
title('Line raw data')
% 
% % figure();
% % plot(line_low);
% % title('Line LPF ')
% 
% %Spiral
% figure();
% plot(x_base, y_base,'.', x_position, y_position ,'o');
% title('before find pair-Spiral');
%%
%·¦¹ÌÆÃ¿ë
% figure()
% plot(line_time,line_xp,'Linewidth',2)
% hold on 
% plot(line_time,line_yp,'Linewidth',2)
% legend('Raw - X axis','Raw - Y axis')

% plot(line_time,line_band_x,'Linewidth',2)
% hold on 
% plot(line_time,line_band_y,'Linewidth',2)
% legend('BPF - X axis','BPF - Y axis')

% plot(line_time,new_l,'Color',[0.7 0.1 0.2],'Linewidth',2)
% legend('PCA - X,Y axis')
 
% plot(line_f,line_fft_result,'Color',[0.2 0.1 0.7],'Linewidth',2);
% hold on
% plot(line_f,m_var4,'Color',[0.5 0.7 0.2],'Linewidth',1.5);
% legend('FFT - PCA(X,Y)','mean + 70*variance')

% plot(line_time,abs(hilbert(new_l)),'Color',[0.8 0.8 0.1],'Linewidth',2)
% legend('Absolute of Hilbert - X,Y axis')

% plot(line_time,line_low_x,'Linewidth',2)
% legend('LPF 1Hz - X axis')

% plot(line_time,r_dis_line11,'Linewidth',2)
% legend('L2err(t)  - X axis')

