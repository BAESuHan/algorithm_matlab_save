clear;
close all;
%%
%read raw data


fs = 60;
cutoff_freq=[3 15];

filename_line = 'C:\Users\qotng\Desktop\GPG\매틀랩\19.12.31코드보완\json_file\tremorquantification-Line RowData_test-export (30).json';

number=1;

%read json file
fid = fopen(filename_line, 'r');

str = fread(fid,'*char').';

fclose(fid);

J = jsondecode(str);

line_time=J.TaskNo00.time';
line_xp= J.TaskNo00.x_position';
line_yp=J.TaskNo00.y_position';


an_bpf_x = J.testing.BPF_x';
an_bpf_y = J.testing.BPF_y';
an_lpf_x = J.testing.LPF_x';
an_pca = J.testing.PCA';
an_hil = J.testing.HIL';
an_fft = J.testing.FFT';
an_measure = J.testing.MEASURE';

B_TM = an_measure(1);
B_TF = an_measure(2);
B_Time = an_measure(3);
B_ED = an_measure(4);
B_Velocity = an_measure(5);

line_time = line_time/1000;


if(length(line_xp) > length(line_yp))
    line_xp = line_xp(1,1:end-1);
    line_time = line_time(1,1:end-1);
end 



%%
%%Filtering
%lowpass(xp_line,3,fs);
%bandpass(xp_line,cutoff_freq,fs);

%Line
line_low_x=lowpass(line_xp,1,fs)-480;
line_low_y=lowpass(line_yp,1,fs)-(1748/2);
[line_band_x,~] = bandpass(line_xp,cutoff_freq,fs);
[line_band_y,~] = bandpass(line_yp,cutoff_freq,fs);



%Line PCA
line_xy = vertcat(line_band_x,line_band_y);
temp_pca = myPCA(line_xy);
new_l= temp_pca(1,4:end);

Csig = line_xy * line_xy'; % Covariance
[V,D] = eig(Csig); % eigenvector
W = V(:,end); % weight 
projection = W'*line_xy;  % projection



%%
%%TM

%Line
[line_envelope,~]= envelope(new_l);
A_TM=rms(line_envelope)/88.18;



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

%%peak값은 Hz
A_TF=line_f(line_index_fft_result);

%%
%%Line error distance
% line_err_distance=rms(line_low_x);
% A_ED = line_err_distance/88.18;


%%
%%time
A_Time = line_time(1,end); 
list_time = line_time;

%%
%%velocity

%Line


for i=1:length(line_xp)
   base_l(1,i) = 0;  
end
for i= 1:length(line_xp)
   line_error_d11(1,i) =sqrt((line_low_x(1,i)-base_l(1,i)).^2) ;
end

%ED

line_error_d=0;
for i= 11:length(line_xp)-7
   line_error_d = line_error_d +sqrt((line_low_x(1,i)-base_l(1,i)).^2) ;
end
A_ED = line_error_d/length(line_xp)/88.18;

% 실제로 그린 line의 길이
real_dis_line=0;
for i= 11:length(line_xp)-7
   real_dis_line =real_dis_line+sqrt((line_low_x(1,i)-line_low_x(1,i-1)).^2) ;
end
A_line_draw_length = real_dis_line/88.18;

A_Velocity = A_line_draw_length/A_Time;

B_an_vel = B_Velocity/A_Velocity*100;

%%
%%plot
% 
% %Line
figure();
% plot(line_base_y,line_base_x);
% hold on
subplot(2,3,1);
plot(line_yp,line_xp,'Linewidth',2);
%  plot(list_time,line_xp);
axis equal;
title('Line raw data')

subplot(2,3,2);
plot(line_band_x);
hold on
plot(an_bpf_x);
title('Line BPF_x ')

subplot(2,3,3);
plot(line_band_y);
hold on
plot(an_bpf_y);
title('Line BPF_y')

subplot(2,3,4);
plot(new_l);
hold on
plot(an_pca);
title('Line PCA ')

subplot(2,3,5);
plot(line_envelope);
hold on
plot(an_hil(1,4:end));
title('Line Hilbert ')

subplot(2,3,6);
plot(line_low_x);
hold on
plot(an_lpf_x);
title('Line LPF_x ')



%%
%랩미팅용
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

