


folder_name = 'samples';
folder_struct = dir(folder_name);       % files will contain '.' and '..', so subtract 2 later.
files_count = length(folder_struct)-2;
samples = [];
for file_index = 1:(length(folder_struct)-2)
    file_name = strcat(folder_struct(file_index+2).folder, '\', folder_struct(file_index+2).name);
    A = readmatrix(file_name);
    samples = [samples, A];
end


folder_name = 'times';
folder_struct = dir(folder_name);       % files will contain '.' and '..', so subtract 2 later.
files_count = length(folder_struct)-2;
sample_count_per_file = 1000;
samples_count = files_count*sample_count_per_file;
times = zeros(1, files_count*sample_count_per_file);
for file_index = 1:(length(folder_struct)-2)
    file_name = strcat(folder_struct(file_index+2).folder, '\', folder_struct(file_index+2).name);
    A = readmatrix(file_name);
    times(:,(1+(file_index-1)*sample_count_per_file):(file_index*sample_count_per_file)) = A;
end





%%%%%%%%%% MSE
y_data = [0.793827101386834;
   2.812652893885026;
   0.469435277216595;
   2.446576623207243;
   0.569933921486705;
   5.155285594063845;
   3.351453978485524;
  -0.880627322676377;
   0.941047112519439;
  -0.605404583359629];


true_value = mean(y_data);
samples_count = length(samples);
Ms = 2.^(3:11);
M_count = length(Ms);
MSEs = zeros(1, M_count);
times_avrg = zeros(1, M_count);
for i = 1:M_count
    M = Ms(i);
    samples_count_tmp = min(samples_count, M*50);
    estimators = zeros(1, floor(samples_count_tmp/M));
    times_tmp = zeros(1, floor(samples_count_tmp/M));
    for j = 0:floor(samples_count_tmp/M)-1
        estimators(j+1) = mean(samples(1,(j*M+1):((j+1)*M)));
        times_tmp(j+1) = sum(times(1,(j*M+1):((j+1)*M)));
    end
    MSEs(i) = mean((estimators-true_value).^2);
    times_avrg(i) = mean(times_tmp);
end

% line_info = polyfit(log2(Ms), log2(MSEs), 1);
% figure;
% plot(log2(Ms),log2(MSEs), '*', 'LineWidth',2)
% hold on
% plot(log2(Ms), line_info(1)*log2(Ms) + line_info(2), '-', 'LineWidth',2);
% legend(['slope = ', num2str(line_info(1))]);
% xlim([2.5;11.5]);
% xlabel('$\log_2$M')
% ylabel('$\log_2$MSE');
% title('MSE for Parameter Estimation');
% fontsize(15, "points");
% %set(gcf,'Color','w');
% hold off;

%%%%%%%%%% Time Histogram
% 
% figure;
% histogram(samples, 'Normalization','probability');
% hold on
% title('Meeting Time for Parameter Estimation');
% fontsize(14, "points");
% %set(gcf,'Color','w');
% xlabel('Meeting Time');
% ylabel('Percentage');
% hold off;

%%%%%%%%%% Level Histogram

% figure;
% histogram(samples, 'Normalization','probability');
% hold on
% title('Levels for Curie-Weiss');
% fontsize(14, "points");
% set(gcf,'Color','w');
% hold off;



line_info = polyfit(log2(Ms), log2(times_avrg), 1);
figure;
plot(log2(Ms),log2(times_avrg), '*', 'LineWidth',2)
hold on
plot(log2(Ms), line_info(1)*log2(Ms) + line_info(2), '-', 'LineWidth',2);
legend(['slope = ', num2str(line_info(1))]);
xlim([2.5;11.5]);
ylim([3,13]);
xlabel('$\log_2$M')
ylabel('$\log_2$Time');
title('Average Run Time for Parameter Estimation');
fontsize(15, "points");
%set(gcf,'Color','w');
hold off;
