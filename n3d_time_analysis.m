
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


line_info = polyfit(log2(Ms), log2(times_avrg), 1);
figure;
plot(log2(Ms),log2(times_avrg), '*', 'LineWidth',2)
hold on
plot(log2(Ms), line_info(1)*log2(Ms) + line_info(2), '-', 'LineWidth',2);
legend(['slope = ', num2str(line_info(1))]);
xlim([2.5;11.5]);
ylim([2,11]);
xlabel('$\log_2$M')
ylabel('$\log_2$Time');
title('Average Run Time for 3D Neuron');
fontsize(15, "points");
%set(gcf,'Color','w');
hold off;
% export_fig figs/toy_mean_time.pdf -q101;


