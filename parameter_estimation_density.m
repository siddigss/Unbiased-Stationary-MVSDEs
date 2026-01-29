
read_data = true;
if read_data

folder_name = 'particles_1';
folder_struct = dir(folder_name);       % files will contain '.' and '..', so subtract 2 later.
files_count = length(folder_struct)-2;
particles_1 = [];
indices_order = zeros(1,files_count);
for file_index = 1:(length(folder_struct)-2)
    file_name = strcat(folder_struct(file_index+2).folder, '\', folder_struct(file_index+2).name);
    A = readmatrix(file_name);
    if size(A,1) ~= 11
        disp(file_name);
    end
    %ps = [ps, [A(1,:);A(12,:)]];
    particles_1 = [particles_1, A];
    index_string = file_name(end-5:end-4);
    if index_string(1) == '_'
        index_string = index_string(2);
    end
    indices_order(file_index) = str2num(index_string);
end 


folder_name = 'particles_2';
folder_struct = dir(folder_name);       % files will contain '.' and '..', so subtract 2 later.
default_name = folder_struct(file_index+2).name;
default_name = default_name(1:end-5);
if default_name(end) ~= '_'
    default_name = default_name(1:end-1);
end
particles_2 = [];
for file_index = 1:(length(folder_struct)-2)
    file_name = strcat(folder_struct(file_index+2).folder, '\', folder_struct(file_index+2).name);
    A = readmatrix(file_name);
    particles_2 = [particles_2, A];
end


folder_name = 'ls';
folder_struct = dir(folder_name);       % files will contain '.' and '..', so subtract 2 later.
default_name = folder_struct(file_index+2).name;
default_name = default_name(1:end-5);
if default_name(end) ~= '_'
    default_name = default_name(1:end-1);
end
ls = [];
for file_index = 1:(length(folder_struct)-2)
    file_name = strcat(folder_struct(file_index+2).folder, '\', folder_struct(file_index+2).name);
    A = readmatrix(file_name);
    ls = [ls, A];
end

folder_name = 'ps';
folder_struct = dir(folder_name);       % files will contain '.' and '..', so subtract 2 later.
default_name = folder_struct(file_index+2).name;
default_name = default_name(1:end-5);
if default_name(end) ~= '_'
    default_name = default_name(1:end-1);
end
ps = [];
for file_index = 1:(length(folder_struct)-2)
    file_name = strcat(folder_struct(file_index+2).folder, '\', folder_struct(file_index+2).name);
    A = readmatrix(file_name);
    ps = [ps, A];
end




end
%%%%%%%%%% Density Estimator
l_min = 3;
l_max = 10;
p_min = 0;
p_max = 7;
T_real = 101-10;
M = length(ls);
step = 0.1;
points = -1:step:4;
density_est = zeros(1, length(points));
coupled_levels = levels(levels ~= l_min);
coupled_level_counter = 1;


l_density = 2.^(-(l_min:l_max)).* ((l_min:l_max)+1) .*log2((l_min:l_max)+2);
l_density = l_density/sum(l_density);
l_density_f = @(l) l_density(l-l_min+1);
p_density = 2.^(-(p_min:p_max)).* ((p_min:p_max)+1) .*log2((p_min:p_max)+2).^2;
p_density = p_density/sum(p_density);
p_density_f = @(p) p_density(p-p_min+1);



% kernel density estimator
kernel = @(x) normpdf(x,0,1);
kernel_window = 0.35;
kernel_step = 0.1;
kernel_points = -2:kernel_step:3;
kernel_density_est = zeros(1, length(kernel_points));

for i = 1:length(kernel_points)
    disp([num2str(i),'/',num2str(length(kernel_points))]);
    phi = @(x) kernel((kernel_points(i) - x(11,:))/kernel_window)/kernel_window;
    %phi = @(x) x(1,:);

    e = 0;
    particles_1_counter = 1;
    particles_2_counter = 1;
    for j = 1:length(ls)
        a = 0;
        I_p = 10*2^ps(j);
        I_p_minus_1 = 5*2^(ps(j));
        particles_1_local = particles_1(:,particles_1_counter:particles_1_counter+I_p);
        particles_1_counter = particles_1_counter + I_p;
        if ls(j) == l_min && ps(j) == p_min
            a = mean(phi(particles_1_local));
        end
        if ls(j) == l_min && ps(j) > p_min
            a = mean(phi(particles_1_local)) - mean(phi(particles_1_local(:,1:I_p_minus_1)));
        end
        if ls(j) > l_min && ps(j) == p_min
            particles_2_local = particles_2(:,particles_2_counter:particles_2_counter+I_p);
            particles_2_counter = particles_2_counter + I_p;
            a = mean(phi(particles_1_local)) - mean(phi(particles_2_local));
        end
        if ls(j) > l_min && ps(j) > p_min
            particles_2_local = particles_2(:,particles_2_counter:particles_2_counter+I_p);
            particles_2_counter = particles_2_counter + I_p;
            a = mean(phi(particles_1_local)) - mean(phi(particles_1_local(:,1:I_p_minus_1))) - mean(phi(particles_2_local)) + mean(phi(particles_2_local(:,1:I_p_minus_1)));
        end
        e = e + a/(l_density_f(ls(j))*p_density_f(ps(j)));
    end
    kernel_density_est(i) = e;
end
kernel_density_est = max(kernel_density_est,0)/M;

true_density_function = @(x) normpdf(x, (mean(y_data)+y_data(10))/2, 0.5);

figure;
plot(kernel_points, kernel_density_est, 'o-', 'LineWidth',2);
hold on;
plot(kernel_points, true_density_function(kernel_points), '-', 'LineWidth',2);
%title('Invariant Distribution for Parameter Estimation');
title('Invariant Distribution of the 10th component of $X_t$');
fontsize(14, "points");
%xlabel('$\theta$');
xlabel('$x$');
%ylabel('Estimated Density');
ylabel('Density');
legend('Estimated Posterior', 'Exact Posterior');
hold off;


