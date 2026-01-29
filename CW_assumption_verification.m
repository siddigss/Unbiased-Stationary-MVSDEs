
read_data = true;
if read_data

folder_name = 'results_6/particles_1';
folder_struct = dir(folder_name);       % files will contain '.' and '..', so subtract 2 later.
files_count = length(folder_struct)-2;
particles_1 = [];
indices_order = zeros(1,files_count);
for file_index = 1:(length(folder_struct)-2)
    file_name = strcat(folder_struct(file_index+2).folder, '\', folder_struct(file_index+2).name);
    A = readmatrix(file_name);
    particles_1 = [particles_1, A];
    index_string = file_name(end-5:end-4);
    if index_string(1) == '_'
        index_string = index_string(2);
    end
    indices_order(file_index) = str2num(index_string);
end

folder_name = 'results_6/particles_2';
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


folder_name = 'results_6/ls';
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


folder_name = 'results_6/ps';
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
levels = ls;
l_min = 3;
l_max = 10;
p_min = 0;
p_max = 7;
T_real = 101-10;
M = length(ls);
step = 0.1;
points = -3:step:3;
density_est = zeros(1, length(points));
coupled_levels = levels(levels ~= l_min);
coupled_level_counter = 1;


% densities
l_density = 2.^(-(l_min:l_max)).* ((l_min:l_max)+1) .*log2((l_min:l_max)+2);
l_density = l_density/sum(l_density);
l_density_f = @(l) l_density(l-l_min+1);
p_density = 2.^(-(p_min:p_max)).* ((p_min:p_max)+1) .*log2((p_min:p_max)+2).^2;
p_density = p_density/sum(p_density);
p_density_f = @(p) p_density(p-p_min+1);

phi = @(x) x.^2;
particles_1_counter = 1;
particles_2_counter = 1;
phi_pre_estimatros = zeros(1, length(ls));
for j = 1:length(ls)
    a = 0;
    I_p = 10*2^ps(j);
    I_p_minus_1 = 5*2^(ps(j));
    particles_1_local = particles_1(particles_1_counter:particles_1_counter+I_p);
    particles_1_counter = particles_1_counter + I_p+1;
    if ls(j) == l_min && ps(j) == p_min
        a = mean(phi(particles_1_local));
    end
    if ls(j) == l_min && ps(j) > p_min
        a = mean(phi(particles_1_local)) - mean(phi(particles_1_local(1:I_p_minus_1)));
    end
    if ls(j) > l_min && ps(j) == p_min
        particles_2_local = particles_2(particles_2_counter:particles_2_counter+I_p);
        particles_2_counter = particles_2_counter + I_p+1;
        a = mean(phi(particles_1_local)) - mean(phi(particles_2_local));
    end
    if ls(j) > l_min && ps(j) > p_min
        particles_2_local = particles_2(particles_2_counter:particles_2_counter+I_p);
        particles_2_counter = particles_2_counter + I_p+1;
        a = mean(phi(particles_1_local)) - mean(phi(particles_1_local(1:I_p_minus_1))) - mean(phi(particles_2_local)) + mean(phi(particles_2_local(1:I_p_minus_1)));
    end
    phi_pre_estimatros(j) = a;
end
assumption_quantity_l_p = zeros(l_max-l_min+1,p_max-p_min+1);
counter_l_p = zeros(l_max-l_min+1,p_max-p_min+1);
for j = 1:length(ls)
    l_index = ls(j) - l_min + 1;
    p_index = ps(j)+1;
    if counter_l_p(l_index, p_index) < 100
        assumption_quantity_l_p(l_index, p_index) = assumption_quantity_l_p(l_index, p_index) + abs(phi_pre_estimatros(j));
        counter_l_p(l_index, p_index) = counter_l_p(l_index, p_index) + 1;
    end
end

assumption_quantity_1 = zeros(1,l_max-l_min+1);
assumption_quantity_2_aux = zeros(1,l_max-l_min+1);
for i = 1:(l_max-l_min+1)
    assumption_quantity_l_local = assumption_quantity_l_p(i,:)./counter_l_p(i,:);
    ratio_matrix = assumption_quantity_l_local(2:end)./assumption_quantity_l_local(1:end-1);
    assumption_quantity_1(i) = max(ratio_matrix(end:end));
    assumption_quantity_2_aux(i) = sum(assumption_quantity_l_local);
end
assumption_quantity_2 = assumption_quantity_2_aux(2:end)./assumption_quantity_2_aux(1:end-1);

figure;
plot(l_min:l_max, assumption_quantity_1, 'o-', 'LineWidth',2, 'MarkerSize', 12);
hold on;
plot(l_min:l_max-1, assumption_quantity_2, 'x-', 'LineWidth',2, 'MarkerSize', 12);
title('Assumption (A3) Verification for Curie-Weiss Model');
fontsize(14, "points");
xlim([2.5,10.5]);
ylim([0.1,0.82]);
xlabel('Level');
ylabel('Estimated Value');
lgs = legend('$\widehat{Q}_1$', '$\widehat{Q}_2$');
lgs.Location = 'southeast';
hold off;
%export_fig figs/CW_assumption.pdf -q101;



function d = level_density(l, l_min, l_max)
    unnorm_density = 2.^(-1.5*(l_min:l_max));
    density = unnorm_density/sum(unnorm_density);
    d = density(l-l_min+1);
end
