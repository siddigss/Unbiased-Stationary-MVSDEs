

job_id = '3';
proc_id = '3';
folder_read = '';
folder_write = sprintf('%s%s', job_id, '/');
results_filename = sprintf('%ssamples_%s_%s.txt', folder_write, job_id, proc_id);
results_filename_3 = sprintf('%sparticles_1_%s_%s.txt', folder_write, job_id, proc_id);
results_filename_4 = sprintf('%sparticles_2_%s_%s.txt', folder_write, job_id, proc_id);
results_filename_6 = sprintf('%slevels_%s_%s.txt', folder_write, job_id, proc_id);
results_filename_7 = sprintf('%sps_%s_%s.txt', folder_write, job_id, proc_id);
results_filename_8 = sprintf('%stimes_%s_%s.txt', folder_write, job_id, proc_id);

T = 1000;
alpha = 0.9;
theta_true = 1;
d = 3;
X0 = [ones(d,1)];

params = [0, 0.5, 0.3, 0.4, 0.4, 0.05];

%phi = @(x) x;
phi = @(x) x(1);

l_min = 3;
l_max = 10;
p_min = 0;
p_max = 7;
sample_count = 1000;
samples = zeros(1,sample_count);
meeting_times_single = [];
meeting_times_coupled_1 = [];
meeting_times_coupled_2 = [];
meeting_times_max = zeros(1,sample_count);
ls = zeros(1,sample_count);
ps = zeros(1,sample_count);
particles_1 = [];
particles_2 = [];
times = zeros(1,sample_count);

for i = 1:sample_count
    [level, l_density] = sample_l(l_max, l_min);
    [p, p_density] = sample_p(p_max, p_min);
    ls(i) = level;
    ps(i) = p;
    delta_t = 2^(-level);
    disc_count = 2^(level);
    law_particles_count = round(10*level);
    I_p = 10*2^p;
    I_p_minus_1 = 10*2^(p-1);

    disp(['level = ', num2str(level)]);
    disp(i);
    iter_time_start = tic;
    if level == l_min

        X_0_array = zeros(d,law_particles_count,1);
        X_0_array(:,:,1) = sample_X0(params, law_particles_count);
        X_unit = zeros(3,I_p+1);
        X_unit(:,1) = X0;
        for j = 1:I_p 
            X_law = simulate_law_model(delta_t, 1, X_0_array, params);
            delta_W = sqrt(delta_t)*randn(d, 1, 1*disc_count);
            X_unit(:,j+1) = simulate_model_given_noise_and_law(delta_t, 1, X_unit(:,j), params, delta_W, X_law);
            X_0_array(:,:,1) = X_law(:,:,end);
        end
        if p ==0
            samples(i) = 1/(l_density*p_density)*mean(phi(X_unit));
        else
            samples(i) = 1/(l_density*p_density)*(mean(phi(X_unit)) - mean(phi(X_unit(:,1:I_p_minus_1))));
        end
        particles_1 = [particles_1, X_unit];        
    end

    if level > l_min
		X_1_0_array = zeros(d,law_particles_count,1);
		X_2_0_array = zeros(d,law_particles_count,1);
        X_1_0_array(1,:,:) = X0(1);
        X_2_0_array(1,:,:) = X0(1);
        X_1_0_array(2:end,:,:) = X0(2);
        X_2_0_array(2:end,:,:) = X0(2);
		X_1_unit = zeros(3,I_p+1);
        X_2_unit = zeros(3,I_p+1);
        X_1_unit(:,1) = X0;
        X_2_unit(:,1) = X0;
        for j = 1:I_p
			delta_W_1 = sqrt(delta_t)*randn(d, law_particles_count, 1/delta_t);
			delta_W_2 = delta_W_1(:,:,1:2:end) + delta_W_1(:,:,2:2:end);
			X_law_1 = simulate_law_model_given_noise(delta_t, 1, X_1_0_array, params, delta_W_1);
			X_law_2 = simulate_law_model_given_noise(2*delta_t, 1, X_2_0_array, params, delta_W_2);
			delta_W_1 = sqrt(delta_t)*randn(d, 1, 1*disc_count);
			delta_W_2 = delta_W_1(:,:,1:2:end) + delta_W_1(:,:,2:2:end);
            X_1_unit(:,j+1) = simulate_model_given_noise_and_law(delta_t, 1 , X_1_unit(:,j), params, delta_W_1, X_law_1);
            X_2_unit(:,j+1) = simulate_model_given_noise_and_law(2*delta_t, 1, X_2_unit(:,j), params, delta_W_2, X_law_2);
			X_1_0_array(:,:,1) = X_law_1(:,:,end);
			X_2_0_array(:,:,1) = X_law_2(:,:,end);
        end
        if p ==0
            samples(i) = 1/(l_density*p_density)*(mean(phi(X_1_unit))-mean(phi(X_2_unit)));
        else
            samples(i) = 1/(l_density*p_density)*(mean(phi(X_1_unit)) - mean(phi(X_1_unit(:,1:I_p_minus_1))) - mean(phi(X_2_unit)) + mean(phi(X_2_unit(:,1:I_p_minus_1))));
        end
        particles_1 = [particles_1, X_1_unit];        
        particles_2 = [particles_2, X_2_unit];              
    end
    iter_time_end = toc(iter_time_start);
    times(i) = iter_time_end;
end

disp(mean(samples));

% writematrix(samples, results_filename);
% writematrix(particles_1, results_filename_3);
% writematrix(particles_2, results_filename_4);
% writematrix(ls, results_filename_6);
% writematrix(ps, results_filename_7);
% writematrix(times, results_filename_8);







% sampling functions
function [l, l_density] = sample_l(l_max, l_min)
    levels = l_min:l_max;
    density_unnorm = 2.^(-levels) .* (levels+1) .*log2(levels+2);
    density = density_unnorm/sum(density_unnorm);
    l = randsample(levels,1,true,density);
    l_density = density(l-l_min+1);
end

function [p, p_density] = sample_p(p_max, p_min)
    levels = p_min:p_max;
    density_unnorm = 2.^(-levels) .* (levels+1) .*log2(levels+2);
    density = density_unnorm/sum(density_unnorm);
    p = randsample(levels,1,true,density);
    p_density = density(p-p_min+1);
end



% the method functions
function [X_1_end, X_2_end] = sample_Q_l_caron(delta_t, X0_1, X0_2, params, laws)
    T = 1;
    particles_count = 1;
    d = size(X0_1,1);
    steps_count = T/delta_t;
    delta_W = sqrt(delta_t)*randn(d, particles_count, steps_count);
    X_1_end = simulate_model_given_noise_and_law(delta_t, T, X0_1, params, delta_W, laws);
    X_2_end = simulate_model_given_noise_and_law(delta_t, T, X0_2, params, delta_W, laws);
end

function [X_1_end, X_2_end] = sample_R_l_caron(delta_t, X0_1, X0_2, params, laws)
    T = 1-delta_t;
    particles_count = 1;
    steps_count = T/delta_t;
    d = size(X0_1,1);
    delta_W = sqrt(delta_t)*randn(d, particles_count, steps_count);
    X_1_end = simulate_model_given_noise_and_law(delta_t, T, X0_1, params, delta_W, laws);
    X_2_end = simulate_model_given_noise_and_law(delta_t, T, X0_2, params, delta_W, laws);
    [X_1_end, X_2_end] = maximal_coupling_single_step_law(delta_t, X_1_end, X_2_end, params, laws(:,steps_count), laws(:,steps_count));
end

function [X_1_end, X_2_end] = sample_K_l_caron(delta_t, X0_1, X0_2, params, alpha, laws)
    a = rand();
    if a < alpha
        [X_1_end, X_2_end] = sample_Q_l_caron(delta_t, X0_1, X0_2, params, laws);
        return;
    end
    [X_1_end, X_2_end] = sample_R_l_caron(delta_t, X0_1, X0_2, params, laws);
end

function [X_1, X_2] = maximal_coupling_single_step_law(delta_t, X0_1, X0_2, params, X_law_1, X_law_2)
    [mean_1, mean_2, var_1, var_2] = calculate_mean_var_single_level_law(delta_t, X0_1, X0_2, params, X_law_1, X_law_2);
    [X_1, X_2] = maximal_coupling_normal(mean_1, mean_2, var_1, var_2);
end


function [X_1, X_2] = maximal_coupling_normal(mean_1, mean_2, var_1, var_2)
    X_1 = mvnrnd(mean_1, var_1)';
    W = rand() * mvnpdf(X_1, mean_1, var_1);
    if W < mvnpdf(X_1, mean_2, var_2)
        X_2 = X_1;
        return;
    end
    while true
        X_2 = mvnrnd(mean_2, var_2)';
        W = rand() * mvnpdf(X_2, mean_2, var_2);
        if W > mvnpdf(X_2, mean_1, var_1)
            return;
        end
    end
end

function [means_1, means_2, var_1, var_2] = calculate_mean_var_single_level_law(delta_t, x_1, x_2, params, X_law_1, X_law_2)
    means_1 = x_1 + drift(params, x_1, X_law_1) * delta_t;
    means_2 = x_2 + drift(params, x_2, X_law_2) * delta_t;
    var_1 = diffusion(params, x_1, X_law_1, delta_t) * diffusion(params, x_1, X_law_1, delta_t)' * delta_t;
    var_2 = diffusion(params, x_2, X_law_2, delta_t) * diffusion(params, x_2, X_law_2, delta_t)' * delta_t;
end


function b = equals_floats(x,y)
    if max(abs(x-y)) < 1e-20
        b = true;
        return;
    end
    b = false;
end
 


% coulped
function [X_1_end, X_2_end, X_3_end, X_4_end] = sample_Q_l_bar(delta_t, X0_1, X0_2, X0_3, X0_4, params, laws_1, laws_2)
    T = 1;
    d = size(X0_1, 1);
    particles_count = 1;
    steps_count = T/delta_t;
    delta_W = sqrt(delta_t)*randn(d, particles_count, steps_count);
    X_1_end = simulate_model_given_noise_and_law(delta_t, T, X0_1, params, delta_W, laws_1);
    X_2_end = simulate_model_given_noise_and_law(delta_t, T, X0_2, params, delta_W, laws_1);
    delta_W_2 = delta_W(:,:,1:2:end) + delta_W(:,:,2:2:end);
    X_3_end = simulate_model_given_noise_and_law(2*delta_t, T, X0_3, params, delta_W_2, laws_2);
    X_4_end = simulate_model_given_noise_and_law(2*delta_t, T, X0_4, params, delta_W_2, laws_2);
end

function [X_1_end, X_2_end, X_3_end, X_4_end] = sample_R_l_bar(delta_t, X0_1, X0_2, X0_3, X0_4, params, laws_1, laws_2)
    T = 1-delta_t;
    particles_count = 1;
    d = size(X0_1, 1);
    steps_count = T/delta_t;
    delta_W = sqrt(delta_t)*randn(d, particles_count, steps_count);
    X_1_end = simulate_model_given_noise_and_law(delta_t, T, X0_1, params, delta_W, laws_1);
    X_2_end = simulate_model_given_noise_and_law(delta_t, T, X0_2, params, delta_W, laws_1);
    delta_W_2 = delta_W(:,:,1:2:end-1) + delta_W(:,:,2:2:end-1);
    X_3_end = simulate_model_given_noise_and_law(2*delta_t, T-delta_t, X0_3, params, delta_W_2, laws_2);
    X_4_end = simulate_model_given_noise_and_law(2*delta_t, T-delta_t, X0_4, params, delta_W_2, laws_2);
    steps_count_2 = (T-delta_t)/(2*delta_t);
    [X_1_end, X_2_end] = maximal_coupling_single_step_law(delta_t, X_1_end, X_2_end, params, laws_1(:,:,steps_count), laws_1(:,:,steps_count));
    [X_3_end, X_4_end] = maximal_coupling_single_step_law(2*delta_t, X_3_end, X_4_end, params, laws_2(:,:,steps_count_2), laws_2(:,:,steps_count_2));
end


function [X_1_end, X_2_end, X_3_end, X_4_end] = sample_K_l_bar(delta_t, X0_1, X0_2, X0_3, X0_4, params, alpha, laws_1, laws_2)
    a = rand();
    if a < alpha
        [X_1_end, X_2_end, X_3_end, X_4_end] = sample_Q_l_bar(delta_t, X0_1, X0_2, X0_3, X0_4, params, laws_1, laws_2);
        return;
    end
    [X_1_end, X_2_end, X_3_end, X_4_end] = sample_R_l_bar(delta_t, X0_1, X0_2, X0_3, X0_4, params, laws_1, laws_2);
end



% simulating functions
% X0 is d*particle_count
function X = simulate_model(delta_t, T, X0, params)
    particle_count = size(X0, 2);
    d = size(X0, 1);
    steps_count = round(T/delta_t);
    X = X0;
    delta_W = sqrt(delta_t)*randn(d, particle_count, steps_count);
    for i = 1:steps_count
        X = X + drift(params, X, X) * delta_t + ...
                          + diffusion(params, X, X, delta_t) * delta_W(:,:,i);
    end
end


function X = simulate_model_given_noise_and_law(delta_t, T, X0, params, delta_W, X_law)
    steps_count = round(T/delta_t);
    X = X0;
	delta_W_rep = zeros(size(delta_W,1),size(delta_W,1),size(delta_W,2), size(delta_W,3));
	for i = 1:size(delta_W,1)
		delta_W_rep(i,:,:,:) = delta_W;
	end
    for i = 1:steps_count
        X = X + drift(params, X, X_law(:,:,i)) * delta_t + ...
                          + squeeze(sum(diffusion(params, X, X_law(:,:,i), delta_t) .* delta_W_rep(:,:,:,i),2));
    end
end

function X = simulate_model_given_noise(delta_t, T, X0, params, delta_W)
    steps_count = round(T/delta_t);
    X = X0;
    for i = 1:steps_count
        X = X + drift(params, X, X) * delta_t + ...
                          + diffusion(params, X, X, delta_t) * delta_W(:,:,i);
    end
end

% simulate laws
function X = simulate_law_model(delta_t, T, X0, params)
    particle_count = size(X0, 2);
    d = size(X0, 1);
    steps_count = round(T/delta_t);
    X = zeros(d,particle_count,steps_count+1);
    X(:,:,1) = X0;
    delta_W = sqrt(delta_t)*randn(d, particle_count, steps_count);
	delta_W_rep = zeros(size(delta_W,1),size(delta_W,1),size(delta_W,2), size(delta_W,3));
	for i = 1:size(delta_W,1)
		delta_W_rep(i,:,:,:) = delta_W;
	end
    for i = 1:steps_count
        X(:,:,i+1) = X(:,:,i) + drift(params, X(:,:,i), X(:,:,i)) * delta_t + ...
                          + squeeze(sum(diffusion(params, X(:,:,i), X(:,:,i), delta_t) .* delta_W_rep(:,:,:,i),2));
    end
end

function X = simulate_law_model_given_noise(delta_t, T, X0, params, delta_W)
    steps_count = round(T/delta_t);
    particle_count = size(X0, 2);
    d = size(X0, 1);
    X = zeros(d,particle_count,steps_count+1);
    X(:,:,1) = X0;
	delta_W_rep = zeros(size(delta_W,1),size(delta_W,1),size(delta_W,2), size(delta_W,3));
	for i = 1:size(delta_W,1)
		delta_W_rep(i,:,:,:) = delta_W;
	end
    for i = 1:steps_count
        X(:,:,i+1) = X(:,:,i) + drift(params, X(:,:,i), X(:,:,i)) * delta_t + ...
                          + squeeze(sum(diffusion(params, X(:,:,i), X(:,:,i), delta_t) .* delta_W_rep(:,:,:,i),2));
    end
end




% X is d*particles_count
function A = drift(params, X, X_law)
	I = 0.5;
	J = 1;
	V_rev = 1;
	c = 0.08;
	b = 0.8;
	a = 0.7;
	T_max = 1;
	lambda = 0.2;
	a_d = 1;
	a_r = 1;
	V_T = 2;
    particle_count = size(X, 2);
    d = size(X, 1);
    A = zeros(d,particle_count);
    A(1,:) = X(1,:) - X(1,:).^3/3 - X(2,:) + I - J*(X(1,:)-V_rev)*mean(X_law(3,:));
    A(2,:) = c*(X(1,:)+a-b*X(2,:));
	A(3,:) = a_r*T_max*(1-X(3,:))./(1+exp(-lambda*(X(1,:)-V_T))) - a_d*X(3,:);
end

% ideally it should return d*d*particles_count, but ignored for now...
function A = diffusion(params, X, X_law, delta_t)
	I = 0.5;
	J = 1;
	V_rev = 1;
	c = 0.08;
	b = 0.8;
	a = 0.7;
	T_max = 1;
	lambda = 0.2;
	a_d = 1;
	a_r = 1;
	V_T = 2;
	b_ext = 0.5;
	b_J = 0.2;
	Gamma = 0.1;
	Lambda = 0.5;
    d = size(X, 1);
    particle_count = size(X, 2);
	A = zeros(d,d,particle_count);
    A(1,1,:) = b_ext;
	A(1,3,:) = -b_J*(X(1,:)-V_rev)*mean(X_law(3,:));
    X_3 = (X(3,:)>0) .* (X(3,:)<1) .* X(3,:);
    A(3,2,:) = (X(3,:)>0) .* (X(3,:)<1) .* Gamma .* exp(-Lambda./(1-(2*X_3-1).^2)) .* sqrt(a_r*T_max*(1-X_3)./(1+exp(-lambda*(X(1,:)-V_T))));
    A(2,3,:) = delta_t;
end

function X0 = sample_X0(params, number)
    X0 = zeros(3, number);
    X0(1,:) = params(1) + sqrt(params(4))*randn(1, number);
    X0(2,:) = params(2) + sqrt(params(5))*randn(1, number);
    X0(3,:) = params(3) + sqrt(params(6))*randn(1, number);
end

