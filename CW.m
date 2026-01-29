
% job_id = getenv('SLURM_JOB_ID');
% proc_id = getenv('SLURM_PROCID');
job_id = '1';
proc_id = '1';
folder_read = '';
folder_write = sprintf('%s%s', job_id, '/');
results_filename = sprintf('%ssamples_%s_%s.txt', folder_write, job_id, proc_id);
results_filename_3 = sprintf('%sparticles_1_%s_%s.txt', folder_write, job_id, proc_id);
results_filename_4 = sprintf('%sparticles_2_%s_%s.txt', folder_write, job_id, proc_id);
results_filename_6 = sprintf('%slevels_%s_%s.txt', folder_write, job_id, proc_id);
results_filename_7 = sprintf('%sps_%s_%s.txt', folder_write, job_id, proc_id);
results_filename_8 = sprintf('%stimes_%s_%s.txt', folder_write, job_id, proc_id);


T = 100;
X0 = 1;
theta = 0;    %useless
sigma = 2;    %useless
alpha = 0.9;
params = [theta, sigma];   %useless
phi = @(x) x.^2;

l_min = 3;
l_max = 10;
p_min = 0;
p_max = 7;
sample_count = 100;
samples = zeros(1,sample_count);
ls = zeros(1,sample_count);
ps = zeros(1,sample_count);
particles_1 = [];
particles_2 = [];
for i = 1:sample_count
    disp(i);
    [level, l_density] = sample_l(l_max, l_min);
    [p, p_density] = sample_p(p_max, p_min);
    ls(i) = level;
    ps(i) = p;
    delta_t = 2^(-level);
    disc_count = 2^(level);
    law_particles_count = round(10*level);
    I_p = 10*2^p;
    I_p_minus_1 = 10*2^(p-1);
    T_initial = 0;

    %disp(['level = ', num2str(level)]);
    iter_time_start = tic;
    if level == l_min
        X_law = simulate_model(delta_t, I_p, X0*ones(law_particles_count,1), params);
        delta_W = sqrt(delta_t)*randn(1, I_p*disc_count);
        X = simulate_model_given_noise_and_law(delta_t, I_p, X0, params, delta_W, X_law);
        unit_indices = 1:disc_count:I_p*disc_count+1;
        X_unit = X(unit_indices);
        if p ==0
            samples(i) = 1/(l_density*p_density)*mean(phi(X_unit));
        else
            samples(i) = 1/(l_density*p_density)*(mean(phi(X_unit)) - mean(phi(X_unit(1:I_p_minus_1))));
        end
        particles_1 = [particles_1, X_unit];        
    end

    if level > l_min
        delta_W_1 = sqrt(delta_t)*randn(law_particles_count, I_p/delta_t);
        delta_W_2 = delta_W_1(:,1:2:end) + delta_W_1(:,2:2:end);
        X_law_1 = simulate_model_given_noise(delta_t, I_p, X0*ones(law_particles_count,1), params, delta_W_1);
        X_law_2 = simulate_model_given_noise(2*delta_t, I_p, X0*ones(law_particles_count,1), params, delta_W_2);
        delta_W_1 = sqrt(delta_t)*randn(1, I_p*disc_count);
        delta_W_2 = delta_W_1(:,1:2:end) + delta_W_1(:,2:2:end);
        X_1 = simulate_model_given_noise_and_law(delta_t, I_p, X0, params, delta_W_1, X_law_1);
        X_2 = simulate_model_given_noise_and_law(2*delta_t, I_p, X0, params, delta_W_2, X_law_2);
        unit_indices = 1:disc_count:I_p*disc_count+1;
        X_1_unit = X_1(unit_indices);
        unit_indices = 1:disc_count/2:I_p*disc_count/2+1;
        X_2_unit = X_2(unit_indices);
        if p ==0
            samples(i) = 1/(l_density*p_density)*(mean(phi(X_1_unit))-mean(phi(X_2_unit)));
        else
            samples(i) = 1/(l_density*p_density)*(mean(phi(X_1_unit)) - mean(phi(X_1_unit(1:I_p_minus_1))) - mean(phi(X_2_unit)) + mean(phi(X_2_unit(1:I_p_minus_1))));
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
    density_unnorm = 2.^(-levels) .* (levels+1) .*log2(levels+2).^2;
    density = density_unnorm/sum(density_unnorm);
    p = randsample(levels,1,true,density);
    p_density = density(p-p_min+1);
end


% the method functions

function [X_1_end, X_2_end] = sample_Q_l_caron(delta_t, X0_1, X0_2, params, laws)
    T = 1;
    particles_count = 1;
    steps_count = T/delta_t;
    delta_W = sqrt(delta_t)*randn(particles_count, steps_count);
    X_1 = simulate_model_given_noise_and_law(delta_t, T, X0_1, params, delta_W, laws);
    X_2 = simulate_model_given_noise_and_law(delta_t, T, X0_2, params, delta_W, laws);
    X_1_end = X_1(:,end);
    X_2_end = X_2(:,end);
end

function [X_1_end, X_2_end] = sample_R_l_caron(delta_t, X0_1, X0_2, params, laws)
    T = 1-delta_t;
    %particles_count = length(X0_1);
    particles_count = 1;
    steps_count = T/delta_t;
    delta_W = sqrt(delta_t)*randn(particles_count, steps_count);
    %X_laws = simulate_discrete_kuramoto(delta_t, T, intial_law, theta, sigma);
    X_1 = simulate_model_given_noise_and_law(delta_t, T, X0_1, params, delta_W, laws);
    X_2 = simulate_model_given_noise_and_law(delta_t, T, X0_2, params, delta_W, laws);
    X_1_end = X_1(:,end);
    X_2_end = X_2(:,end);
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
    X_1 = sqrt(var_1)*randn() + mean_1;
    W = rand() * normpdf(X_1, mean_1, sqrt(var_1));
    if W < normpdf(X_1, mean_2, sqrt(var_2))
        X_2 = X_1;
        return;
    end
    while true
        X_2 = sqrt(var_2)*randn() + mean_2;
        W = rand() * normpdf(X_2, mean_2, sqrt(var_2));
        if W > normpdf(X_2, mean_1, sqrt(var_1))
            return;
        end
    end
end

function [means_1, means_2, var_1, var_2] = calculate_mean_var_single_level_law(delta_t, x_1, x_2, params, X_law_1, X_law_2)
    means_1 = x_1 + drift(params, x_1, X_law_1) * delta_t;
    means_2 = x_2 + drift(params, x_2, X_law_2) * delta_t;
    var_1 = diffusion(params, x_1, X_law_1)^2 * delta_t;
    var_2 = diffusion(params, x_2, X_law_2)^2*delta_t;
end

function [unit_sequence_1, unit_sequence_2, meeting_time] = sample_coupled_single_level_till_meeting_given_laws(T, delta_t, X0, params, alpha, X_law)
    meeting_time = T+1;
    law_particles_count = size(X_law,1);
    particles_count = 1;
    unit_sequence_1 = zeros(particles_count, T+1);
    unit_sequence_2 = zeros(particles_count, T+1);
    delta_W = sqrt(delta_t)*randn(1, 1/delta_t);
    %X_1s = simulate_model(delta_t, 1, X0*ones(law_particles_count,1), params);
    X_1s = simulate_model_given_noise_and_law(delta_t, 1, X0, params, delta_W, X_law(:,1:1/delta_t));
    unit_sequence_1(:,1) = X_1s(1,end);
    unit_sequence_2(:,1) = X0;
    for i = 1:T
        [unit_sequence_1(:,i+1), unit_sequence_2(:,i+1)] = sample_K_l_caron(delta_t, unit_sequence_1(:,i), unit_sequence_2(:,i), params, alpha, X_law(:,(i/delta_t+1):((i+1)/delta_t)));
        %if unit_sequence_1(:,i+1) == unit_sequence_2(:,i)
        if equals_floats(unit_sequence_1(:,i+1), unit_sequence_2(:,i+1))
            meeting_time = i+1;
            break;
        end
    end
end

function b = equals_floats(x,y)
    if abs(x-y) < 1e-20
        b = true;
        return;
    end
    b = false;
end
 


% coulped
function [X_1_end, X_2_end, X_3_end, X_4_end] = sample_Q_l_bar(delta_t, X0_1, X0_2, X0_3, X0_4, params, laws_1, laws_2)
    T = 1;
    %particles_count = length(X0_1);
    particles_count = 1;
    steps_count = T/delta_t;
    delta_W = sqrt(delta_t)*randn(particles_count, steps_count);
    X_1 = simulate_model_given_noise_and_law(delta_t, T, X0_1, params, delta_W, laws_1);
    X_2 = simulate_model_given_noise_and_law(delta_t, T, X0_2, params, delta_W, laws_1);
    X_1_end = X_1(:,end);
    X_2_end = X_2(:,end);
    delta_W_2 = delta_W(:,1:2:end) + delta_W(:,2:2:end);
    X_3 = simulate_model_given_noise_and_law(2*delta_t, T, X0_3, params, delta_W_2, laws_2);
    X_4 = simulate_model_given_noise_and_law(2*delta_t, T, X0_4, params, delta_W_2, laws_2);
    X_3_end = X_3(:,end);
    X_4_end = X_4(:,end);
end

function [X_1_end, X_2_end, X_3_end, X_4_end] = sample_R_l_bar(delta_t, X0_1, X0_2, X0_3, X0_4, params, laws_1, laws_2)
    T = 1-delta_t;
    %particles_count = length(X0_1);
    particles_count = 1;
    steps_count = T/delta_t;
    delta_W = sqrt(delta_t)*randn(particles_count, steps_count);
    X_1 = simulate_model_given_noise_and_law(delta_t, T, X0_1, params, delta_W, laws_1);
    X_2 = simulate_model_given_noise_and_law(delta_t, T, X0_2, params, delta_W, laws_1);
    X_1_end = X_1(:,end);
    X_2_end = X_2(:,end);
    delta_W_2 = delta_W(:,1:2:end-1) + delta_W(:,2:2:end-1);
    X_3 = simulate_model_given_noise_and_law(2*delta_t, T-delta_t, X0_3, params, delta_W_2, laws_2);
    X_4 = simulate_model_given_noise_and_law(2*delta_t, T-delta_t, X0_4, params, delta_W_2, laws_2);
    X_3_end = X_3(:,end);
    X_4_end = X_4(:,end);
    steps_count_2 = (T-delta_t)/(2*delta_t);
    [X_1_end, X_2_end] = maximal_coupling_single_step_law(delta_t, X_1_end, X_2_end, params, laws_1(:,steps_count), laws_1(:,steps_count));
    [X_3_end, X_4_end] = maximal_coupling_single_step_law(2*delta_t, X_3_end, X_4_end, params, laws_2(:,steps_count_2), laws_2(:,steps_count_2));
end


function [X_1_end, X_2_end, X_3_end, X_4_end] = sample_K_l_bar(delta_t, X0_1, X0_2, X0_3, X0_4, params, alpha, laws_1, laws_2)
    a = rand();
    if a < alpha
        [X_1_end, X_2_end, X_3_end, X_4_end] = sample_Q_l_bar(delta_t, X0_1, X0_2, X0_3, X0_4, params, laws_1, laws_2);
        return;
    end
    [X_1_end, X_2_end, X_3_end, X_4_end] = sample_R_l_bar(delta_t, X0_1, X0_2, X0_3, X0_4, params, laws_1, laws_2);
end

function [unit_sequence_1, unit_sequence_2, unit_sequence_3, unit_sequence_4, meeting_time_1, meeting_time_2] = sample_coupled_coupled_level_till_meeting_given_laws(T, delta_t, X0_1, X0_2, params, alpha, X_law_1, X_law_2)
    meeting_time_1 = T+1;
    meeting_time_2 = T+1;
    law_particles_count = size(X_law_1,1);
    particles_count = 1;
    unit_sequence_1 = zeros(particles_count, T+1);
    unit_sequence_2 = zeros(particles_count, T+1);
    unit_sequence_3 = zeros(particles_count, T+1);
    unit_sequence_4 = zeros(particles_count, T+1);
    delta_W_1 = sqrt(delta_t)*randn(1, 1/delta_t);
    delta_W_2 = delta_W_1(:,1:2:end) + delta_W_1(:,2:2:end);
    X_1s = simulate_model_given_noise_and_law(delta_t, 1, X0_1, params, delta_W_1, X_law_1(:,1:1/delta_t));
    X_3s = simulate_model_given_noise_and_law(2*delta_t, 1, X0_2, params, delta_W_2, X_law_2(:,1:1/(2*delta_t)));
    unit_sequence_1(:,1) = X_1s(1,end);
    unit_sequence_2(:,1) = X0_1;
    unit_sequence_3(:,1) = X_3s(1,end);
    unit_sequence_4(:,1) = X0_2;
    break_1 = false;
    break_2 = false;
    for i = 1:T
        [unit_sequence_1(:,i+1), unit_sequence_2(:,i+1), unit_sequence_3(:,i+1), unit_sequence_4(:,i+1)] = sample_K_l_bar(delta_t, ...
                    unit_sequence_1(:,i), unit_sequence_2(:,i), unit_sequence_3(:,i), unit_sequence_4(:,i), params, alpha, ...
                    X_law_1(:,(i/delta_t+1):((i+1)/delta_t)), X_law_2(:,(i/(2*delta_t)+1):((i+1)/(2*delta_t))));
        %if unit_sequence_1(:,i+1) == unit_sequence_2(:,i)
        if equals_floats(unit_sequence_1(:,i+1), unit_sequence_2(:,i+1)) && ~break_1
            meeting_time_1 = i+1;
            break_1 = true;
        end
        if equals_floats(unit_sequence_3(:,i+1), unit_sequence_4(:,i+1)) && ~break_2
            meeting_time_2 = i+1;
            break_2 = true;
        end
        if break_1 && break_2
            break;
        end
    end
end



% simulating functions
function X = simulate_model(delta_t, T, X0, params)
    particle_count = length(X0);     % is this good?
    steps_count = round(T/delta_t);
    X = zeros(particle_count, steps_count+1);
    X(:,1) = X0;
    delta_W = sqrt(delta_t)*randn(particle_count, steps_count);
    for i = 1:steps_count
        X(:,i+1) = X(:,i) + drift(params, X(:,i), X(:,i)) * delta_t + ...
                          + diffusion(params, X(:,i), X(:,i)) .* delta_W(:,i);
    end
end


function X = simulate_model_given_noise_and_law(delta_t, T, X0, params, delta_W, X_law)
    particle_count = length(X0);     % is this good?
    steps_count = round(T/delta_t);
    X = zeros(particle_count, steps_count+1);
    X(:,1) = X0;
    for i = 1:steps_count
        X(:,i+1) = X(:,i) + drift(params, X(:,i), X_law(:,i)) * delta_t + ...
                          + diffusion(params, X(:,i), X_law(:,i)) .* delta_W(:,i);
    end
end

function X = simulate_model_given_noise(delta_t, T, X0, params, delta_W)
    particle_count = length(X0);     % is this good?
    steps_count = round(T/delta_t);
    X = zeros(particle_count, steps_count+1);
    X(:,1) = X0;
    for i = 1:steps_count
        X(:,i+1) = X(:,i) + drift(params, X(:,i), X(:,i)) * delta_t + ...
                          + diffusion(params, X(:,i), X(:,i)) .* delta_W(:,i);
    end
end


function A = drift(params, X, X_law)
    K = 0.25;
    A = -X.^3 + X + K*mean(X_law);
end

function A = diffusion(params, X, X_law)
    A = 1;
end








