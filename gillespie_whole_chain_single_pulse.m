clear all;

%%% Parameters %%%%%%%%
%num_enhancers = 55;
%num_molecules = num_enhancers + 2;

%slopes = [0.95, 0.99];
f = 1.0;
RC = 1.0;
%r = 0.8;
%g = 0.8;
%ma = 1.05;
%mr = 100.5;
for num_enhancers = 20
    num_molecules = num_enhancers + 2;
    %slopes_1_values = [0.89,0.9,0.95,0.99,0.999,1.001,1.01,1.05,1.1,1.5];%0.99;
    slopes_1_values = 0.95;
    slopes_2_values = 1.01;
    %slopes_2_values = [0.89,0.9,0.95,0.99,0.999,1.001,1.01,1.05,1.1,1.5];
    
    for s1 = slopes_1_values
        for s2 = slopes_2_values
            slopes(1) = s1;
            slopes(2) = s2;
            r_values = 5.0;
            g_values = [1.0, 2.0, 3.0, 4.0, 5.0];
            for r = r_values
                for g = g_values
                    ma_values = [1.05,5.5,10.5,100.5];
                    for ma = ma_values
                        mr_values = [1.05,5.5,10.5,100.5];
                        for mr = mr_values
                            for s = 1:40
                                runs = 100000;

                                %%% Defining state %%%%
                                enh = zeros(3, num_molecules);
                                rate = cell(3, num_molecules);

                                for i = 1:size(rate, 1)
                                    for j = 2:size(rate, 2) - 1
                                        rate{1, j} = zeros(2); % Initialize rate matrix for activator
                                        rate{2, j} = zeros(2); % Initialize rate matrix for enhancer
                                        rate{3, j} = zeros(2); % Initialize rate matrix for repressor
                                    end
                                end

                                enh(2, 2) = 1; % Initial state
                                time = 0;
                                % fileID1 = fopen('test1.txt', 'w');
                                fileID2 = fopen(sprintf('enh_num_%d_r_%.1f_g_%.1f_s1_%.3f_s2_%.3f_ma_%.1f_mr_%.1f_AC_100_t%d.txt', num_enhancers, r, g, slopes(1), slopes(2), ma, mr, s), 'w');
                                %fileID2 = fopen(sprintf('enh_num_%d_s1_%d_s2_%d_ma_%d_mr_%d_AC_100_t%d.txt', num_enhancers,slopes(1),slopes(2),ma,mr,s), 'w');
                                avg_col = 0.0;

                                for i = 1:runs
                                    % Generate new rate values
                                    AC = rect(i); % This could be updated based on your specific logic
                                    [a, b, c, d, u, v, e_to_ea, ea_to_e, e_to_er, er_to_e] = generate_values(slopes, num_molecules, r, g, ma, mr, AC, RC, f);

                                    % Update rate matrices
                                    for j = 2:size(rate, 2) - 1
                                        rate{1, j}(1, 1) = 0; % up rate
                                        rate{1, j}(1, 2) = ea_to_e(j); % down rate
                                        rate{1, j}(2, 1) = c(j); % right rate
                                        rate{1, j}(2, 2) = d(j - 1); % left rate

                                        rate{2, j}(1, 1) = e_to_ea(j); % up rate
                                        rate{2, j}(1, 2) = e_to_er(j); % down rate
                                        rate{2, j}(2, 1) = a(j); % right rate
                                        rate{2, j}(2, 2) = b(j - 1); % left rate

                                        rate{3, j}(1, 1) = er_to_e(j); % up rate
                                        rate{3, j}(1, 2) = 0; % down rate
                                        rate{3, j}(2, 1) = u(j); % right rate
                                        rate{3, j}(2, 2) = v(j - 1); % left rate
                                    end

                                    % Simulation
                                    [row, col] = find(enh == 1);
                                    index_i = row;
                                    index_j = col;

                                    avg_col = avg_col + col;

                                    c1 = rate{index_i, index_j}(1, 1); % up
                                    c2 = rate{index_i, index_j}(1, 2); % down
                                    c3 = rate{index_i, index_j}(2, 1); % right
                                    c4 = rate{index_i, index_j}(2, 2); % left

                                    R = c1 + c2 + c3 + c4; % sum rates of all the reactions
                                    r2 = R * rand(); % generate random number between 0 and R
                                    r1 = rand(); % another random number for time
                                    t = (1 / R) * log(1 / r1); % time from exponential distribution
                                    time = time + t; % time at which next event happens

                                    % Update enhancer states based on reaction rates
                                    if index_i > 0 && r2 > 0 && r2 <= c1
                                        enh(index_i, index_j) = 0;
                                        enh(index_i - 1, index_j) = 1;
                                    elseif index_i < 3 && r2 > c1 && r2 <= c1 + c2
                                        enh(index_i, index_j) = 0;
                                        enh(index_i + 1, index_j) = 1;
                                    elseif index_j < num_molecules && r2 > c1 + c2 && r2 <= c1 + c2 + c3
                                        enh(index_i, index_j) = 0;
                                        enh(index_i, index_j + 1) = 1;
                                    elseif index_j > 1 && r2 > c1 + c2 + c3 && r2 <= c1 + c2 + c3 + c4
                                        enh(index_i, index_j) = 0;
                                        enh(index_i, index_j - 1) = 1;
                                    end

                                    % Logging the state
                                    % fprintf(fileID1, 'Step %f:\n', time);
                                    % for row = 1:size(enh, 1)
                                    %     fprintf(fileID1, '%f ', enh(row, :));
                                    %     fprintf(fileID1, '\n');
                                    % end
                                    % fprintf(fileID1, '\n');

                                    fprintf(fileID2, '%f %f %d\n', time, AC, col - 1);

                                    %disp(AC);
                                    %disp(avg_col);
                                end

                                cl = avg_col / runs;
                                %disp(cl);

                                % fclose(fileID1);
                                fclose(fileID2);
                            end
                        end
                    end
                end
            end
        end
    end
end

%%% Function to generate the rates
function [a, b, c, d, u, v, e_to_ea, ea_to_e, e_to_er, er_to_e] = generate_values(slopes, num_molecules, r, g, ma, mr, AC, RC, f)
    % Initialize arrays for each variable with initial values
    a = zeros(1, num_molecules);
    b = zeros(1, num_molecules);
    c = zeros(1, num_molecules);
    d = zeros(1, num_molecules);
    u = zeros(1, num_molecules);
    v = zeros(1, num_molecules);
    e_to_ea = zeros(1, num_molecules);
    ea_to_e = zeros(1, num_molecules);
    e_to_er = zeros(1, num_molecules);
    er_to_e = zeros(1, num_molecules);

    for i = 2:(num_molecules - 2)
        a(i) = f * max(0, r * (slopes(1) ^ ((i - 1))));
        b(i) = f * max(0, g * (slopes(2) ^ ((i - 1))));
    end

    for i = 2:(num_molecules - 2)
        c(i) = f * max(0, (ma * r) * (slopes(1) ^ ((i - 1))));
        d(i) = f * max(0, g * (slopes(2) ^ (i - 1)));
        u(i) = f * max(0, (r) * (slopes(1) ^ ((i - 1))));
        v(i) = f * max(0, (mr * g) * (slopes(2) ^ ((i - 1))));
    end

    for i = 2:(num_molecules - 1)
        e_to_ea(i) = f * AC * r;
        ea_to_e(i) = f * g;
        e_to_er(i) = f * RC * r;
        er_to_e(i) = f * g;
    end
end

function AC = rect(i)
    T = 100000; % Cycle duration
    amplitude = 100; % Amplitude of the rectangular pulse
    AC = amplitude * (mod(i, T) > T / 2);
end
