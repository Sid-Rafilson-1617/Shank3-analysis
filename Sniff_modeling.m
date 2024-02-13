%% ________________________ CONFIG _____________________________


main_config = main_config();

configToWorkspace(main_config)


%% ____________________________ PRE-PROCESSING ____________________________

%% reading signal for continuous time instantaneous sniff signal

data = cell(Nmice,Nsignals);
raw_sniff = cell(Nmice, Nsignals);
for i = 1:Nmice

    %data path directory
    if i == 1
        number = "3001";
    elseif i == 2
        number = "3005";
    elseif i == 3
        number = "3005";
    elseif i == 4
        number = "3006";
    end

    path = ("E:\Shank3\data\" + number + "\Spontaneous");
    % reading in the data. Propogates a 1:N cell array for N signals.
    
    for session = 1:Nsignals
        file_path = path + '\'+ num2str(session);
        [time, freqs, sniff] = read_sniff_continuous(file_path, f, visualization); 
        data{i, session} = [time, freqs];
        raw_sniff{i, session} = sniff;
    end
end

%% reading signal for discrete time smoothed  sniff signal

%window size in seconds
window_size = 6;

%number of bins per window
steps = 6;

data = cell(Nmice,Nsignals);
for i = 1:Nmice

    %data path directory
    if i == 1
        number = "3001";
    elseif i == 2
        number = "3002";
    elseif i == 3
        number = "3005";
    elseif i == 4
        number = "3006";
    end

    path = ("E:\Shank3\data\" + number + "\Spontaneous");
    % reading in the data. Propogates a 1:N cell array for N signals.
    
    for session = 1:Nsignals
        file_path = path + '\'+ num2str(session);
        [time, freqs] = read_sniff_discrete(file_path, f, visualization, window_size, steps); 
        data{i, session} = [time, freqs];
        %time is in seconds * window_size(seconds)
    end
end


%% Visualize raw sniff accross time
for i = 1:Nmice
    % Visualizing the sniff frequency time series
    figure;
    hold on
    for session = 1:Nsignals
        plot(raw_sniff{i, session})
    end
    title(mouse_names{i})
    xlabel('Time (s)')
    ylabel('Voltage')
    hold off
end


%% Visualize sniff freqs accross time

for mouse = 1:Nmice
    for session = 1:Nsignals
        figure;
        hold on;
        for bin = 1:Nbins
            x = data_bins{mouse, session, bin}(:,1); % x-values
            y = data_bins{mouse, session, bin}(:,2); % y-values
            scatter(x, y, [], y, '.'); % Color-code by y-value
        end
        colormap('winter');
        hold off;
        set(gca, 'YScale', 'log');
        yticks([2 4 8 14]);
        yticklabels({'2', '4', '8', '14'});
        yline(4);
        yline(8);
        title(mouse_names{mouse});
        subtitle(signal_names{session});
        xlabel('Time (s)');
        ylabel('Sniff Frequency (Hz)');
    end
end

%% subdiving signals based on frequency bin membership

bins = {[0,4],[4,8],[8,inf]};
Nbins = length(bins);
data_bins = cell(Nmice,Nsignals,Nbins);
for mouse = 1:Nmice
    for session = 1:Nsignals
        for bin = 1:Nbins
            indices_in_bin = find(data{mouse, session}(:,2) > bins{bin}(1) & data{mouse, session}(:,2) < bins{bin}(2));
            extracted_data = data{mouse, session}(indices_in_bin, :);
            data_bins{mouse, session, bin} = extracted_data;
        end
    end
end

%% Plotting result

colors = copper(Nbins);

for mouse = 1:Nmice
    for session = 1:Nsignals
        figure;
        hold on
        for bin = 1:Nbins
            scatter(data_bins{mouse, session, bin}(:,1), data_bins{mouse, session, bin}(:,2), '.','CData', colors(bin, :));
        end
        hold off
        set(gca, 'YScale', 'log')
        yticks([2 4 8 14])
        yticklabels([2 4 8 14])
        title(mouse_names{mouse})
        subtitle(signal_names{session})
        xlabel('Time (s)')
        ylabel('Sniff Frequency (Hz)')
    end
end


%% _______________________ANALYSIS___________________


sniff_config = sniff_analysis_config();

configToWorkspace(sniff_config)

%% Estimating the Autocovariance function

%maximum lag (seconds)
max_lag = 5;

ACVs = cell(Nmice, Nsignals);
lags = cell(Nmice, Nsignals);

for mouse = mouse_0:mouse_1
    for trace = signal_0:signal_1
        [ACVs{mouse,trace}, lags{mouse,trace}] = ...
            xcov(data{mouse, trace}(:,2), floor(max_lag * f), 'none');
        lags{mouse,trace} = lags{mouse,trace}/f;
    end
end


%% Plotting the ACV

for mouse = mouse_0:mouse_1
    figure;
    t = tiledlayout(signal_1 - signal_0 + 1, 1);
    for trace = signal_0:signal_1
        nexttile;
        plot(lags{mouse,trace}, ACVs{mouse,trace});
        title(sprintf('%s',...
            signal_names{trace}))
        xlabel('Time Lag (s)')
        ylabel('ACV')
    end
    sgtitle(mouse_names{mouse})
end

%% AutoCorrelation

lag0 = zeros(1,1,1,1);
ACRs = cell(Nmice, Nsignals, Nmice, Nsignals);
for mouse = mouse_0:mouse_1
    for trace = signal_0:signal_1
                lag0(mouse, trace) = ...
                    find(lags{mouse, trace} == 0);

                ACRs{mouse, trace} = ...
                    ACVs{mouse, trace} ...
                    / ACVs{mouse, trace}(lag0(mouse, trace));
    end
end

% plotting
for mouse = 1:Nmice
    figure;
    t = tiledlayout(Nsignals, 1);
    for trace = signal_0:signal_1
            nexttile;
            plot(lags{mouse,trace}, ACRs{mouse,trace});
            yline(mean(ACRs{mouse,trace}), 'r', 'LineWidth', 2)
            title(sprintf('%s',...
                signal_names{trace}))
            xlabel('Time Lag (s)')
            ylabel('Correlation')
            legend('Mean Correlation')
    end
    sgtitle(mouse_names{mouse})
end

%% Power Spectral Density

noverlap = 10*f;
L = 20*f;
PSD = cell(Nmice, Nsignals, Nmice, Nsignals);
angular_frequencies = cell(Nmice, Nsignals, Nmice, Nsignals);


for mouse = mouse_0:mouse_1
    for trace = signal_0:signal_1
        for cross_mouse = mouse_0:mouse_1
            for cross_trace = signal_0:signal_1

                signal1 = raw_sniff{mouse, trace};
                signal2 = raw_sniff{cross_mouse, cross_trace};
                
                % Check the lengths of the signals and zero-pad the shorter signal if necessary
                len1 = length(signal1);
                len2 = length(signal2);
                if len1 < len2
                    signal1 = [signal1; zeros(len2-len1, 1)];
                elseif len2 < len1
                    signal2 = [signal2; zeros(len1-len2, 1)];
                end

                [PSD{mouse, trace, cross_mouse, cross_trace},...
                angular_frequencies{mouse, trace, cross_mouse, cross_trace}] =...
                cpsd(signal1, signal2, hamming(L), noverlap, [], f);
   
            end
        end
    end
end

for mouse = mouse_0:mouse_1
    figure;
    t = tiledlayout(Nsignals, Nsignals);
    for trace = signal_0:signal_1
        for cross_trace = signal_0:signal_1
            nexttile;
            hold on
            plot(angular_frequencies{mouse, trace, mouse, cross_trace}/(2*pi), log(abs(PSD{mouse, trace, cross_mouse, cross_trace})))
            title(sprintf('%s vs %s (S_{%d,%d})',...
                signal_names{trace},signal_names{cross_trace},trace,cross_trace))
            ylabel("log Power Spectrum")
            xlabel('Frequency')
            xlim([0 14])
            hold off
        
        end
    end
    PSDtitle = {mouse_names{mouse}, 'instantaneous Frequencies', "Welch's Method" , '20s windows'};
    sgtitle(PSDtitle)
end

%% __________________________________Modeling__________________________________


%% fit linear regression model
for i = 1:Nmice
    figure;
    hold on
    for session = 1:Nsignals
        lin_reg = fitlm(data{i, session}(:,1), data{i, session}(:,2));
        scatter(data{i, session}(:,1), data{i, session}(:,2))
        plot(lin_reg)
    end
    xlabel('Time (s)')
    ylabel('Sniff Frequency (Hz)')
    hold off
end


%% detrend signal (for stationarity assumption)
data_detrended = cell(Nmice,Nsignals);
for i = 1:Nmice
    for session = 1:Nsignals
        signal_detrended = detrend(data{i, session}(:,2));
        data_detrended{i, session}(:,1) = data{i,session}(:,1);
        data_detrended{i, session}(:,2) = signal_detrended + mean(data{i, session}(:,2));
    end
end



%% Plotting detrended signal

for i = 1:Nmice
    % Visualizing the sniff frequency time series
    figure;
    hold on
    for session = 1:Nsignals
        scatter(data_detrended{i, session}(:,1), data_detrended{i, session}(:,2), '.')
        plot(data{i, session}(:,1), data{i, session}(:,2) - data_detrended{i,session}(:,2) + mean(data{i, session}(:,2)))
    end
    xlabel('Time (s)')
    ylabel('Sniff Frequency (Hz)')
    hold off
end






%% training a CTMC for all signals in parallel

CTMC_config = CTMC_config(main_config);

for i = 1:Nmice
    [best_fitted_Q, empirical_proportions, sigma] = fit_CTMC_sniff(NS, data(i,:), Ntrials, 12, true);
    Q{i} = best_fitted_Q;
    empirical{i} = empirical_proportions;
    stationary{i} = sigma;
end




%% training a DTMC for all signals in parallel
[best_fitted_P, empirical_proportions, sigma] = fit_DTMC_sniff(NS, min_frequency, max_frequency, data, Ntrials, 12, false);
P{i} = best_fitted_P;
empirical{i} = empirical_proportions;
stationary{i} = sigma;




%% visual for transition matrix CTMC
figure;
heat_matrix = Q{1} - diag(Q{1});
heatmap(heat_matrix, 'Colormap', copper);




%% visual for transition matrix DTMC
figure;
heatmap(P{1}, 'Colormap', copper);




%% simulating a trajectory CTMC
max_time = 1;
[transition_times, states] = simulate_CTMC(best_fitted_Q, max_time);

figure;
plot(transition_times, states);



%% Graphical Analysis for Stationary Distribution

mouse = 1;


total_proportion = mean(empirical{mouse}, 2);

figure;
hold on
sigma_bar = bar(stationary{mouse});

colors = summer(size(empirical{mouse},2));
empirical_bar = bar(empirical{mouse});

for i = 1:size(empirical{mouse}, 2)
    set(empirical_bar(i), 'FaceColor', colors(i, :));
end

ylabel('State Proportion')
xlabel('Sniff frequency bins')
title('CTMC stationary distribution')
for i = 1:size(empirical{mouse},2)
    h(i) = patch(NaN, NaN, colors(i), "Visible", "off");
end

legend([sigma_bar, h], {'stationary distribution', 'empirical proportions'});
hold off




%%
series = data{1,5}(:,1);
save sniff_frequencies series

%%
load Q_matrix best_fitted_Q



