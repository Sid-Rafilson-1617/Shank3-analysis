function [times, mean_frqs] = read_sniff_discrete(file_path, f, visualize, window_size, steps)
    % reads in NiDAQ sniff data file and finds inhalation times and uses simple moving average with window_size (seconds) to discretize sniff frequency. Arguements are file_path, sampling frequency, and
    % visualization of processing

    % degree of polynomial used for savitsky golay smoothing
    degree = 25;

    % read in file
    path_back = pwd;
    cd(file_path)
    sniffpile = dir('NiDAQ_sniff.dat');
    cd(sniffpile.folder)
    sniff_file = fopen(sniffpile.name,'r');
    sniff=fread(sniff_file,'float64');
    sniff=sniff(1:4:length(sniff)); %analog sniff trace is every fourth entry
    cd(path_back)
    % smooth signal and asign time in seconds
    sniff = smooth(sniff, degree,'sgolay');
    times = (1:length(sniff))./f;
    
    % plot smoothed signal
    if visualize == true
        figure;
        plot(times, sniff,'k-','LineWidth',1)
        ylabel('Temperature (AU)')
        xlabel('Time (s)')
        title('sgolay smoothed signal');
        hold off
    end

    %peak detection 
    fprintf('begining peak detection for %s...\n', file_path)
    inhale_times = [];
    if ~isempty(sniff)
        windows = round(1:(f*window_size):length(sniff));
        zniff = zeros(1, length(windows));
        for scan = 2:length(windows)
            time_stamp=(windows(scan-1)):windows(scan);
            zniff=zscore(sniff( time_stamp));

            [sniff_pks,in_locs] = findpeaks(zniff, ...
                    'MinPeakDistance',degree,'MinPeakProminence',0.5);

                %visual validation of peak finder
                if (scan == 5) && (visualize == true)
                    figure;
                    hold on
                    plot(zniff)
                    plot(sniff(time_stamp), "Color" , "Black")
                    plot(in_locs,zniff(in_locs),'ro')
                    plot(ex_locs,zniff(ex_locs),'go')
                    title('Peak Finder Validation')
                    hold off
                end

                for peak= 1:length(in_locs)
                    inhale=in_locs(peak);
                    if ~isempty(inhale)
                        inhale_times=[inhale_times; inhale + windows(scan-1)];
                        sniff_frequencies = f./diff(inhale_times);
                        inhalation_times_seconds = inhale_times(2:end)./f;                       
                    end
                end
        end
    end


    %smoothing and discretizing the signal by averaging the instantaneous sniff
    %frequency in sliding bins of size window_size/steps
    num_chunks = ceil(inhalation_times_seconds(end) / (window_size));
    mean_frqs = zeros(1, num_chunks*steps);

    width = window_size/steps;
    count = 0;
    for chunk = 1:num_chunks
        for window = 1: steps
            count = count + 1;
            lower_bound = (chunk + ((window - 1 )/steps) - 1) * width;
            upper_bound =(chunk + ((window - 1)/steps)) * width;
            indices_in_window = find(inhalation_times_seconds >= lower_bound &...
            inhalation_times_seconds <= upper_bound);
            if ~isempty(indices_in_window)
                mean_frqs(count) = mean(sniff_frequencies(indices_in_window(1):indices_in_window(end)));
            else
                mean_frqs(count) = NaN;
                fprintf('empty window durring chunk %g \n', count)
            end
        end
    end

    mean_frqs = mean_frqs(1:count);
    mean_frqs = mean_frqs';
    fprintf('Peak detection complete\n\n');
    times = linspace(1, length(mean_frqs), length(mean_frqs));
    times = times';
end

    
