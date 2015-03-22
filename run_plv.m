function PLV = run_plv(corrcoeff, filter, noise, phase_method, freq, do_plot)

% in this variable corrcoeff or Noiseampli can't be provided lower and
% upper range at a time.
% provision two add noise and filter the desired signal has been kept in
% this script. for this 'noiseadd' is selected to add noise in the
% desired signal. for example noiseadd=[0,1] will add noise to 2nd signal
% only and similarly [0,1]= for 2nd signal noise, [1,0]= for 1st signal
% noise, [1,1]= noise in both signals and [0,0] = noise in none
% similar pattern for sleection of signal for filteration is opted in 'filtersignal'

%example    run_plv([0.1 0], [1 3],[1,0],[1,1], 'hilbert',[0:50])
close all

N = 100000; % time series length (number of time samples)

r = corrcoeff;
[s1, s2] = corr2_signal(r, N);
[s1_f, s2_f] = filter_signal(s1, s2, filter);
[s1_n, s2_n] = add_noise(s1_f, s2_f, noise, N);
fprintf('SNR: S1: %.02f, S2: %.02f\n', snr(s1_n, s1_f), snr(s2_n, s2_f))

corr([s1_f;s2_f]')

if strcmp(phase_method,'both')
    plvs_hilbert = hilbert_plv(s1_n, s2_n, freq);
    plvs_morlet = morlet_plv(s1_n, s2_n, freq);
    PLV = [cell2mat(plvs_hilbert); cell2mat(plvs_morlet)];
elseif strcmp(phase_method,'morlet')
    plvs_morlet = morlet_plv(s1_n, s2_n, freq);
    PLV = cell2mat(plvs_morlet);
else
    plvs_hilbert = hilbert_plv(s1_n, s2_n, freq);
    PLV = cell2mat(plvs_hilbert);
end

% PLOTTING

if do_plot
    plot_raw_signals(s1, s2, r);
    plot_filtered_signal(s1_f, s2_f)
    plot_with_noise(s1_n, s2_n)
end

if strcmp(phase_method,'both') && do_plot
    plot_plv_dual(plvs_hilbert, plvs_morlet, freq)
elseif strcmp(phase_method,'morlet') && do_plot
    plot_plv(plvs_morlet, freq, 'morlet')
elseif strcmp(phase_method,'hilbert') && do_plot
    plot_plv(plvs_hilbert, freq, 'hilbert')
end


    function [s1, s2] = corr2_signal(r, N)
        
        m1=0;m2=0; % means
        sd1=1; sd2=1; % standard deviations
        u=randn(1,N); % Gaussian time series, mean=0, sd=1;
        v=randn(1,N); % 2nd Gaussian time series (independent of u)
        
        s1=sd1*u+m1; % random time series with mean=m1, sd=s1;
        s2=sd2*(r*u+sqrt(1-r^2)*v)+m2; % random time series correlated with x (if |r|>0)
    end


    function [s1_n, s2_n] = add_noise(s1, s2, noise, N)
        
        noise1 = noise(1);
        noise2 = noise(2);
        
        N1 = (noise1 .* rand(1,N));
        N2 = (noise2 * rand(1,N));
        
        N1 = N1 - mean(N1);
        N2 = N2 - mean(N2);
        
        s1_n = s1 + N1;
        s2_n = s2 + N2;
    end

    function [s1_f, s2_f] = filter_signal(s1, s2, filter)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% FILTER CUTOFF FREQS %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lowFreq = 20;
        hiFreq = 25;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fs = 100;
        order = 5;
        [b,a] = butter(order, [lowFreq hiFreq]/(fs/2), 'bandpass');
        
        if filter(1)==0
            s1_f = s1;
        else
            s1_f = filtfilt(b, a, s1);
        end
        
        if filter(2)==0
            s2_f = s2;
        else
            s2_f = filtfilt(b, a, s2);
        end
    end

    function plvs = hilbert_plv(s1, s2, freq)
        plvs = hilbert_compute_plv(s1, s2, freq);
    end

    function plv = hilbert_compute_plv(s1, s2, freq)
        
        plv = cell(0);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% HILBERT WINDOW WIDTH %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        width = 8;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        half_width = width/2;
        
        for j = freq
            
            [b,a] = get_filter(j-half_width, j+half_width);
            s1_filt = filtfilt(b, a, s1);
            s2_filt = filtfilt(b, a, s2);
            complex1=(hilbert(s1_filt));
            complex2=(hilbert(s2_filt));
            
            p1 = conj( complex1 ./ sqrt(sum(abs(complex1).^2,1)) );
            
            p2 = complex2 ./ sqrt(sum(abs(complex2).^2,1));
            
            diff_vectors            = p1 .* p2;
            int_PLV                 = abs(mean(diff_vectors,2));
            plv{j}                  = int_PLV;
        end
    end

    function [b,a] = get_filter(lowFreq, hiFreq)
%         lowFreq = 20;
%         hiFreq = 25;
        fs = 100;
        order = 1;
        
        if hiFreq >= fs/2
            hiFreq = (fs/2) - 0.00001
        end
        
        if lowFreq <= 0
            lowFreq = 0.001
        end
        
        [b,a] = butter(order, [lowFreq hiFreq]/(fs/2), 'bandpass');
    end

    function plvs = morlet_plv(s1, s2, freq)
        plvs = plvm.PLVM.resting_plv_pairwise_stacked([s1; s2], 100, freq);
    end


    function plot_raw_signals(s1, s2, r)
        figure(1);
        subplot(2,2,1)
        plot(s1)
        grid on
        title('Signal 1 in Time Domain')
        subplot(2,2,2)
        %  pwelch(s1,2,1,'twosided')
        pwelch(s1, [], [], [], 100)
        subplot(2,2,3)
        plot(s2(1,:))
        grid on
        title(['Signal 2 for Corr coeff ' num2str(4)])
        subplot(2,2,4)
        pwelch(s2(1,:), [], [], [], 100)
        
        figure(2)
        plot(s2(1,:))
        hold on
        plot(s1,'r')
        grid on
        title('Time doamin signals')
        legend('signal2', 'signal1')
    end

    function plot_filtered_signal(s1_f, s2_f)
        figure(3)
        subplot(2,2,1)
        plot(s1_f), grid on
        title('BandPass Filtered Signal1 in Time Domain ')
        subplot(2,2,2)
        pwelch(s1_f, [], [], [], 100)
        subplot(2,2,3)
        plot(s2_f), grid on
        title('BandPass Filtered Signal2 in Time Domain ')
        subplot(2,2,4)
        pwelch(s2_f, [], [], [], 100)
    end

    function plot_with_noise(s1_n, s2_n)
        figure(4)
        subplot(2,2,1)
        plot(s1_n),grid on
        title('Addition of Noise in signal1 ')
        subplot(2,2,2)
        pwelch(s1_n, [], [], [], 100)
        
        subplot(2,2,3)
        plot(s2_n),grid on
        title('Addition of Noise in signal2 ')
        subplot(2,2,4)
        pwelch(s2_n, [], [], [], 100)
        
    end


    function plot_plv(plvs, freq, phase_type)
        
        x = cell2mat(plvs);
        figure(6)
        plot(freq, x)
        title('PLV')
        xlabel('frequency')
        ylabel('PLV Value')
        legend(phase_type)
        
    end

    function plot_plv_dual(plvs_hilbert, plvs_morlet, freq)
        
        xhilbert = cell2mat(plvs_hilbert);
        xmorlet = cell2mat(plvs_morlet);
        figure(6)
        plot(freq, xhilbert)
        hold on;
        plot(freq, xmorlet)
        title('PLV')
        xlabel('frequency')
        ylabel('PLV Value')
        legend('Hilbert', 'Morlet')
        
    end

end