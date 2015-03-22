function PLV = run_plv( corrcoeff, Noiseampli,noiseadd, filtersignl,filter,freq)

%%%% in this variable corrcoeff or Noiseampli can't be provided lower and
%%%% upper range at a time.
%%% provision two add noise and filter the desired signal has been kept in
%%% this script. for this 'noiseadd' is selected to add noise in the
%%% desired signal. for example noiseadd=[0,1] will add noise to 2nd signal
%%% only and similarly [0,1]= for 2nd signal noise, [1,0]= for 1st signal
%%% noise, [1,1]= noise in both signals and [0,0] = noise in none
%%% similar pattern for sleection of signal for filteration is opted in 'filtersignal'

%%%example    run_plv([0.1 0], [1 3],[1,0],[1,1], 'hilbert',[0:50])
close all

if corrcoeff(2)~=0 && Noiseampli(2)~=0
    msgbox('Enter the range of either corr coeff or noise amplitude','Error','error')
    return
end

if corrcoeff(2)~=0 && corrcoeff(2)< corrcoeff(1) || Noiseampli(2)~=0 && Noiseampli(2)< Noiseampli(1)
    msgbox('Enter the range in correct order','Error','error')
    return
end

%% Signal generation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r1=corrcoeff(1); % lower range value
r2=corrcoeff(2); % upper range value

N=10000; % time series length (number of time samples)
m1=0;m2=0; % means
sd1=1; sd2=2; % standard deviations
u=randn(1,N); % Gaussian time series, mean=0, sd=1;
v=randn(1,N); % 2nd Gaussian time series (independent of u)
temp=1;
if r1<r2
    for ii=r1:0.1:r2
        s1=sd1*u+m1; % random time series with mean=m1, sd=s1;
        s2(temp,:)=sd2*(ii*u+sqrt(1-ii^2)*v)+m2; % random time series correlated with x (if |r|>0)
        temp=temp+1;
    end
else
    s1=sd1*u+m1; % random time series with mean=m1, sd=s1;
    s2=sd2*(r1*u+sqrt(1-r1^2)*v)+m2; % random time series correlated with x (if |r|>0)
end
figure(1);subplot(2,2,1)
plot(s1),grid on
title('Signal 1 in Time Domain')
subplot(2,2,2)
%  pwelch(s1,2,1,'twosided')
pwelch(s1)
subplot(2,2,3)
plot(s2(1,:)),grid on
  title(['Signal 2 for Corr coeff ' num2str(r1)])
  subplot(2,2,4)
  pwelch(s2(1,:))
  
figure(2)
plot(s2(1,:))
hold on
plot(s1,'r')
grid on
title('Time doamin signals')
legend('signal2', 'signal1')
  %% Noise addition
temp=1;
Lmag=Noiseampli(1);
Hmag=Noiseampli(2);
if noiseadd(1)==0 && noiseadd(2)==1
    if Lmag<Hmag
        for ii=Lmag:1:Hmag
            s11(temp,:)=s1; %%% signal1
            s22(temp,:)=s2+ii*rand(1,N); %%% signal2
            temp=temp+1;
        end
    else
        s11=s1; %%% signal1
        for jj=1:size(s2,1)
            s22(temp,:)=s2(temp,:)+Lmag*rand(1,N); %%% signal2
            temp=temp+1;
        end
    end
end
if noiseadd(1)==1 && noiseadd(2)==1
    if Lmag<Hmag
        for ii=Lmag:1:Hmag
            s11(temp,:)=s1+ii*rand(1,N); %%% signal1
            s22(temp,:)=s2+ii*rand(1,N); %%% signal2
            temp=temp+1;
        end
    else
        s11=s1+Lmag*rand(1,N); %%% signal1
        for jj=1:size(s2,1)
            s22(temp,:)=s2(temp,:)+Lmag*rand(1,N); %%% signal2
            temp=temp+1;
        end
    end
end
if noiseadd(1)==1 && noiseadd(2)==0
    if Lmag<Hmag
        for ii=Lmag:1:Hmag
            s11(temp,:)=s1+ii*rand(1,N); %%% signal1
            s22(temp,:)=s2; %%% signal2
            temp=temp+1;
        end
    else
        s11=s1+Lmag*rand(1,N); %%% signal1
        s22=s2;
    end
end
if noiseadd(1)==0 && noiseadd(2)==0
    s11=s1;s22=s2;
end

figure(3)
subplot(2,2,1)
plot(s11(1,:)),grid on
title('Addition of Noise in signal1 ')
subplot(2,2,2)
pwelch(s11(1,:),2,1,'twosided')

subplot(2,2,3)
plot(s22(1,:)),grid on
title('Addition of Noise in signal2 ')
subplot(2,2,4)
pwelch(s22(1,:),2,1,'twosided')

s1=s11;s2=s22;
  %% Filter
%%% for signal 1  
lowFreq = 20;
hiFreq = 25;
fs = 100;
order = 5;
[b,a] = butter(order, [lowFreq hiFreq]/(fs/2), 'bandpass');

if filtersignl(1)==0
    out_signal=s1;
else
    for ii=1:size(s1,1)
        out_signal(ii,:) = filtfilt(b,a,s1(ii,:));
    end
end
figure(4)
subplot(2,2,1)
plot(out_signal(ii,:)), grid on
title('BandPass Filtered Signal1 in Time Domain ')
subplot(2,2,2)  
pwelch(out_signal(ii,:))
%%% for signal 2  
if filtersignl(2)==0 
    out_signal=s2;
else
    for ii=1: size(s2,1)
        out_signal(ii,:) = filtfilt(b,a,s2(ii,:));
    end
end
subplot(2,2,3)
plot(out_signal(1,:)), grid on
title('BandPass Filtered Signal2 in Time Domain ')
subplot(2,2,4)  
pwelch(out_signal(1,:))

%% phase angle

if strcmp(filter,'hilbert')
    phase1=(angle(hilbert(s1)));
    phase2=(angle(hilbert(s2)));
    figure(5)
    plot(phase1(1,:),'LineWidth',2)
    hold on
    plot(phase2(1,:),'r','LineWidth',2)
    grid on
    title('Phase angle of signal in radians')
    legend('Signal1','Signal2')
end
if strcmp(filter,'morlet')
    h = waitbar(0,'Please wait...');
        steps=size(s2,1);
        for ii=1:size(s2,1)    
            for fb=freq(1):freq(2)
                [PSI,X]=cmorwavf(-8,8,N,fb,1);
                phase2(fb,:)=unwrap(angle(conv(PSI,s2(ii,:))));                
            end
            temp.phase2(ii)={phase2};
            waitbar(ii/(steps))
        end
        steps=size(s1,1);
        for ii=1:size(s1,1)    
            for fb=freq(1):freq(2)
                [PSI,X]=cmorwavf(-8,8,N,fb,1);
                phase1(fb,:)=unwrap(angle(conv(PSI,s1(ii,:))));
            end
            temp.phase1(ii)={phase1};
            waitbar(ii/steps)
        end
        close(h)
        figure(5)
        plot(X,PSI)
        grid on
        title('morlet transform at frequency 50 Hz')
        figure(6)
        subplot(2,1,1)
        plot(phase1(fb,:))
        hold on
        plot(phase2(fb,:),'r')
        title('Phase angle of signal in radians')
        legend('Signal1','Signal2')
        subplot(2,1,2)
        plot(real(conv(PSI,s1(1,:))))
        hold on
        plot(real(conv(PSI,s2(1,:))),'r')
        grid on
        title('Envelope of Signal')
        legend('Signal1','Signal2')
        phase2=temp.phase2;
        phase1=temp.phase1;
end

%% PLV calculation
color=['r','g','b','c','m','y','k'];
if strcmp(filter,'hilbert')
    for ii=1:size(phase2,1)
        if size(phase2,1)==size(phase1,1)
            diff=phase1(ii,:)-phase2(ii,:);
        else
            diff=phase1(size(phase1,1),:)-phase2(ii,:);
        end
        plv(ii,:)=abs(1i*exp(diff))/(length(diff));
        figure(7)
        plot(plv(ii,1:100),color(ii),'linewidth',2)
        hold on
        grid on
    end
    title('Phase Locking Value for Hilbert Transform')
    xlabel('no of samples')
    ylabel('PLV')
end
if strcmp(filter,'morlet')
    for ii=1:length(phase2)
        if length(phase1)==length(phase2)
            diff=phase1{ii}-phase2{ii};
        else
            diff=phase1{length(phase1)}-phase2{ii};
        end
        plv(ii,:)=abs(sum(1i*exp(diff)))/(50*length(diff));
        figure(7)
        plot(plv(ii,1:100),color(ii),'linewidth',2)
        hold on
        grid on
    end
title('Phase Locking Value for Morlet Wavelet for 1-50 Hz frequency')
xlabel('no of samples')
    ylabel('PLV')
end
legend('-dynamiclegend')
end