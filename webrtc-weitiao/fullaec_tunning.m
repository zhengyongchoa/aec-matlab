% Partitioned block frequency domain adaptive filtering NLMS and
% standard time-domain sample-based NLMS
%near is micphone captured signal
% fid=fopen('near.pcm', 'rb'); % Load far end
% ssin=fread(fid,inf,'float32');
% fclose(fid);
% %far is speaker played music
% fid=fopen('far.pcm', 'rb'); % Load fnear end
% rrin=fread(fid,inf,'float32');
% fclose(fid);

[rrin,fs] = audioread('ref16k.wav');
[ssin,fs] = audioread('mic16k.wav');


ssin = ssin.*32768;
rrin = rrin.*32768;

rand('state',13);
fs=16000;
mult=fs/8000;
if fs == 8000
cohRange = 2:3;
elseif fs==16000
cohRange = 2;
end

% Flags
NLPon=1; % NLP on
CNon=0; % Comfort noise on
PLTon=0; % Plotting on

M = 16; % Number of partitions
N = 64; % Partition length
L = M*N; % Filter length
if fs == 8000
    mufb = 0.6;
else
    mufb = 0.8;
end
VADtd=48;
alp = 0.15; % Power estimation factor 
alc = 0.1; % Coherence estimation factor
beta = 0.9; % Plotting factor
%% Changed a little %%
step = 0.1875;%0.1875; % Downward step size
%%
if fs == 8000
    threshold=2e-6; % DTrob threshold
else
    %threshold=0.7e-6;
    threshold=1.5e-6; 
end

if fs == 8000
    echoBandRange = ceil(300*2/fs*N):floor(1800*2/fs*N);
else
    echoBandRange = ceil(60*2/fs*N):floor(1500*2/fs*N);
end
suppState = 1;
transCtr = 0;

Nt=1;
vt=1;

ramp = 1.0003; % Upward ramp
rampd = 0.999; % Downward ramp
cvt = 20; % Subband VAD threshold;
nnthres = 20; % Noise threshold

shh=logspace(-1.3,-2.2,N+1)';
sh=[shh;flipud(shh(2:end-1))]; % Suppression profile

len=length(ssin);
w=zeros(L,1); % Sample-based TD(time domain) NLMS
WFb=zeros(N+1,M); % Block-based FD(frequency domain) NLMS
WFbOld=zeros(N+1,M); % Block-based FD NLMS
YFb=zeros(N+1,M);
erfb=zeros(len,1);
erfb3=zeros(len,1);

ercn=zeros(len,1);
zm=zeros(N,1);
XFm=zeros(N+1,M);
YFm=zeros(N+1,M);
pn0=10*ones(N+1,1);
pn=zeros(N+1,1);
NN=len;
Nb=floor(NN/N)-M;
erifb=zeros(Nb+1,1)+0.1;
erifb3=zeros(Nb+1,1)+0.1;
ericn=zeros(Nb+1,1)+0.1;
dri=zeros(Nb+1,1)+0.1;
start=1;
xo=zeros(N,1);
do=xo;
eo=xo;

echoBands=zeros(Nb+1,1);
cohxdAvg=zeros(Nb+1,1);
cohxdSlow=zeros(Nb+1,N+1);
cohedSlow=zeros(Nb+1,N+1);
%overdriveM=zeros(Nb+1,N+1);
cohxdFastAvg=zeros(Nb+1,1);
cohxdAvgBad=zeros(Nb+1,1);
cohedAvg=zeros(Nb+1,1);
cohedFastAvg=zeros(Nb+1,1);
hnledAvg=zeros(Nb+1,1);
hnlxdAvg=zeros(Nb+1,1);
ovrdV=zeros(Nb+1,1);
dIdxV=zeros(Nb+1,1);
SLxV=zeros(Nb+1,1);
hnlSortQV=zeros(Nb+1,1);
hnlPrefAvgV=zeros(Nb+1,1);
mutInfAvg=zeros(Nb+1,1);
%overdrive=zeros(Nb+1,1);
hnled = zeros(N+1, 1);
weight=zeros(N+1,1);
hnlMax = zeros(N+1, 1);
hnl = zeros(N+1, 1);
overdrive = ones(1, N+1);
xfwm=zeros(N+1,M);
dfm=zeros(N+1,M);
WFbD=ones(N+1,1);

fbSupp = 0;
hnlLocalMin = 1;
cohxdLocalMin = 1;
hnlLocalMinV=zeros(Nb+1,1);
cohxdLocalMinV=zeros(Nb+1,1);
hnlMinV=zeros(Nb+1,1);
dkEnV=zeros(Nb+1,1);
ekEnV=zeros(Nb+1,1);
ovrd = 2;
ovrdPos = floor((N+1)/4);
ovrdSm = 2;
hnlMin = 1;
minCtr = 0;
SeMin = 0;
SdMin = 0;
SeLocalAvg = 0;
SeMinSm = 0;
divergeFact = 1;
dIdx = 1;
hnlMinCtr = 0;
hnlNewMin = 0;
divergeState = 0;

Sy=ones(N+1,1);
Sym=1e7*ones(N+1,1);

wins=[0;sqrt(hanning(2*N-1))];
ubufn=zeros(2*N,1);
ebuf=zeros(2*N,1);
ebuf2=zeros(2*N,1);
ebuf4=zeros(2*N,1);
mbuf=zeros(2*N,1);

cohedFast = zeros(N+1,1);
cohxdFast = zeros(N+1,1);
cohxd = zeros(N+1,1);
Se = zeros(N+1,1);
Sd = zeros(N+1,1);
Sx = zeros(N+1,1);
SxBad = zeros(N+1,1);
Sed = zeros(N+1,1);
Sxd = zeros(N+1,1);
SxdBad = zeros(N+1,1);
hnledp=[];

cohxdMax = 0;

for kk=1:Nb
    pos = N * (kk-1) + start;
    
    % FD block method
    % ---------------------- Organize data
    
    %far is speaker played music
    xk = rrin(pos:pos+N-1);
    %near is micphone captured signal
    dk = ssin(pos:pos+N-1);
    
    %----------------------- far end signal process
    xx = [xo;xk];
    xo = xk;
    tmp = fft(xx);
    XX = tmp(1:N+1);

    dd = [do;dk]; % Overlap
    do = dk;
    tmp = fft(dd); % Frequency domain
    DD = tmp(1:N+1);
    
    % ------------------------far end Power estimation
    pn0 = (1 - alp) * pn0 + alp * real(XX.* conj(XX));
    pn = pn0;
%   pn = (1 - alp) * pn + alp * M * pn0;
    
    % ---------------------- Filtering
    XFm(:,1) = XX;
    for mm=0:(M-1)
        m=mm+1;
        YFb(:,m) = XFm(:,m) .* WFb(:,m);
    end
    yfk = sum(YFb,2);
    tmp = [yfk ; flipud(conj(yfk(2:N)))];
    ykt = real(ifft(tmp));
    ykfb = ykt(end-N+1:end);
    
    % ---------------------- Error estimation
    ekfb = dk - ykfb;
    %if sum(abs(ekfb)) < sum(abs(dk))
        %ekfb = dk - ykfb;
    % erfb(pos:pos+N-1) = ekfb;
    %else
        %ekfb = dk;
    % erfb(pos:pos+N-1) = dk;
    %end
%(kk-1)*(N*2)+1
    erfb(pos:pos+N-1) = ekfb;
    tmp = fft([zm;ekfb]); % FD version for cancelling part (overlap-save)
    Ek = tmp(1:N+1);

    % ------------------------ Adaptation
    %Ek2 = Ek ./(M*pn + 0.001); % Normalized error
    Ek2 = Ek ./(pn + 0.001); % Normalized error
    
    absEf = max(abs(Ek2), threshold);
    absEf = ones(N+1,1)*threshold./absEf;
    Ek2 = Ek2.*absEf;

    mEk = mufb.*Ek2;
    PP = conj(XFm).*(ones(M,1) * mEk')';
    tmp = [PP ; flipud(conj(PP(2:N,:)))];
    IFPP = real(ifft(tmp));
    PH = IFPP(1:N,:);
    tmp = fft([PH;zeros(N,M)]);
    FPH = tmp(1:N+1,:);
    WFb = WFb + FPH;

    %if mod(kk, 10*mult) == 0
    WFbEn = sum(real(WFb.*conj(WFb)));
    %WFbEn = sum(abs(WFb));
    [tmp, dIdx] = max(WFbEn);

    WFbD = sum(abs(WFb(:, dIdx)),2);
    %WFbD = WFbD / (mean(WFbD) + 1e-10);
    WFbD = min(max(WFbD, 0.5), 4);
    %     end
    dIdxV(kk) = dIdx;
    
    % NLP
    if (NLPon)  
        ee = [eo;ekfb];
        eo = ekfb;
        window = wins;
        if fs == 8000
            gamma = 0.9;
        else
        gamma = 0.93;
        end

        tmp = fft(xx.*window);
        xf = tmp(1:N+1);
        tmp = fft(dd.*window);
        df = tmp(1:N+1);
        tmp = fft(ee.*window);
        ef = tmp(1:N+1);

        xfwm(:,1) = xf;
        xf = xfwm(:,dIdx);
        %fprintf(1,'%d: %f\n', kk, xf(4));
        dfm(:,1) = df;
        
        SxOld = Sx;

        Se = gamma*Se + (1-gamma)*real(ef.*conj(ef));
        Sd = gamma*Sd + (1-gamma)*real(df.*conj(df));
        Sx = gamma*Sx + (1 - gamma)*real(xf.*conj(xf));



        % coherence
        Sxd = gamma*Sxd + (1 - gamma)*xf.*conj(df);
        Sed = gamma*Sed + (1-gamma)*ef.*conj(df);

        cohed = real(Sed.*conj(Sed))./(Se.*Sd + 1e-10);
        cohedAvg(kk) = mean(cohed(echoBandRange));

        cohxd = real(Sxd.*conj(Sxd))./(Sx.*Sd + 1e-10);
        freqSm = 0.6;
        cohxd(2:end) = filter(freqSm, [1 -(1-freqSm)], cohxd(2:end));
        cohxd(end:2) = filter(freqSm, [1 -(1-freqSm)], cohxd(end:2));
        cohxdAvg(kk) = mean(cohxd(echoBandRange));

        hnled = min(1 - cohxd, cohed);


        if kk > 1
            cohxdSlow(kk,:) = 0.99*cohxdSlow(kk-1,:) + 0.01*cohxd';
            cohedSlow(kk,:) = 0.99*cohedSlow(kk-1,:) + 0.01*(1-cohed)';
        end


        if 0
            hnlMax = hnlMax*0.9999;
            hnlMax = max(hnlMax, hnled);
            avgHnl = mean(hnlMax(echoBandRange));
            if avgHnl > 0.3
                overdrive = max(log(avgHnl)/log(0.99), 1);
            end
            weight(4:end) = max(hnlMax) - hnlMax(4:end);
        end
        
        

        cohedMean = mean(cohed(echoBandRange));
        [hnlSort, 	hnlSortIdx] = sort(1-cohxd(echoBandRange));
        [xSort, xSortIdx] = sort(Sx);
        hnlSortQ = mean(1 - cohxd(echoBandRange));

        [hnlSort2, hnlSortIdx2] = sort(hnled(echoBandRange));
        hnlQuant = 0.75;
        hnlQuantLow = 0.5;
        qIdx = floor(hnlQuant*length(hnlSort2));
        qIdxLow = floor(hnlQuantLow*length(hnlSort2));
        hnlPrefAvg = hnlSort2(qIdx);
        hnlPrefAvgLow = hnlSort2(qIdxLow);
   

        if cohedMean > 0.98 && hnlSortQ > 0.9
            suppState = 0;
        elseif cohedMean < 0.95 || hnlSortQ < 0.8
            suppState = 1;
        end

        if hnlSortQ < cohxdLocalMin && hnlSortQ < 0.75
            cohxdLocalMin = hnlSortQ;
        end

        if cohxdLocalMin == 1
            ovrd = 3;
            hnled = 1-cohxd;
            hnlPrefAvg = hnlSortQ;
            hnlPrefAvgLow = hnlSortQ;
        end

        if suppState == 0
            hnled = cohed;
            hnlPrefAvg = cohedMean;
            hnlPrefAvgLow = cohedMean;
        end
            
        if hnlPrefAvgLow < hnlLocalMin && hnlPrefAvgLow < 0.6
            hnlLocalMin = hnlPrefAvgLow;
            hnlMin = hnlPrefAvgLow;
            hnlNewMin = 1;
            hnlMinCtr = 0;
            if hnlMinCtr == 0
                hnlMinCtr = hnlMinCtr + 1;
            else
                hnlMinCtr = 0;
                hnlMin = hnlLocalMin;
                SeLocalMin = SeQ;
                SdLocalMin = SdQ;
                SeLocalAvg = 0;
                minCtr = 0;
                ovrd = max(log(0.0001)/log(hnlMin), 2);
                divergeFact = hnlLocalMin;
            end
        end

        if hnlNewMin == 1
            hnlMinCtr = hnlMinCtr + 1;
        end
        if hnlMinCtr == 2
            hnlNewMin = 0;
            hnlMinCtr = 0;
            ovrd = max(log(0.00000001)/(log(hnlMin + 1e-10) + 1e-10), 5);
        end
        hnlLocalMin = min(hnlLocalMin + 0.0008/mult, 1);
        cohxdLocalMin = min(cohxdLocalMin + 0.0004/mult, 1);

        if ovrd < ovrdSm
            ovrdSm = 0.99*ovrdSm + 0.01*ovrd;
        else
            ovrdSm = 0.9*ovrdSm + 0.1*ovrd;
        end

        ekEn = sum(Se);
        dkEn = sum(Sd);

        if divergeState == 0
            if ekEn > dkEn
                ef = df;
                divergeState = 1;
            end
        else
            if ekEn*1.05 < dkEn
                divergeState = 0;
            else
                ef = df;
            end
        end

        if ekEn > dkEn*19.95
            WFb=zeros(N+1,M); % Block-based FD NLMS
        end

        ekEnV(kk) = ekEn;
        dkEnV(kk) = dkEn;

        hnlLocalMinV(kk) = hnlLocalMin;
        cohxdLocalMinV(kk) = cohxdLocalMin;
        hnlMinV(kk) = hnlMin;
    
        aggrFact = 0.3;
  
        wCurve = [0; aggrFact*sqrt(linspace(0,1,N))' + 0.1];

        weight = wCurve;
  

        hnled = weight.*min(hnlPrefAvg, hnled) + (1 - weight).*hnled;

    

        od = ovrdSm*(sqrt(linspace(0,1,N+1))' + 1);
   
        sshift = ones(N+1,1);

        hnled = hnled.^(od.*sshift);

    
        hnl = hnled;

        ef = ef.*(hnl);


        
        ovrdV(kk) = ovrdSm;

        hnledAvg(kk) = 1-mean(1-cohed(echoBandRange));
        hnlxdAvg(kk) = 1-mean(cohxd(echoBandRange));

        hnlSortQV(kk) = hnlPrefAvgLow;
        hnlPrefAvgV(kk) = hnlPrefAvg;


        % Comfort noise
        if (CNon)
            snn=sqrt(Sym);
            snn(1)=0; % Reject LF noise
            Un=snn.*exp(j*2*pi.*[0;rand(N-1,1);0]);

            % Weight comfort noise by suppression
            Un = sqrt(1-hnled.^2).*Un;
            Fmix = ef + Un;
        else
            Fmix = ef;
        end

        % Overlap and add in time domain for smoothness
        tmp = [Fmix ; flipud(conj(Fmix(2:N)))];
        mixw = wins.*real(ifft(tmp));
        mola = mbuf(end-N+1:end) + mixw(1:N);
        mbuf = mixw;
        ercn(pos:pos+N-1) = mola;%%%%%-------------you can hear the effect by sound(10*ercn,16000),add by Shichaog
    
    end % NLPon

    % Filter update
    Ek2 = Ek ./(pn + 0.001); % Normalized error
    %Ek2 = Ek ./(100*pn + 0.001); % Normalized error




% Shift old FFTs
    XFm(:,2:end) = XFm(:,1:end-1);
    YFm(:,2:end) = YFm(:,1:end-1);
    xfwm(:,2:end) = xfwm(:,1:end-1);
    dfm(:,2:end) = dfm(:,1:end-1);

end
audiowrite('full_out_linear.wav',erfb./32768,16000);
audiowrite('full_out.wav',ercn./32768,16000);
%figure(2);
%plot([feat(:,1) feat(:,2)+1 feat(:,3)+2 mfeat+3]);
%plot([feat(:,1) mfeat+1]);

%figure(3);
%plot(10*log10([dri erifb erifb3 ericn]));
%legend('Near-end','Error','Post NLP','Final',4);
% Compensate for delay
%ercn=[ercn(N+1:end);zeros(N,1)];
%ercn_=[ercn_(N+1:end);zeros(N,1)];

%figure(11);
%plot(cohxdSlow);

%figure(12);
%surf(cohxdSlow);
%shading interp;

%figure(13);
%plot(overdriveM);

%figure(14);
%surf(overdriveM);
%shading interp;

figure(10);
t = (0:Nb)*N/fs;
rrinSubSamp = rrin(N*(1:(Nb+1)));
plot(t, rrinSubSamp/max(abs(rrinSubSamp)),'b');
hold on
plot(t, hnledAvg, 'r');
plot(t, hnlxdAvg, 'g');
plot(t, hnlSortQV, 'y');
plot(t, hnlLocalMinV, 'k');
plot(t, cohxdLocalMinV, 'c');
plot(t, hnlPrefAvgV, 'm');
%plot(t, cohxdAvg, 'r');
%plot(cohxdFastAvg, 'r');
%plot(cohxdAvgBad, 'k');
%plot(t, cohedAvg, 'k');
%plot(t, 1-cohedFastAvg, 'k');
%plot(ssin(N*(1:floor(length(ssin)/N)))/max(abs(ssin)));
%plot(echoBands,'r');
%plot(overdrive, 'g');
%plot(erfb(N*(1:floor(length(erfb)/N)))/max(abs(erfb)));
hold off
%tight x;

% figure(11)
% plot(t, ovrdV);
%tightx;
%plot(mfeat,'r');
%plot(1-cohxyp_,'r');
%plot(Hnlxydp,'y');
%plot(hnledp,'k');
%plot(Hnlxydp, 'c');
%plot(ccohpd_,'k');
%plot(supplot_, 'g');
%plot(ones(length(mfeat),1)*rr1_, 'k');
%plot(ones(length(mfeat),1)*rr2_, 'k');
%plot(N*(1:length(feat)), feat);
%plot(Sep_,'r');
%axis([1 floor(length(erfb)/N) -1 1])
%hold off
%plot(10*log10([Se_, Sx_, Seu_, real(sf_.*conj(sf_))]));
%legend('Se','Sx','Seu','S');
%figure(5)
%plot([ercn ercn_]);

% figure(12)
% plot(t, dIdxV);
%plot(t, SLxV);
%tightx;

%figure(13)
%plot(t, [ekEnV dkEnV]);
%plot(t, dkEnV./(ekEnV+1e-10));
%tightx;

%close(hh);
%spclab(fs,ssin,erfb,ercn,'outxd.pcm');
%spclab(fs,rrin,ssin,erfb,1.78*ercn,'vqeOut-1.pcm');
%spclab(fs,erfb,'aecOutLp.pcm');
%spclab(fs,rrin,ssin,erfb,1.78*ercn,'aecOut25.pcm','vqeOut-1.pcm');
%spclab(fs,rrin,ssin,erfb,ercn,'aecOut-mba.pcm');
%spclab(fs,rrin,ssin,erfb,ercn,'aecOut.pcm');
%spclab(fs, ssin, erfb, ercn, 'out0.pcm');