% In this demo, we assume that sdrrx and sdrtx are open, and are fully
% calibrated. Look at the calibration demo to make sure this is done. In
% the minimum, the timing and phase offsets need to be calibrated.

sdrtx = sdr0;
sdrrx = sdr1;

sdrtx.nch = 8;
sdrrx.nch = 8;
nFFT = 1024;
nread = nFFT;
nskip = nFFT*5;
ntimes = 50;
txfd = zeros(nFFT, 1);
constellation = [1+1j 1-1j -1+1j -1-1j];

sdrrx.set_switches('off');

scMin = 100; scMax = 100;
for scIndex = scMin:scMax
    txfd(nFFT/2 + 1 + scIndex) = constellation(randi(4));
end
txfd = fftshift(txfd);
txtd = ifft(txfd);
m = max(abs(txtd));
txtd = txtd / m * 50;

naod = 41;
aods = linspace(-1, 1, naod);

pArray = zeros(1, naod);

for iaod = 1:naod
    p = 0;
    fprintf('.');
    txtdMod = zeros(nFFT, sdrrx.nch);
    aod = aods(iaod);
    rad2deg(aod)
    for txIndex=1:sdrrx.nch
        txtdMod(:, txIndex) = txtd * exp(1j*txIndex*pi*sin(aod)); % Apply BF
    end
    txtdMod = sdrrx.applyCalTxArray(txtdMod);
    sdrrx.send(txtdMod);

    rxtd = sdrtx.recv(nread, nskip, ntimes);
    rxtd = sdrtx.applyCalRxArray(rxtd);
        
    for itimes = 1:ntimes
        refRxIndex = 1;
        fd = fftshift(fft(rxtd(:, itimes, refRxIndex)));
        p = p + sum(abs(fd( nFFT/2 + 1 + scMin : nFFT/2 + 1 + scMax)));
    end %itimes
    pArray(iaod) = p;
end % iaoa

% Plot
pArray = pArray / max(pArray);
figure(3); clf;
plot(rad2deg(aods), mag2db(pArray));
xlabel('Angle of Departure (Deg)');
ylabel('Power (dB)');
grid on; grid minor;
ylim([-20 0])

% Stop transmitting and do a dummy read
txtd = zeros(nFFT, sdrrx.nch);
sdrrx.send(txtd);
sdrrx.recv(nread, nskip, ntimes);
sdrrx.recv(nread, nskip, ntimes);

% Clear workspace variables
clear aoa aoas fd iaoa naoa p pArray refTxIndex td tdbf txtdMod;
clear ans itimes m nFFT nread nskip ntimes rxIndex rxtd;
clear scIndex txfd txtd constellation scMax scMin;
