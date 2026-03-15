clear; clc; close all;

load("raytraceSARdata.mat")
disp("Example rayData entry:")
disp(rayData(1,1))

%constants
c = physconst('LightSpeed');
lambda = c/fc;

%radar parameters
rangeResolution = 3;
crossRangeResolution = 3;

bw = c/(2*rangeResolution);

prf = 1000;
tpd = 3e-6;
fs = 120e6;

maxRange = 200;

Nr = ceil((2*maxRange/c)*fs);

fastTime = (0:Nr-1)/fs;

%waveform
waveform = phased.LinearFMWaveform( ...
    'SampleRate',fs,...
    'PulseWidth',tpd,...
    'PRF',prf,...
    'SweepBandwidth',bw);

%generating raw SAR data
rxsig = zeros(Nr,Np);

numTargets = numel(latTargets);

for n = 1:Np

    txpulse = waveform();
    txpulse = txpulse(1:Nr);

    for t = 1:numTargets

        if ~isfield(rayData(n,t),'R')
            continue
        end

        for k = 1:length(rayData(n,t).R)

            R = rayData(n,t).R(k);

            tau = 2*R/c;

            delaySamples = round(tau*fs);
            if n==1 && t==1
                disp("delaysamples:")
                disp(delaySamples)
            end

            if delaySamples > Nr
                continue
            end

            A = rayData(n,t).amplitude(k);
            ph = rayData(n,t).phase(k);

            echo = A * ph * txpulse;

            idx = delaySamples + (1:length(echo));
            idx(idx>Nr) = [];

            rxsig(idx,n) = rxsig(idx,n) + echo(1:length(idx));
            if n==1 && t==1
                disp("max echo")
                disp(max(abs(echo)))
            end

        end

    end

end
disp("max raydata")
disp(max([rayData(1,:).R]))

% disp("Maximum raw signal amplitude:")
% disp(max(abs(rxsig(:))))
% figure
% imagesc(abs(rxsig))
% title("Raw SAR Data")
% xlabel("Pulse")
% ylabel("Range bin")
% colorbar
figure
plot(abs(rxsig(:,1)))
title("Range profile first pulse")


%range compression
pulseCompression = phased.RangeResponse( ...
    'RangeMethod','Matched filter',...
    'PropagationSpeed',c,...
    'SampleRate',fs);

matchingCoeff = getMatchedFilter(waveform);

[cdata,rnggrid] = pulseCompression(rxsig,matchingCoeff);
figure
imagesc(abs(cdata))
colorbar

figure
imagesc(abs(cdata))
title("Range Compressed Data")
xlabel("Pulse")
ylabel("Range")

%backprojection SAR image

Nx = 100;
Ny = 100;

x = linspace(-50,50,Nx);
y = linspace(50,150,Ny);

img = zeros(Ny,Nx);

for ix = 1:Nx
for iy = 1:Ny

    pixel = [x(ix) y(iy) 0];

    for n = 1:Np

        R = norm(platformPos(n,:) - pixel);

        tau = 2*R/c;

        bin = round(tau*fs);

        if bin>0 && bin<=Nr

            img(iy,ix) = img(iy,ix) + ...
                cdata(bin,n)*exp(1j*4*pi*R/lambda);

        end

    end

end
end

figure
imagesc(x,y,abs(img))
axis xy
colormap hot
title("Backprojection SAR Image")

%platform velocity
dt = 1/prf;

vx = diff(platformPos(:,1))/dt;
vy = diff(platformPos(:,2))/dt;

platformSpeed = mean(sqrt(vx.^2 + vy.^2));

%scene/target center
targetLat = mean(latTargets(:));
targetLon = mean(lonTargets(:));

[xc,yc,zc] = geodetic2enu(targetLat,targetLon,0,lat0,lon0,h0,wgs84Ellipsoid);

Rc = sqrt(xc^2 + yc^2 + zc^2);
disp("max cdata")
disp(max(abs(cdata(:))))
%range Migration Algorithm

% rma_processed = helperSquintRangeMigration( ...
%     cdata,fastTime,fc,fs,prf,...
%     platformSpeed,Np,c,Rc,0);
% 
% figure
% imagesc(abs(rma_processed.'))
% title("RMA Focused Image")
% xlabel("Azimuth")
% ylabel("Range")
% colorbar