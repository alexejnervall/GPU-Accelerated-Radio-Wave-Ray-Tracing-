clear; clc; close all;

%constants
c = physconst('LightSpeed');
fc = 10e9;
lambda = c/fc;

%scene viewer
viewer = siteviewer(Buildings="kista.osm");

%radar parameters
Np = 20;           % pulses
aperture = 80;      % meters

%reference location (Kista)
lat0 = 59.4039;
lon0 = 17.9448;
h0 = 30;

meters_per_deg_lon = 111320*cosd(lat0);

deltaLon = (aperture/2)/meters_per_deg_lon;

dLon = linspace(-deltaLon,deltaLon,Np);

lat = lat0*ones(1,Np);
lon = lon0 + dLon;
h   = h0*ones(1,Np);

%target grid
Nx = 10;
Ny = 10;

targetLat = 59.4034;
targetLon = 17.9457;

latGrid = linspace(targetLat-0.00009,targetLat+0.00009,Nx);
lonGrid = linspace(targetLon-0.00009,targetLon+0.00009,Ny);

[latTargets,lonTargets] = meshgrid(latGrid,lonGrid);

numTargets = numel(latTargets);

%propagation model
pm = propagationModel("raytracing", ...
    Method="sbr", ...
    MaxNumReflections=2, ...
    MaxNumDiffractions=0);

%storage
platformPos = zeros(Np,3);
rayData = struct;
totalRays = 0;
%ray tracing loop
for n = 1:Np
    fprintf("Pulse %d / %d\n", n, Np);
    tx = txsite( ...
        Latitude = lat(n), ...
        Longitude = lon(n), ...
        AntennaHeight = h(n), ...
        TransmitterFrequency = fc);

    for t = 1:numTargets
        rayData(Np,numTargets).R = [];

        if isfield(rayData(n,t),'R')
            totalRays = totalRays + length(rayData(n,t).R);
        end
        

        tgt = rxsite( ...
            Latitude = latTargets(t), ...
            Longitude = lonTargets(t), ...
            AntennaHeight = 0);

        rays = raytrace(tx,tgt,pm,Type="pathloss");
        if n==1 && t==1
            disp("Example ranges:")
            disp([rays{1}.PropagationDistance])
        end

        if isempty(rays{1})
            continue
        end

        for k = 1:length(rays{1})

            pathLossDB = rays{1}(k).PathLoss;
            R = rays{1}(k).PropagationDistance;

            amplitude = 10^(-pathLossDB/20);
            phase = exp(-1j*4*pi*R/lambda);

            rayData(n,t).R(k) = R;
            rayData(n,t).amplitude(k) = amplitude;
            rayData(n,t).phase(k) = phase;

        end

    end

    %convert UAV position to ENU
    [x,y,z] = geodetic2enu(lat(n),lon(n),h(n),lat0,lon0,h0,wgs84Ellipsoid);
    platformPos(n,:) = [x y z];

end
disp("Total rays found:")
disp(totalRays)

%plot UAV path
figure
plot(platformPos(:,1),platformPos(:,2),'o-')
xlabel('x (m)')
ylabel('y (m)')
title('UAV Trajectory')
axis equal

%save data
save("raytraceSARdata.mat")

disp("Ray tracing complete")











