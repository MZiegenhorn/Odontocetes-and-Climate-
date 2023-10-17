%created by MAZ on 10/11/22 to plot Hawaii maps using MATLAB version of
%Generic Mapping Tools software
%adapted from code by SBP 

clear all
close all
clc

addpath(genpath('E:\code'))
addpath(genpath('C:\programs\gmt6'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
site = 'Hawaii'; %set site

% setup for colormap
mindepth = -5000; %depth in meters
maxheight = 100; %height in meters
step = 10; %step size in meters

%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
depthstr = [num2str(mindepth),'/',num2str(maxheight),'/',num2str(step)];

if strcmp(site,'Manawai')
    %read in gridded data
    G = gmt ('read', '-Tg E:\NCSS_bathymetric.nc'); %file needs to be in current directory
    % Evaluate an asymmetrical color palette with hinge at sealevel
    C = gmt('makecpt', ['-Cgeo -T',depthstr]);
    % C = gmt('makecpt', '-Cbukavu -T-4000/3000/20');
    % Cgray = gmt('makecpt', '-Cgray -T-4000/0/100');

    %set NaNs to zero
    zvals = G.z;
    zvals(isnan(zvals)) = 100;
    %set above sea level to 100
    zvals(zvals>=0) = 100;
    G.z = zvals;

    %change colormap such that land is all grey
    gr = [0.5 0.5 0.5];
    idx = abs(mindepth)/step-6;
    C.colormap(idx+1:end,:) = repmat(gr,length(C.colormap(idx+1:end,1)),1);
    C.cpt(idx+1:end,:) = repmat(gr,length(C.colormap(idx+1:end,1)),2);

    % to create a png for matlab display
    clear map
    % Make the GMT plot
    % map = gmt('psbasemap', '-B50m50m/50ma50m -JM8i -R-122/-116.87/32.22/34.75 -K -X1 -Y1');
    map = gmt('psbasemap', '-B10m10m/10ma10m -JM8i -R-176.2/-175.5/27.64/28.1 -K -X2 -Y1');

elseif strcmp(site,'Hawaii')
    %read in gridded data
    G = gmt ('read', '-Tg E:\UHI_bathymetric.nc'); %file needs to be in current directory
    % Evaluate an asymmetrical color palette with hinge at sealevel
    C = gmt('makecpt', ['-Cgeo -T',depthstr]);
    % C = gmt('makecpt', '-Cbukavu -T-4000/3000/20');
    % Cgray = gmt('makecpt', '-Cgray -T-4000/0/100');

    %set NaNs to zero
    zvals = G.z;
    zvals(isnan(zvals)) = 100;
    %set above sea level to 100
    zvals(zvals>=0) = 100;
    G.z = zvals;

    %change colormap such that land is all grey
    gr = [0.5 0.5 0.5];
    idx = abs(mindepth)/step-10;
    C.colormap(idx+1:end,:) = repmat(gr,length(C.colormap(idx+1:end,1)),1);
    C.cpt(idx+1:end,:) = repmat(gr,length(C.colormap(idx+1:end,1)),2);

    % to create a png for matlab display
    clear map
    % Make the GMT plot
    % map = gmt('psbasemap', '-B50m50m/50ma50m -JM8i -R-122/-116.87/32.22/34.75 -K -X1 -Y1');
    map = gmt('psbasemap', '-B10m10m/10ma10m -JM8i -R-156.3/-155.8/19.35/19.8 -K -X2 -Y1');

elseif strcmp(site,'all')
    %read in gridded data
    G = gmt ('read', '-Tg E:\gebco_2022_sub_ice_topo\gebco_2022_sub_ice_topo.nc'); %file needs to be in current directory
    % Evaluate an asymmetrical color palette with hinge at sealevel
    C = gmt('makecpt', ['-Cgeo -T',depthstr]);
    % C = gmt('makecpt', '-Cbukavu -T-4000/3000/20');
    % Cgray = gmt('makecpt', '-Cgray -T-4000/0/100');

    %set NaNs to zero
    zvals = G.z;
    zvals(isnan(zvals)) = 100;
    %set above sea level to 100
    zvals(zvals>=0) = 100;
    G.z = zvals;

    %change colormap such that land is all grey
    gr = [0.5 0.5 0.5];
    idx = abs(mindepth)/step-6;
    C.colormap(idx+1:end,:) = repmat(gr,length(C.colormap(idx+1:end,1)),1);
    C.cpt(idx+1:end,:) = repmat(gr,length(C.colormap(idx+1:end,1)),2);

    % to create a png for matlab display
    clear map
    % Make the GMT plot
    map = gmt('psbasemap', '-B500m500m/500ma500m -JM8i -R-179/-154/10/30 -K -X2 -Y1');
end

map = gmt('grdimage', '-R -JM -V -X0 -Y0 -O -C -K',G,C);
map = gmt('grdcontour', '-R -C100 -L-8000/0 -Q50 -JM -V -X0 -Y0 -O -K',G);
% map = gmt('pscoast', '-R -J -O -W0.25p -Dh -LjTR+w200k+u+f+c38:30N+o0.5i/0.2i -F+gwhite+p0.5p');
map = gmt('psscale', '-D9i/1i/5c/0.5c -B1000/:m: -C -O -X -Y -V',C);

% Convert plot to a PNG and display in MATLAB as Figure
I = gmt('psconvert -TG -P -E300 -A', map); %transparent PNG

figure()
h = imshow (I.image);
file = ['E:\',site,'_map.png'];
saveas(h,file,'png')
