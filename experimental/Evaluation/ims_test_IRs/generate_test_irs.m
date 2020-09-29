clear all, close all %#ok
% requires:
% addpath('../../../resources/shoebox-roomsim')

% PATH SETUP 
path_out = './';
if ~exist(path_out, 'dir')
    mkdir(path_out)
end

% GLOBALS
fs = 48000;
c = 343;
sh_order = 20;
roomType = 'medium'; 
APPLY_AIR_ABSORPTION = 1;


%% SCENE PARAMETERS
% Here set up you scene parameters, source locations, microphone locations,
% room boundaries and target reverberation times (or directly wall
% absorption coefficients). If there is no room (anechoic rendering) the
% global origin is arbitrary (e.g. can be at one of the microphones),
% however if there is a room (reverberant rendering), all receiver and
% source positions should be given with respect to the bottom left corner
% of the room (top view), with positive x+ extending to the east, and
% positive y+ extending to the north, while z+ is extending purpendicular
% to them towards the viewer (right-hand rule)
%
%   length/width
%   |----------|
%   ^ y           .
%   |    ^z      /height
%   |   /       /
%   .__/_______.           _
%   | /        |           |
%   |/         |           | width/length
%   o__________.------> x  _
%
% Note that there is no checking for the source-microphone coordinates
% falling inside the boundaries of the room.

%%%%%%%%%%% ROOM 
switch roomType
    case 'medium'
        room = [10 7 3];
        rt60 = [1.0 0.8 0.7 0.6 0.5 0.4].*0.666;
        abs_wall_ratios = [0.75 0.86 0.56 0.95 0.88 1];
        
    case 'anechoic'
        room = [10 7 3];
        rt60 = 0.1;
        abs_wall_ratios = 1;
        % NOTE: abs_wall is forced to be 1 below
end
        
% source position 
clear src 
src(1,:) = [5.1 6 1.1];  

% receiver position
clear rec 
rec(1,:) = [5 4 1]; 
nRec = size(rec, 1);

% reverberation time per octave band 
nBands = length(rt60);

% lowest octave band
band_centerfreqs(1) = 125;
for nb=2:nBands, band_centerfreqs(nb) = 2*band_centerfreqs(nb-1); end % octave band centerfreqs for RT60
% absorption for approximately achieving the RT60 above - row per band
abs_wall = findAbsCoeffsFromRT(room, rt60,abs_wall_ratios);
% critical distance for the room
[~,  d_critical] = room_stats(room, abs_wall);
if strcmp(roomType, 'anechoic'), abs_wall = ones(size(abs_wall)); end

% store scene parameters
scene.room = room;
scene.src = src;
scene.rec = rec;
scene.rt60 = rt60;
scene.abs_wall = abs_wall;
scene.abs_bands = band_centerfreqs;
scene.abs_wall_ratios = abs_wall_ratios;
scene.AIR_ABS = APPLY_AIR_ABSORPTION;
scene.fs = fs;
scene.c = c;
scene.path_out = path_out;

 
%% Generate room IRs 
% limit the RIR by reflection order or by time-limit
nBands = length(scene.rt60);
type = 'maxTime'; % 'maxTime' 'maxOrder'

maxlim = max(scene.rt60); % just cut if it's longer than that ( or set to max(rt60) )
for nb = 1:nBands
    if (scene.rt60(nb)<maxlim) limits(nb,1) = scene.rt60(nb);
    else limits(nb,1) = maxlim;
    end
end
         
% compute echograms
nCH = (sh_order+1)^2;
src2 = [scene.src(:,1) scene.room(2)-scene.src(:,2) scene.src(:,3)]; % change y coord for src/rec due to convention inside the IMS function
rec2 = [scene.rec(:,1) scene.room(2)-scene.rec(:,2) scene.rec(:,3)]; % change y coord for src/rec due to convention inside the IMS function
abs_echograms = compute_echograms_sh(scene.room, src2, rec2, scene.abs_wall, limits, sh_order);
%abs_echograms_array = compute_echograms_arrays(scene.room, src2, rec2, scene.abs_wall, limits);
% time   - reflection propagation time
% value  - reflection propagation attenuation
% order  - reflection order
% coords - image source coordinates with respect to receiver
 
% render RIRs
sh_rirs = render_sh_rirs(abs_echograms, scene.abs_bands, scene.fs); 
sh_rirs = sh_rirs*sqrt(4*pi);
nSrc = size(src2,1); 

disp('Applying air absorption - SH')
if scene.AIR_ABS
    for ns=1:nSrc 
        for nr=1:nRec
            for nc=1:nCH
                sh_rirs_air(:,nc,nr,ns) = applyAirAbsorption(sh_rirs(:,nc,nr,ns),scene.fs);
            end 
        end
    end
    sh_rirs = sh_rirs_air;
end


%% Save 
audiowrite(['ref_o' num2str(sh_order) '_' roomType '_room.wav'], sh_rirs(:,:,1,1), fs);
 

%% Plot  
figure
ni=1;
for ns = 1:nSrc
    subplot(1,nSrc,ni)
    plot(sh_rirs(:,1,1,ns)) 
    title(['pressure, source: ' num2str(ns)])
    ni=ni+1;
end



