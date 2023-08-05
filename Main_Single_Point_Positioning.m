clc
clearvars
close all
format longg

%=========================================================================%
%=====                                                               =====%
%=====          Single Point Positioning with GPS                    =====%
%=====                                                               =====%
%=====        MSR-04 Methods for quality assurance                   =====%
%=====        Task A and B - Jeanne Paulie Yap                       =====%
%=========================================================================%

%% ----- Definition of constant and directory--------------------------- %%
% speed of light (global variable)
global v_light; 
v_light = 299792458;  % [m/s]

% working directory
path = pwd;
% add functions
addpath([path,'/02_functions']);

%% ----- Import data --------------------------------------------------- %%
% navigation file
nav = '01_data/master_igg.16n';
% observation RINEX file
rnx = '01_data/master_igg_short.16o';
[errNum, header] = readRinexHeader('01_data/master_igg_short.16o');
% load *.mat files if they already exist 
if isfile('01_data/eph.mat') && isfile('01_data/observations.mat') && ...
        isfile('01_data/satellites.mat') && isfile('01_data/time.mat')
    
    load('01_data/eph.mat');
    load('01_data/observations.mat');
    load('01_data/time.mat');
    load('01_data/satellites.mat');

else
    % ephemeris
    [~, eph, ephTime, ~] = readEphGps(nav,'prepareData',40);
    % observations
    [~, obs, date, time] = readRinex(rnx,'prepareData',40,1);
    %read header
    % [errNum, header] = readRinexHeader(rnx,'prepareData',40,1);
    % disp(header)


    % reallocate dataset
    n = size(date,1);
    observations = cell(n,1);
    satellites = cell(n,1);
    for i = 1:n
        observations{i} = obs(i).obsG(obs(i).satG,1);
        satellites{i} = obs(i).satG;
    end
    % save mat files
    % save('01_data/eph.mat','eph')
    % save('01_data/time.mat','time')
    % save('01_data/satellites.mat','satellites')
    % save('01_data/observations.mat','observations')
end


% check if 'posCode' field exists in'header'
if isfield(header, 'posCode')
        rec_pos = header.posCode
else
        error("The 'posCode' field does not exist in the 'header' structure.");
end

    % disp(recCor)

%% ----- Single Point Positioning -------------------------------------- %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO: initialize vector/matrix for parameter and results
position_ = zeros(length(time),3);
PDOP_ = zeros(length(time),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rec_pos = header.posCode;
% % ----- positioning
    % az = cell(length(time),1);
    % el = cell(length(time),1);
    % dist = cell(length(time),1);
    % SatPosition = cell(length(time),1);

  
for i = 1:length(time)

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% TODO: write/complete the function 'calcSPP' to calculate the receiver
%     %%       position based on measurements of a single epoch at time(i)

    fprintf("Observations:");
    celldisp(observations(i,:));
    %disp(observations(i,:))
    [num_of_sats, position,PDOP,design_mat_H,cov_mat_W] = RANSAC(eph, time(i,:), satellites(i,:), observations(i,:),rec_pos');
    
    % fprintf("Design Matrix:");
    % disp(design_mat_H);
    % celldisp(satellites(i,:))
    % 
    % fprintf("Covariance Matrix:");
    % disp(cov_mat_W);
    % 
    % fprintf("Number of Satellites:");
    % disp(num_of_sats);

    num_of_sats_(i,:) = num_of_sats %task A


    % fprintf("Num of sats (cell)");
    % celldisp(num_of_sats_(i,:));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %% TODO: save output of 'calcSPP' in the initialized result vector/matrix
    position_(i,:) = position;

    disp(position_);
    % celldisp(position_(i,:));
    PDOP_(i,:) = PDOP;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%comment out:     
% [num_of_sats,position,PDOP, azimuth, elevation , disti, Sat_Pos] = calcSPP(eph, time(i,:), satellites(i,:), observations(i,:),rec_pos');

    % az{i} = azimuth;
    % el{i} = elevation;
    % dist{i} = disti;
    % SatPosition{i} = Sat_Pos;


%% ----- Plot results -------------------------------------------------- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO: plot results in cartesian coordinates with corresponding time

% figure(1)
% x = (time - time(1,1))/3600;
% X_ = position_(:,1) - mean(position_(:,1));
% Y_ = position_(:,2) - mean(position_(:,2));
% Z_ = position_(:,3) - mean(position_(:,3));
% plot(x,X_,'r',x,Y_,'b',x,Z_,'g')
% legend('X','Y','Z')
% title('Coordinates')
% ylabel('Position (m)')
% xlabel('Time(h)')

%% TODO: plot PDOP values and satellites with corresponding time

% figure(2)
% x = (time - time(1,1))/3600;
% plot(x,PDOP_,'r',x,num_of_sats_,'b')
% legend('PDOP','No. of Satellites')
% title('PDOP Values and Satellites Count')
% ylabel('Value')
% xlabel('Time(h)')
% 
% figure(3)
% x = (time - time(1,1))/3600;
% plot(x,PDOP_)
% title('PDOP Values')
% ylabel('PDOP')
% xlabel('Time(h)')
% 
% figure(4)
% x = (time - time(1,1))/3600;
% plot(x,num_of_sats_)
% title('Satellites Count')
% ylabel('Observed Satellites')
% xlabel('Time(h)')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part IV - Analysis of SPP Results 

% 1. Plot the calculated coordinates minus the mean as a time series separated by X-, Y- and Zcoordinates depending on the observation epochs. 
%     Seen on Figure 1

% 2. Implement different elevation masks for the single point positioning: Is the application of elevation masks recommended in practice?
%     Seen on Line 82 of calcSPP 
%     Elevation masks are still recommended in practice since they can
%     filter out satellite signal that are deemed obstructed or
%     unreliable.This depends on the application and environment though. 


% 3. Plot the number of satellites and corresponding PDOPs. Are there any correlations between PDOPs and corresponding position?
%     Seen on Figures 2 to 4. Yes, there is a correlation between the PDOPs
%     and corresponding position. The PDOP can provide information about
%     the precision and accuracy of the position solution. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TASK B - P RANSAC Application 
% for i = 1:length(time)
% 
% %missing: initialize best_model, inliers_idx 
% %missing: save outputs of the RANSAC function to the initializaed matrices!
% 
% 
% %RANSAC testing 
%     num_iterations = 5;
%     threshold = 2.0;
%     [best_model, inliers_idx] = RANSAC(design_mat_H, position_(i,:), cov_mat_W, num_iterations, threshold)
% 
% end
