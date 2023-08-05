function [num_of_sats,rec_pos, PDOP, design_mat_H,cov_mat_W] = RANSAC(eph, time, satellites, observations,rec_pos)
% function [num_of_sats,rec_pos, PDOP, azimuth , elevation , disti, Sat_Pos] = calcSPP(eph, time, satellites, observations,rec_pos)
% ----------------------------------------------------------------------- %
% function for calculation of single point positioning for a single (time)epoch
% ----------------------------------------------------------------------- %

% ----- Single point positioning -------------------------------------- %%
% Initialize variables
iteration = 1;
deltax = 9999;
c = 299792458;
sv = cell2mat(satellites);
obs = cell2mat(observations);

% TODO: initialize vector/matrix for parameter and results
G = zeros(length(sv), 4);
dist_ = [];
eleva_ = [];
B = 0;
count = 1;

%insert subset selection and RANSAC here... 

while max(abs(deltax)) > 1e-3
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    out_aspects = [];
    % For-each observation
    m = length(sv);

    for i = 1:m
        if iteration<2
        t_travel = obs(i,1)/c;
        else
            t_travel = dist_(i,1)/c;
        end
   
        %% (1) Satellite Position and Clock for Transmission Time
        % Correction of transmission time with calculated signal travel time based on code observation or calculated distance
        % Calculate satellite position and clock offset

        transmTime = time - t_travel;
        satellites_ = sv(i,1);
        [Sat_Pos, Sat_Clk] = satellitePosition(transmTime, eph,satellites_ );


        %% (2) Correction due to earth rotation
        [Sat_Pos_Rot] = e_r_corr(t_travel,Sat_Pos);
       

        %% (3) Impact of Troposphere 
        [az,el,dist] = topocent(rec_pos(1:3,:),Sat_Pos_Rot-rec_pos(1:3,:));    
        [troposhere] = tropo(sin(el),0,1013,293,50,0,0,0);
        dist_(i,1) = dist;
        eleva_(i,1) = rad2deg(el);
        

        %% (4) Estimate Receiver - Satellite Distance
        fxyz = norm(Sat_Pos_Rot-rec_pos, 'fro');   %geometric distance between sat_pos_rot and re_pos_0 in m
        C = Sat_Clk * c;  %satellite clock offset in m
        D = troposhere;   %tropo error in m
        sat_dist(i,1) = obs(i,1) - (fxyz + B - C + D); %sat_dist = measured - calc

        %% (5) LSQ Adjustment - Gauss Markov Model
        %% TODO: write the calculation of designmatrix A (G) , 
        %% cofactor matrix Q_ll and delta_l 

        G(i, :) =  [-(Sat_Pos_Rot(1) - rec_pos(1)) / fxyz ,...
                    -(Sat_Pos_Rot(2) - rec_pos(2)) / fxyz ,...
                    -(Sat_Pos_Rot(3) - rec_pos(3)) / fxyz ,...
                     1];                 

        P = eye(i,i);
        Sat_Pos(i,1:3) = Sat_Pos_Rot';
        azimuth(i,1) = az;
        elevation(i,1) = el;
        disti(i,1) = dist;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% (6) Additional Aspects- perform a satellite selection based on elevation mask
        %%       - calculate delta_x in a least square adjustment 
        
        %% Integration of elevation mask 
        for k = 1:length(sv)
        if eleva_(k) < 5 
            out_aspects(count) = k;
            count = count +1;
        end
        end
        
        if isempty(out_aspects) == 0
            G(out_aspects,:) = [];
            P = eye(size(G,1),size(G,1));
            sat_dist(out_aspects,:) = [];
        end

        deltax = pinv(G'*P*G)*G'*P*sat_dist;

        rec_pos =  rec_pos + deltax(1:3,:);


        %Insert inlier vs outlier decision making here 
        %iterate again
        B = B + deltax(4,:);        

        %% update position for next iteration 
        iteration = iteration + 1; 
        count = 1;
end
fprintf("Design Matrix:");
disp(G);

design_mat_H = G;

num_of_sats = size(G,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (6) calculate PDOP

 covariancw = pinv(G'*P*G);
 % fprintf("covariance mat:");
 % disp(covariancw);
 cov_mat_W = covariancw;
 PDOP = sqrt(covariancw(1,1) + covariancw(2,2)+ covariancw(3,3));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

