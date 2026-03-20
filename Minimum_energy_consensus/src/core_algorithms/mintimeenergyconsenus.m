%function [xbar, tfbar] = mintimeenergyconsenus(IC,beta)

%% === INPUT PARAMETERS ===

%IC = [-2609.81,54.676; 1768.64,30.106; 996.061,-22.228; 0,20; 2000,0; 4000,-30]%---CORRECT RESULT X =[2000,10]
%IC = [-3376.3,73.002; 299.312,60.102; -1201.54,20.212; 0,40; 2000,0; 3000,-20]%---CORRECT RESULT X=[3118,29]
%IC = [0.04, 0.1; 0.39, 1.05; 0.3, -0.525; 0.5, -0.525; 0.2, -0.2 ; 0.1, 0.1];

%IC = [0.04, 0.1; 0.39, 1.05; 0.3, -0.525; 0.5, -0.525];%--correct result X=[3118,29]
%IC = [0.04, 0.1; 0.39, 1.05; -0.9, 0.72];

beta =30; 
IC = [3,-1.1; 
      -3, 1.1;  
      3.1, -3;
      -3.1 ,3; 
       6.5, 6.5 ;
       -7.5,-6.5];
% tf = 8.158576;
%xc = -0.5 ;  
% yc = 0;

tolerance = 1e-7;
numAgents = size(IC, 1);
numCases = 4;

% === GENERATE ALL TRIPLETS ===-
triplet_inds = nchoosek(1:numAgents, 3);
numTriplets = size(triplet_inds, 1);

% === INITIALIZE STORAGE ===
tf_matrix = nan(numCases, numTriplets);  % 18 x 20
flag_matrix = zeros(numCases, numTriplets);
sw_matrix = nan(numTriplets * numCases, 12);
xc_yc_matrix = nan(numTriplets * numCases, 2);

% === MAIN PROCESS LOOP ===
for tripIdx = 1:numTriplets
    idxs = triplet_inds(tripIdx, :);
    IC_trip = IC(idxs, :);
    flag = 0
    c = 1; %case number
    while flag == 0
        fname = sprintf('case%d', c);
        x110 = IC_trip(1, 1); x120 = IC_trip(1, 2);
        x210 = IC_trip(2, 1); x220 = IC_trip(2, 2);
        x310 = IC_trip(3, 1); x320 = IC_trip(3, 2);
        try
            [xc, yc, tf, flag] = ...
                feval(fname, x110, x120, x210, x220, x310, x320, beta);

            idx = (tripIdx - 1) * numCases + c;  % keep for sw/xc_yc
              % write by (row=case, col=triplet)
            if flag == 0
                tf_matrix(c, tripIdx) = NaN;
            else 
                tf_matrix(c, tripIdx) = tf;
            end
            xc_yc_matrix(idx, :) = [xc, yc];
        catch ME
            warning('Error in %s, triplet %d, case %d: %s', fname, tripIdx, c, ME.message);
        end
        c = c+1
    end
end

tf_matrix_reshaped = reshape(tf_matrix, numCases, numTriplets);
flag_matrix_reshaped = reshape(flag_matrix, numCases, numTriplets);

%% === COMPUTE MINIMUM BETAS PER TRIPLET ===
min_tf_per_triplet = min(tf_matrix_reshaped, [], 1, 'omitnan'); % 1 x numTriplets vector

fprintf('Minimum tf per triplet:\n');
disp(min_tf_per_triplet);

%% === CREATE MATRIX FOR MEAN SW1 (start) AND SW2 (end) ===

%% === FIND MAXIMUM OF THE MINIMUM BETAS AND CORRESPONDING SWITCHING TIMES AND STATES ===
[max_of_min_tf, idx_max] = max(min_tf_per_triplet);

[~, min_case_idx] = min(tf_matrix_reshaped(:, idx_max), [], 'omitnan');

xc_yc_matrix_reshaped = reshape(xc_yc_matrix, numCases, numTriplets, 2);
states_for_max_min_tf = squeeze(xc_yc_matrix_reshaped(min_case_idx, idx_max, :));

fprintf('Maximum of minimum tf: %.6f found in triplet %d, case %d\n', max_of_min_tf, idx_max, min_case_idx);
fprintf('Corresponding [xc, yc] states:\n');
disp(states_for_max_min_tf');

tfbar = max_of_min_tf;
xbar = states_for_max_min_tf;
%end