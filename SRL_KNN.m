%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K-Nearest Neighbour (KNN)
% Input : 
% - Mean_Current_RSSI : Array of Mean RSSI for each APs 
% - Cluster_Array:  X Y Z MAC1 Mean_RSSI1 MAC2 Mean_RSSI2 ...   
% - Num_Mac (Number of APs in each vector)
% - Value of K (1,2,3)
% Output: 
% X Y Z

% Euclidean Distance
% D = sqrt((ss-ss1)^2+(ss-ss1)^2+...)
% Weight:
% Wi = exp(|X-Xi-1|^2/(2*sigma^2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Position = SRL_KNN(OldPosition, sigma, Current_RSSI, ...
                    Cluster_Vector,Length_Cluster, Num_Mac, K)
    
    Position = zeros(1,2);
    Distance_Array = zeros(1,Length_Cluster); 
    Distance_Array_NoWeight = Distance_Array;
    
    Weight_Sum = 0;
    % Calculate the distance for Cluster
    for ii = 1:Length_Cluster
        Count_Distance = 0;
        % Get X,Y of RPs 
        Pos_X_Temp = Cluster_Vector(ii,1); % X
        Pos_Y_Temp = Cluster_Vector(ii,2); % Y
        % Get X,Y of previous position
        Pos_X = OldPosition(1);
        Pos_Y = OldPosition(2); 
        % Distance between 2 locations
        Distance_Temp = (Pos_X_Temp-Pos_X)^2 + (Pos_Y_Temp-Pos_Y)^2;
        
        % Calculate the Weight based on the Gaussian Window
        Weight_Point = exp(Distance_Temp/(2*sigma^2)); % Based on Gaussian
        Weight_Sum = Weight_Sum + Weight_Point;
        
        % Scan all MAC
        for jj=1:Num_Mac
            RSSI_Temp = Cluster_Vector(ii,2+jj); % RSSI position in cluster: 3,4...
            if Current_RSSI(jj)>=-99 % RSSI is strong enough
                Distance_Array(ii) = Distance_Array(ii)+(Current_RSSI(jj)-RSSI_Temp)^2;
                Count_Distance = Count_Distance + 1;
            end
        end  
        % No Apply Weight
        Distance_Array_NoWeight(ii) = (Distance_Array(ii)/Count_Distance); % Calculate Distance 
        % Apply Weight
        Distance_Array(ii) = Weight_Point*(Distance_Array(ii)/Count_Distance); % Calculate Distance      
    end
    Distance_Array = Distance_Array/Weight_Sum; 
    % Calculate Mean of distance
   % Sort Distance Array:
    [~, Idx_sort] = sort(Distance_Array);
%     [Distance_Array_Sort2, index2] = sort(Distance_Array_NoWeight);
   
   % Get the weighted average position
   beta_total = 0;
   for ii=1:K
        MinPos = Idx_sort(ii);
        X_temp = Cluster_Vector(MinPos,1); % X
        Y_temp = Cluster_Vector(MinPos,2); % Y
        beta = Distance_Array(MinPos);   
        beta = 1/beta;
        beta_total = beta_total + beta;
        Position(1) =  Position(1) + X_temp*beta; % X
        Position(2) =  Position(2) + Y_temp*beta; % Y
   end
   Position(1) =  Position(1)/beta_total; % X
   Position(2) =  Position(2)/beta_total;
   
end
        
        

