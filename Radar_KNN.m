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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Position = Radar_KNN(Current_RSSI, Cluster_Vector,Length_Cluster, Num_Mac, K)

    Position = zeros(1,2);
    Distance_Array = zeros(1,Length_Cluster);
    % Calculate the distance
    for ii = 1:Length_Cluster
        
        Count_Distance = 0;
        % Scan all MAC
        for jj=1:Num_Mac
              RSSI_Temp = Cluster_Vector(ii,2+jj); % RSSI position in cluster: 3,4...
              if (RSSI_Temp >= -99) && (Current_RSSI(jj)>=-99) % RSSI is strong enough
                Distance_Array(ii) = Distance_Array(ii)+(Current_RSSI(jj)-RSSI_Temp)^2;
                Count_Distance = Count_Distance + 1;
              end
        end  
        Distance_Array(ii) = (Distance_Array(ii)/Count_Distance); % Calculate Mean 
 %       Mean_Distance_Array = Mean_Distance_Array + Distance_Array(ii);
    end
    [~, Idx_sort] = sort(Distance_Array);
%     [Distance_Array_Sort2, index2] = sort(Distance_Array_NoWeight);
   
   % Get the weighted average position
   for ii=1:K
        MinPos = Idx_sort(ii);
        X_temp = Cluster_Vector(MinPos,1); % X
        Y_temp = Cluster_Vector(MinPos,2); % Y
        Position(1) =  Position(1) + X_temp; % X
        Position(2) =  Position(2) + Y_temp; % Y
   end
   Position(1) =  Position(1)/K; % X
   Position(2) =  Position(2)/K;
end
