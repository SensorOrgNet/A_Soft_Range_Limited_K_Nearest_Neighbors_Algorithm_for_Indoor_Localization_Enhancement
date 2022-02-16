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

function Position = SRL_KNN_Mean_Rank(OldPosition, sigma, Mean_Current_RSSI, Rank_Current_RSSI, ...
                    Cluster_Vector, Cluster_Rank, Length_Cluster, Num_Mac,Num_Rank)

    % Parameter:
    Num_Min_Selection = 5;
    Num_Nearest_Neighbours = 3;
    
    Distance_Array = zeros(1,Length_Cluster);
%% Using the Fingerprint Mean
    % Calculate the distance for Cluster
    for ii = 1:Length_Cluster
        Count_Distance = 0;
        Pos_X_Temp = Cluster_Vector(ii,1); % X
        Pos_Y_Temp = Cluster_Vector(ii,2); % Y
        Pos_X = OldPosition(1);
        Pos_Y = OldPosition(2); 
        Distance_Temp = (Pos_X_Temp-Pos_X)^2 + (Pos_Y_Temp-Pos_Y)^2;
        % Calculate Weight
        Weight_Point = exp(Distance_Temp/(2*sigma^2)); % Based on Gaussian
    
        % Scan all MAC
        for jj=1:Num_Mac
            RSSI_Temp = Cluster_Vector(ii,2+jj); % RSSI position in cluster: 3,4...
            if (RSSI_Temp > -99) && (Mean_Current_RSSI(jj)>-99) % RSSI is strong enough
                Distance_Array(ii) = Distance_Array(ii)+(Mean_Current_RSSI(jj)-RSSI_Temp)^2;
                Count_Distance = Count_Distance + 1;
            end
        end  
        Distance_Array(ii) = Weight_Point*(Distance_Array(ii)/Count_Distance); % Calculate Distance
    end
    
%% Using the Fingerprint RSSI Difference 
   % Sort Distance Array:
    Distance_Array_Sort = sort(Distance_Array);
    
    % Take 7 Min values from the Array
    Min_Array_Pos = Distance_Array_Sort(1,Num_Min_Selection);
    for ii=1:Num_Min_Selection
        for jj = 1:Length_Cluster 
            if Distance_Array_Sort(ii) == Distance_Array(jj)
                Min_Array_Pos(ii) = jj;
            end
        end
    end
    
    % Choose 3 most matched points using RSSI differences
    Distance_Rank_Array = zeros(1,Num_Min_Selection);
    for ii=1:Num_Min_Selection
        Pos_X_Temp = Cluster_Vector(Min_Array_Pos(ii),1); % X
        Pos_Y_Temp = Cluster_Vector(Min_Array_Pos(ii),2); % Y
        Pos_X = OldPosition(1);
        Pos_Y = OldPosition(2); 
        Distance_Temp = (Pos_X_Temp-Pos_X)^2 + (Pos_Y_Temp-Pos_Y)^2;
        % Calculate Weight
%         Weight_Point = exp(Distance_Temp/(2*sigma^2)); % Based on Gaussian
    
        % Scan all MAC
        % Euclidean
         % Scan all MAC
        for jj=1:Num_Rank
            Rank_Temp = Cluster_Rank(Min_Array_Pos(ii),2+jj); % RSSI position in cluster: 3,4...
            Distance_Rank_Array(ii) = Distance_Rank_Array(ii)+(Rank_Current_RSSI(jj)-Rank_Temp)^2;
        end  
      end
    
     % Sort Distance Array:
    Distance_Rank_Array_Sort = sort(Distance_Rank_Array);
    
   % Take Num_Nearest_Neighbours Min values from the Array
    Min_Rank_Array_Pos = zeros(1,Num_Nearest_Neighbours)-1;
    for ii=1:Num_Nearest_Neighbours
        for jj = 1:Num_Min_Selection 
            if (Distance_Rank_Array_Sort(ii) == Distance_Rank_Array(jj)) && (Min_Array_Pos(jj) ~= Min_Rank_Array_Pos(1)) ...
                          &&(Min_Array_Pos(jj) ~= Min_Rank_Array_Pos(2)) && (Min_Array_Pos(jj) ~= Min_Rank_Array_Pos(3))
                Min_Rank_Array_Pos(ii) = Min_Array_Pos(jj);
                break;
            end
        end
    end
    
    Min_Position = Min_Rank_Array_Pos(1);
    Second_Min_Position = Min_Rank_Array_Pos(2);
    Third_Min_Position = Min_Rank_Array_Pos(3);
    
    X1 = Cluster_Vector(Min_Position,1); % X
    Y1 = Cluster_Vector(Min_Position,2); % Y

     % Calculate Combine Probability: beta = (RSSI-RSSI1)^2+...+(RSSI-RSSI6)^2
    beta1 = 0; 
    Count = 0;
    for jj=1:Num_Mac
        RSSI_Temp = Cluster_Vector(Min_Position,2+jj); % RSSI position in cluster: 5,7,9 ...
        if RSSI_Temp > -99 % RSSI is strong enough
            beta1 = beta1+(Mean_Current_RSSI(jj)-RSSI_Temp)^2;
            Count = Count + 1;
        end
    end  
    beta1 = 1/(beta1/Count+1); 
    
    % K = 3   
    X2 = Cluster_Vector(Second_Min_Position,1); % X
    Y2 = Cluster_Vector(Second_Min_Position,2); % Y
     % Calculate Combine Probability: beta = (RSSI-RSSI1)^2+...+(RSSI-RSSI6)^2
    beta2 = 0; 
    Count = 0;
    for jj=1:Num_Mac
        RSSI_Temp = Cluster_Vector(Second_Min_Position,2+jj); % RSSI position in cluster: 5,7,9 ...
        if RSSI_Temp > -99 % RSSI is strong enough
            beta2 = beta2+(Mean_Current_RSSI(jj)-RSSI_Temp)^2;
            Count = Count + 1;
        end
    end  
    beta2 = 1/(beta2/Count+1); 
    
    X3 = Cluster_Vector(Third_Min_Position,1); % X
    Y3 = Cluster_Vector(Third_Min_Position,2); % Y
         % Calculate Combine Probability: beta = (RSSI-RSSI1)^2+...+(RSSI-RSSI6)^2
    beta3 = 0; 
    Count = 0;
    for jj=1:Num_Mac
         RSSI_Temp = Cluster_Vector(Third_Min_Position,2+jj); % RSSI position in cluster: 5,7,9 ...
        if RSSI_Temp > -99 % RSSI is strong enough
            beta3 = beta3+(Mean_Current_RSSI(jj)-RSSI_Temp)^2;
            Count = Count + 1;
        end
    end  
    beta3 = 1/(beta3/Count+1);  
   
    beta_total = beta1+beta2+beta3;

    % Calculate the Combination Probability
    % beta_i = sum((Current_RSSI-RSSIi)^2)
    % X = X1*(beta1/beta_total)+ X2*(beta2/beta_total)+ X3*(beta3/beta_total)
    
    Position(1) =  X1*(beta1/beta_total) + X2*(beta2/beta_total)+ X3*(beta3/beta_total); % X
    Position(2) =  Y1*(beta1/beta_total) + Y2*(beta2/beta_total)+ Y3*(beta3/beta_total); % X

   
end
        
        

