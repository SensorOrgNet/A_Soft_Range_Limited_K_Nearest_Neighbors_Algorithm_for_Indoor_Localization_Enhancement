%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Histogram Database
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

function Position = SRL_KNN_Histogram(OldPosition, sigma, Mean_Current_RSSI, ...
                    Cluster_Vector,Length_Cluster, Num_Mac, NumProbRSSI,K)
    
    Position = zeros(1,2);
    Distance_Array = zeros(1,Length_Cluster/NumProbRSSI);
    
    % Find the strongest router  
    % Calculate the distance for Cluster
    CountDis = 1;
    for ii = 1:NumProbRSSI:Length_Cluster-1
        Pos_X_Temp = Cluster_Vector(ii,1); % X
        Pos_Y_Temp = Cluster_Vector(ii,2); % Y
        Pos_X = OldPosition(1);
        Pos_Y = OldPosition(2); 
        Distance_Temp = (Pos_X_Temp-Pos_X)^2 + (Pos_Y_Temp-Pos_Y)^2;
        % Calculate Weight
        Weight_Point = exp(Distance_Temp/(2*sigma^2)); % Based on Gaussian
    
        % Scan all MAC - Calculate Eucliean distance
        for jj=1:Num_Mac
            if Mean_Current_RSSI(jj)>-99 % ignore -100    
                % Scan all Histogram Data
    %             AverageRSSI = 0;
                for CntRSSI = 0:NumProbRSSI-1
                    RSSI_Temp = Cluster_Vector(ii+CntRSSI,3+(jj-1)*2); % RSSI position in cluster: 5,7,9 ...
                    if RSSI_Temp > -99 
                        Prob_Temp = Cluster_Vector(ii+CntRSSI,4+(jj-1)*2);
                        % Modified Euclidean Distance
                        Distance_Array(CountDis) = Distance_Array(CountDis)+Prob_Temp*(Mean_Current_RSSI(jj)-RSSI_Temp)^2;
                    end                   
                end
            end
        end  
        Distance_Array(CountDis) = Weight_Point*(Distance_Array(CountDis)); % Calculate Distance
        CountDis = CountDis + 1;
    end
    
   % Calculate Mean of distance
   % Sort Distance Array:
   % Sort Distance Array:
   [Distance_Array_Sort Idx_sort] = sort(Distance_Array);
   % Distance_Array_Sort = sort(Distance_Array);
    
   if length(Idx_sort) < K % in case not enough K neighbours
        K = length(Idx_sort);
   end
   beta_total = 0;
   for ii=1:K
        MinPos = Idx_sort(ii);
        X_temp = Cluster_Vector((MinPos-1)*NumProbRSSI+1,1); % X
        Y_temp = Cluster_Vector((MinPos-1)*NumProbRSSI+1,2); % Y
        beta = Distance_Array(MinPos);   
        beta = 1/beta;
        beta_total = beta_total + beta;
        Position(1) =  Position(1) + X_temp*beta; % X
        Position(2) =  Position(2) + Y_temp*beta; % Y
   end

   Position(1) =  Position(1)/beta_total; % X
   Position(2) =  Position(2)/beta_total;
   
end
        
        

