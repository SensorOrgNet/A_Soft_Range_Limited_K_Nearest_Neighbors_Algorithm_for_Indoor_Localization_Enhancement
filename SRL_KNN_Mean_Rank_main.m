%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wifi Indoor Localization
% Minhtu 
% Fingerprint Mean & Rank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K-Nearest_Neighbour
% Euclidean Distance
% D = sqrt((ss-ss1)^2+(ss-ss1)^2+...)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Database: Mean RSSI
% X Y Z MAC1 Mean_RSSI1 MAC2 Mean_RSSI2 ...    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter - 1 Unit = 40 inches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = 3; 
Num_Mac = 11; % Number of APs per vector in database
Num_Rank = 11;

% Load Database 
myFolder = 'C:\Users\minh_\Desktop\CSI_RSSI_Database\RSSI_6AP_Experiments\Nexus4_Data\'; % Database Folder
% Database
Input = importdata([myFolder 'MeanDatabase_8June2018_11MAC.csv']);
Database = Input;
Input_Rank = importdata([myFolder 'RankDatabase_15June2018_updated.csv']);
Database_Rank = Input_Rank;

% Test
InputTest = importdata([myFolder 'Long_Traj_172Locations_RSSI.csv']);
TestPoint = InputTest;
Temp = size(TestPoint);
NumTestingPoint = Temp(1); % Number of Test points

% History Buffer Option
IdealFlag = 0; % 1: Using Ideal History
                % 0: Using Real History
       
IdealHistory_Array = TestPoint(:,1:2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Init
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting Position
Starting_Position = TestPoint(1,1:2);

% Assumption
v_user_max = 2; % m/s
t_request = 1; % Period of request time is 5 seconds
distance_max = v_user_max*t_request; % meters

sigma = 2*distance_max; % 1 Unit = 40 inches
Cluster_Max_Distance = 5*distance_max; % Maximum of Possible distance 
Cluster_Min_Distance = 0;

% Estimated Location Buffer
Estimated_Position        = zeros(NumTestingPoint,2);
Distance_Error_Meters     = zeros(1,NumTestingPoint);

% History Buffer
% Num_Buff_His = 2;
History = Starting_Position;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop to find the location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for CountPoint = 1:NumTestingPoint
    
    % RSSI of the current point
    Current_RSSI = TestPoint(CountPoint,3:end);
    
    % Choose 2.4 GHz for ranking
    RSSI_2_4GHz = TestPoint(CountPoint,3:end);
    [Sort Index] = sort(RSSI_2_4GHz,'descend');
    % Rank_Current_RSSI = Index;
    Rank_Current_RSSI = zeros(1,Num_Rank);
    for jj = 1:Num_Rank
       Rank_Current_RSSI(Index(jj)) = jj; 
    end
    
    Core_Weight_Point =  History; % Most Recent Previous Point
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 1: Clustering (Choose Possible Zone)
    % Input: 
    % - Database: X Y Z MAC1 Mean_RSSI1 MAC2 Mean_RSSI2 MAC3 Mean_RSSI3  
    % - Mean_Current_RSSI
    % Output: Cluster includes the vectors which have RSSI values are nearest to Mean_Current_RSSI  
    % - Cluster_Array:  X Y Z MAC1 Mean_RSSI1 MAC2 Mean_RSSI2 MAC3 Mean_RSSI3  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    NumCol_Cluster = 2 + Num_Mac;
    NumCol_Rank = 2 + Num_Rank;
    Size_Database = size(Database);
    Length_Database = Size_Database(1);
    NumRow_Cluster = round(Length_Database/2);
    Cluster_Vector = zeros(NumRow_Cluster,NumCol_Cluster);
    Cluster_Rank = zeros(NumRow_Cluster,NumCol_Rank);
    
    % Possible Zone
    Pos_Array = NearDistanceNeighbour(Core_Weight_Point, Database, Length_Database, Cluster_Min_Distance, Cluster_Max_Distance);

    Length_Cluster = length(Pos_Array);
    for ii=1:Length_Cluster
        Cluster_Vector(ii,:) = Database(Pos_Array(ii),:);
        Cluster_Rank(ii,:) = Database_Rank(Pos_Array(ii),:);
    end
    Cluster_Vector = Cluster_Vector(1:Length_Cluster,:);
    Cluster_Rank = Cluster_Rank(1:Length_Cluster,:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 2: Locate Position - Using Gaussian Weight
    % Input : 
    % - Mean_Current_RSSI : Array of Mean RSSI for each APs 
    % - Cluster_Array:  X Y Z MAC1 Mean_RSSI1 MAC2 Mean_RSSI2 ...   
    % - Num_Mac (Number of APs in each vector)
    % - Value of K (1,2,3)
    % Output: 
    % X Y Z
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % KNN without Prediction
    EstPosition  = SRL_KNN_Mean_Rank(Core_Weight_Point, sigma, Current_RSSI, Rank_Current_RSSI, ...
                        Cluster_Vector,Cluster_Rank, Length_Cluster, Num_Mac,Num_Rank);

    Estimated_Position(CountPoint,:) =   EstPosition;            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update History
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if IdealFlag == 0 % Real History
        History = EstPosition; % Update recent point
    else
        History = IdealHistory_Array(CountPoint,:); % Update recent point
    end
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Distance Error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    True_Position = [TestPoint(CountPoint,1) TestPoint(CountPoint,2)];
    Distance_Error = (EstPosition(1)-True_Position(1))^2+(EstPosition(2)-True_Position(2))^2;
    Distance_Error_Meters(CountPoint) = sqrt(Distance_Error);
end
Average_Error = sum(Distance_Error_Meters)/length(Distance_Error_Meters)

% Variance Calculation 
Var_Dis_Error = 0;
for jj=1:length(Distance_Error_Meters)
    Var_Dis_Error = Var_Dis_Error + (Distance_Error_Meters(jj)-Average_Error)^2;
end
Var_Dis_Error = Var_Dis_Error/length(Distance_Error_Meters);
Std_Dis_Error = sqrt(Var_Dis_Error)

% CDF 
Theshold_Array = 0:0.5:max(Distance_Error_Meters)+1;
CDF_Array = zeros(1,length(Theshold_Array));

for ii= 1:length(Theshold_Array)
    Count_CDF = 0;
    for jj=1:length(Distance_Error_Meters)
        if Distance_Error_Meters(jj) <= Theshold_Array(ii)
            Count_CDF = Count_CDF + 1;
        end
    end
    CDF_Array(ii) = Count_CDF/length(Distance_Error_Meters);
end
figure,plot(Theshold_Array,CDF_Array);
title('KNN With Mean Fingerprint');
xlabel('Distance Error (meter)');
ylabel('CDF'); 

% Draw Ground Truth
figure,plot(Database(:,1),Database(:,2),'b*');
hold on;
plot(IdealHistory_Array(:,1),IdealHistory_Array(:,2),'r-');
hold on;
plot(Estimated_Position(:,1),Estimated_Position(:,2),'k-');
legend('RP','Ground Truth','Predicted Trajectory');
