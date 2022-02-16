%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get The Set of Possible Locations from the Original Database 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Pos_Array = NearDistanceNeighbour(OldPosition, Cluster_Vector, Length_Cluster, Threshold_Min_Dis, Threshold_Max_Dis)

    Min_Dis = Threshold_Min_Dis;
    Max_Dis = Threshold_Max_Dis;
    Threshold_Min_DisSquare = Threshold_Min_Dis^2;
    Threshold_Max_DisSquare = Threshold_Max_Dis^2;
    
    Pos_X = OldPosition(1);
    Pos_Y = OldPosition(2);
    
    Distance_Array = 100+zeros(1,Length_Cluster);
    Pos_Array_Temp = zeros(1,Length_Cluster);
    CountPos = 1;
    for ii = 1:Length_Cluster % Calculate Distance
        Pos_X_Temp = Cluster_Vector(ii,1); % X
        Pos_Y_Temp = Cluster_Vector(ii,2); % Y
        
        %%%%% Get the subset of possible locations from original database  
        if (abs(Pos_X_Temp-Pos_X)<=Max_Dis) ... % Near Neighbour
                || (abs(Pos_Y_Temp-Pos_Y)<=Max_Dis)
            Distance_Array(CountPos) = (Pos_X_Temp-Pos_X)^2 + (Pos_Y_Temp-Pos_Y)^2;
            
            if((Distance_Array(CountPos) >= Threshold_Min_DisSquare)&&(Distance_Array(CountPos) <= Threshold_Max_DisSquare))
                Pos_Array_Temp(CountPos) = ii;
                CountPos = CountPos + 1;
            end
        end
    end
    
    Pos_Array = zeros(1, CountPos-1);
    for jj = 1:CountPos-1
        Pos_Array(jj) = Pos_Array_Temp(jj);
    end
end
