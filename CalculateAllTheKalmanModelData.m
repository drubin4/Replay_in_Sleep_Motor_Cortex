%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021


%% Generate Kalman Data output as a function of input data and time constant

function Kalman_Data = CalculateAllTheKalmanModelData(A_adj);

S = load('KalmanModelData.mat');
%A_KM = KalmanModel.A;
H_KM = S.KalmanModel.H;
K_KM = S.KalmanModel.K;

A_KM = A_adj;

ncTXInds = S.KalmanModel.ncTXInds;
spikePowerInds = S.KalmanModel.spikePowerInds;

S = load('ncTX_Both_Array.mat');
ncTX_KM = S.ncTX_Both_Arrays(:,ncTXInds);
clear S

S = load('SpikePower_Both_Arrays.mat');
spikePower_KM = S.SpikePower_Both_Arrays(:,spikePowerInds);
TimeAxis = S.TimeAxis;
clear S

KM_features = (cat(2,ncTX_KM,spikePower_KM));

 timeStart = datetime(2020,11,23,23,0,0);
 timeEnd = datetime(2020,11,25,9,0,0);
 
 duration_in_minutes = minutes(timeEnd - timeStart);
 analysis_interval = 1;                                  %unit is minutes
 number_of_analysis_segments = ceil(duration_in_minutes/analysis_interval);
 
DateNum_IC_TimeAxis = datenum(TimeAxis);

Kalman_Data = zeros(3,length(DateNum_IC_TimeAxis));

for i = 1:number_of_analysis_segments
     
     segmentTimeStart = timeStart + minutes(analysis_interval)*(i-1);
     segmentTimeEnd = segmentTimeStart+minutes(2);
    
     %%% First find if there is IC data present during the time interval
     
    StartTimeIndex_IC = find(TimeAxis > segmentTimeStart);
    StartTimeIndex_IC = StartTimeIndex_IC(1);

    EndTimeIndex_IC = find(TimeAxis > segmentTimeEnd);
    EndTimeIndex_IC = EndTimeIndex_IC(1);
    
    if EndTimeIndex_IC > StartTimeIndex_IC      %This ensures there is some data that exists to analyze
        
        Xlim1 = segmentTimeStart;
        Xlim2 = segmentTimeEnd;
        
         StartTimeIndex_IC = find(TimeAxis > Xlim1);
         StartTimeIndex_IC = StartTimeIndex_IC(1);

        EndTimeIndex_IC = find(TimeAxis > Xlim2);
        EndTimeIndex_IC = EndTimeIndex_IC(1);

            %%%%%Run the data through the Kalman Decoder Model
           CursorPositionOutput = PlayNeuralDataThroughKalmanEML(KM_features,StartTimeIndex_IC,EndTimeIndex_IC,A_KM,H_KM,K_KM);    

            Kalman_Data(:,StartTimeIndex_IC:EndTimeIndex_IC) = CursorPositionOutput;
  
    end  
   i/number_of_analysis_segments  
end

end