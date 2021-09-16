%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021


function CursorPositionOutput = PlayNeuralDataThroughKalmanEML(KM_features,StartTimeIndex_IC,EndTimeIndex_IC,A_KM,H_KM,K_KM)
%%%%%Run the neural data through the Kalman Decoder Model

%z-score each of the features
numfeatures = size(KM_features,2);
local_KM_features = double(KM_features(StartTimeIndex_IC:EndTimeIndex_IC,:));
z_KM_features = 0*local_KM_features;

for i = 1:numfeatures
    z_KM_features(:,i) = local_KM_features(:,i) - mean(local_KM_features(:,i));
    z_KM_features(:,i) = z_KM_features(:,i)./std(z_KM_features(:,i));
end

%flip dimensions to fit Kalman EML script
z_KM_features = z_KM_features';

CursorPositionOutput = zeros(3,(EndTimeIndex_IC - StartTimeIndex_IC + 1));

%Run kalman data over time t
for t = 2:(EndTimeIndex_IC - StartTimeIndex_IC + 1)    
    z = z_KM_features(:,t);    
    x_est = kalmanDecode_EML(z,CursorPositionOutput(:,t-1),A_KM,H_KM,K_KM);        
    CursorPositionOutput(:,t) = x_est;    
end

end
