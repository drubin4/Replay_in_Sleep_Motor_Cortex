%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021


function [DistanceVector, CrossCorrVector] = CalculateDistanceFromTemplateOverTime(KalmanData_adj,KM_One_Template_adj,KM_Two_Template_adj)
% At each time step, how far is the current Kalman Model output from the
% idealized template, in terms of normalized euclidean distance:

VectorLength = size(KalmanData_adj,2);
TemplateLength = length(KM_One_Template_adj);

DistanceVector1 = zeros(VectorLength,1);
DistanceVector2 = zeros(VectorLength,1);

CrossCorr_TD1 = zeros(VectorLength,1);
CrossCorr_TD2 = zeros(VectorLength,1);

tic
% 
% normed_KM_One_Template_adj = KM_One_Template_adj/norm(KM_One_Template_adj);
% normed_KM_Two_Template_adj = KM_Two_Template_adj/norm(KM_Two_Template_adj);

parfor j=1:(VectorLength - TemplateLength)
    ThisSegmentOfData = KalmanData_adj(1:2,j:(j+TemplateLength-1));
    
    
    %Find the cosine between the activity and the template vectors
    
    DistanceVector1(j) = dot(KM_One_Template_adj,ThisSegmentOfData(1,:))./(norm(KM_One_Template_adj)*norm(ThisSegmentOfData(1,:)))
    
    DistanceVector2(j) = dot(KM_Two_Template_adj,ThisSegmentOfData(2,:))./(norm(KM_Two_Template_adj)*norm(ThisSegmentOfData(2,:)))
    
    
    %For comparison also find the correlation coefficent as per the primary
    %analysis:
    
     [r p] = corrcoef(KM_One_Template_adj,ThisSegmentOfData(1,:));
    CrossCorr_TD1(j) = r(1,2);
        
     [r p] = corrcoef(KM_Two_Template_adj,ThisSegmentOfData(2,:));
    CrossCorr_TD2(j) = r(1,2);
    

end
toc

DistanceVector = [DistanceVector1 DistanceVector2];
CrossCorrVector = [CrossCorr_TD1 CrossCorr_TD2];

end

