%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021


function [KM_One_Template_adj,KM_Two_Template_adj] = GenerateIdealizedTemplates(thisTemplateCompressionFactor,KM_One_Template,KM_Two_Template)

%% Generate some code to create splines fitted to templates, and use these splines to squeeze and dilate the templates for matching process

SampleRate = 50;

x = (1:length(KM_One_Template)) - 1;
x = x./SampleRate;

%bigger x is just for plotting
bigger_x = (1:10000) - 1;
bigger_x = bigger_x./SampleRate;

%fit the template to a cubic spline function
cs1 = csapi(x,KM_One_Template);

%In order to shrink/expand the projetion from space x with samplerate y,
%need to first evaluate the function over the normal range ([0 end(x)])

%Need to chose a sample rate for the new x axis such that when divided by
%the new end time, the final sample rate ends up matching the original

%for example, to speed up the template by 2x, would want the curve to end
%at 2 seconds instead of four seconds (i.e. end(x)/2). Thus, would want the
%sample rate to be twice big as before:
    % new_x = x(1):SampleRate*2:x(end)
    % new_Template = fnval(cs,new_x)
    % plot(x(1:length(new_x)),new_Template)

    %When projected into the orignal X space, the new template will cover
    %the same trajectory but do so over fewer or more points, effectively
    %dilating/compressing the template
    
new_x = x(1):(1/SampleRate)*thisTemplateCompressionFactor:x(end);
KM_One_Template_adj = fnval(cs1,new_x);
    

cs2 = csapi(x,KM_Two_Template);

new_x = x(1):(1/SampleRate)*thisTemplateCompressionFactor:x(end);
KM_Two_Template_adj = fnval(cs2,new_x);
    

% 
% figure;
% hold on
% 
% fnplt(cs1);
% plot(x,KM_One_Template,'bo')
% plot(bigger_x(1:length(new_x)),KM_One_Template_adj)
% 
% fnplt(cs2);
% plot(x,KM_Two_Template,'ro')
% plot(bigger_x(1:length(new_x)),KM_Two_Template_adj)
% 
% hold off

end

