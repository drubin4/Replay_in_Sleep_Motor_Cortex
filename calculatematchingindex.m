%% Code Supporting the publication "Learned motor patterns are replayed in human motor cortex during sleep" by Rubin et al.
%% Copyright Daniel B. Rubin 2021


function Imatch = calculatematchingindex(template,frame)
% Calculate matching index per method of Ji and Wilson (2007)


%M is the total number of channels (96 here) in a given replay event (e.g.
%all of them).  There are a total of M(M-1)/2 possible cell pairs
% if "m" is the number of pairs that have the same order between a template
% and a frame, and "n" is the number of pairs with opposite order, then

% I = (m - n)/(m + n) 

% Where I is a matching index bounded by [-1 to 1], where 1 is a perfect
% match and -1 is an inverted match

m = 0;
n = 0;

for i = 1:length(template)
   for j = 1:length(template)
       if i ~= j
          
           a = find(frame == template(i));
           b = find(frame == template(j));
           
           if (i < j) && (a < b)
               m = m+1;
           end
           if (i > j) && (a < b)
               n = n+1;
           end
           if (i < j) && (a > b)
               n = n+1;
           end
           if (i > j) && (a > b)
               m = m+1;
           end
       end
       
   end
end
    
Imatch = (m - n)/(m + n);


end

