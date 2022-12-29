% Function counts number of nodes in the wrong community. 
% It is used to calculate the accuracy of a partition.
% Input: E vector expected communities
%        P vector predicted communities
% Output: n number of nodes in the wrong community 

function [n] = wrong(E,P)

%[T] = confusion_matrix(E,P); %confusion matrix, coded by Sara Venturini
[T] = confusionmat(E,P); %confusion matrix, Matlab function

n=0; %counts nodes in wrong community
while ~isempty(T)
    M = max(max(T)); %maximum value of T
    [x,y] = find(T==M); %indices of M in T
    %use x(1) and y(1) because maybe more entries correspond to maximum value 
    n = n + sum(T(x(1),:)) - T(x(1),y(1));
    n = n + sum(T(:,y(1))) - T(x(1),y(1));
    T(x(1),:) = [];
    T(:,y(1)) = [];
end

end

