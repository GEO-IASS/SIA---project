function [] = loopoverant( varargin )
%UNTITLED this function loops over the ant script and allows user to
%iteratively change parameters such as temperature, grid size, time
%incrementation, etc

%   Input: 
% specifier followed by value spread: 
% e.g. 
% 'time blocks',[10:100] % in years
% 'grid resolution',[10:100] % in km

changingVal = varargin{2};
N = length(changingVal);

for n = 1:N;
    ant('trend','hotter', 'variable temperature','on',...
        varargin{1},changingVal(n));
end

end

