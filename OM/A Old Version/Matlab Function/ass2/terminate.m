function [lookfor stop direction] = terminate(t,s)
% This function specifies the event at which  terminates.
% -----------------------------------------------------------------------
lookfor = norm(s(1:3)) - 100 - 6371; % ¼ 0 when altitude ¼ 100 km
stop = 1; % 1 means terminate at lookfor ¼ 0; Otherwise 0
direction = -1; % -1 means zero crossing is from above
end %terminate