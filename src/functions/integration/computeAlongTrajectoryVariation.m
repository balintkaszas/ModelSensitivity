function [t,uncertAtTime] = computeAlongTrajectoryVariation(derivative, initialPoint, timeInterval, stepSize )
t = timeInterval(1):stepSize:timeInterval(2);
n = length(t);
uncertAtTime = zeros(size(t));

for i=1:n
    %% calculate over growing intervals
    iv = [timeInterval(1), timeInterval(1) + i*stepSize];
    % 10 steps per itnerval, resolution becomes [1,1] since we only follow
    % 1 point
    temp = uncertEstimateVariational(derivative, initialPoint,[1,1], iv, stepSize);
    %disp([iv, temp]);

    uncertAtTime(i) = temp;
end
uncertAtTime = transpose(uncertAtTime);

end