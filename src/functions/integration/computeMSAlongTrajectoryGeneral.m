function [t,uncertAtTime] = computeMSAlongTrajectoryGeneral(derivative, derivativeEOV, initialPoint, timeInterval, timeStep, toleranceFD, toleranceInt, Deltas, useReduced, r)
t = timeInterval(1):timeStep:timeInterval(2);
n = length(t);
uncertAtTime = zeros(size(t));
for i=1:n
    %% calculate over growing intervals
    iv = [timeInterval(1), timeInterval(1) + i*timeStep];
    % 10 steps per itnerval, resolution becomes [1,1] since we only follow
    % 1 point
    disp(i)
    temp = computeModelSensitivity(derivative, derivativeEOV, initialPoint, iv, toleranceFD, toleranceInt, Deltas, false, useReduced, r);
    %disp([iv, temp]);

    uncertAtTime(i) = temp;
end
uncertAtTime = transpose(uncertAtTime);

end