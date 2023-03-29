function targetPosition_after = classPos(targetPosition_before)

targetPosition_after(targetPosition_before <= 9) = 0; % LVF
targetPosition_after(targetPosition_before > 9) = 1; % RVF

end