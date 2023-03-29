function [corrResp, idxCorrResp] = corrResponse(targetSide,subjectRes)

    % by comparing the response data and the tilt of target, will generate an
    % array of "correct response" in x3 (trial dimension).
    % 0 = wrong; 1 = correct
    corrResp = subjectRes == targetSide; % the correct trials

    % using the information above, the correct response indices are obtained:
    idxCorrResp = find(corrResp == 1); % according the correct response
    
% thisCR = subjectRes(idxCorrResp);
%     



end

