function [h,hnorm] = angularMomentum(rr,vv)
% Function returns angular momemntum at each data point rr and vv

[row,col] = size(rr);

if row == 3
    % initialize
    h = zeros(3,col);
    hnorm = zeros(col,1);
    for i = 1:col
        % compute angular momentum
        rtmp = rr(:,i);
        vtmp = vv(:,i);
        htmp = cross(rtmp, vtmp);
        % store
        h(:,i) = htmp;
        hnorm(i,1) = norm(htmp); 
    end
elseif col == 3
    % initialize
    h = zeros(row,3);
    hnorm = zeros(row,1);
    for i = 1:row
        % compute angular momentum
        rtmp = rr(i,:)';
        vtmp = vv(i,:)';
        htmp = cross(rtmp, vtmp);
        htmp = htmp';
        % store
        h(i,:) = htmp;
        hnorm(i,1) = norm(htmp);
    end
end




end