function [err_, bound_, Lintersect_] = boundErrPts(fits, data)
% 
% fits are either:
%   [0 y0 20 y1] for horizontal  bounds
%   [x0 0 x1 0.75] for vertical bounds

% data is matrix, rows are trials, columns are:
%   bound lower limit
%   bound upper limit
%   bound t1

% z-scored
% TMIN = -1;
% TMAX = 15;
% BMIN = -1.3;
% BMAX = 3.3;

% raw
TMIN = 0; 
TMAX = 20;
BMIN = 0;
BMAX = 0.8;

if fits(1) < 1
    bx0 = TMIN+(1-fits(1))*(TMAX-TMIN);
    by0 = BMIN;
else % 1-2
    bx0 = TMIN;
    by0 = BMIN + (fits(1)-1)*(BMAX-BMIN);
end

if fits(2) < 1
    bx1 = TMIN+fits(2)*(TMAX-TMIN);
    by1 = BMAX;
else % 1-2
    bx1 = TMAX;
    by1 = BMIN + (2-fits(2))*(BMAX-BMIN);
end

% x1, y1 -> bx0, by0
% x2, y2 -> bx1, by1
% x3, y3 -> data(:,3), data(:,1)
% x4, y4 -> data(:,3)+1, data(:,2)

% Line segments intersect parameters
% us = ((x1-x3)*(y1-y2) - (y1-y3)*(x1-x2)) / ((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
% ts = ((x1-x3)*(y3-y4) - (y1-y3)*(x3-x4)) / ((x1-x2)*(y3-y4)-(y1-y2)*(x3-x4));
us = ((bx0-data(:,3)).*(by0-by1) - (by0-data(:,1)).*(bx0-bx1)) ./ ((bx0-bx1).*(data(:,1)-data(:,2))-(by0-by1).*(data(:,3)-data(:,3)-1));
ts = ((bx0-data(:,3)).*(data(:,1)-data(:,2)) - (by0-data(:,1)).*(data(:,3)-data(:,3)-1)) ./ ((bx0-bx1).*(data(:,1)-data(:,2))-(by0-by1).*(data(:,3)-data(:,3)-1));


% Line segments intersect parameters

% Check if intersection exists
Lintersect = (us >= 0 & us <= 1.0) & (ts >= 0 & ts <= 1.0);

% err_ = sum(~Lintersect);

% disp(err_)
% cla reset; hold on;
% plot([bx0 bx1], [by0 by1], 'r-')
% for xx = 1:size(data,1)
%     h= plot([data(xx,3) data(xx,3)+1], data(xx,1:2), 'k-');
%     if Lintersect(xx)
%         set(h, 'Color', 'g')
%     end
% end
% axis([0 20 0 0.75])
% r = input('next')

% Compute shortest distance between the two segments
% point1s = [bx0 by0];
% point1e = [bx1 by1];
% point2s = [data(:,3) data(:,1)]
% point2e = [data(:,3)+1 data(:,2)]

if any(~Lintersect)

    v1 = repmat([bx0 by0 0], sum(~Lintersect), 1);
    v2 = repmat([bx1 by1 0], sum(~Lintersect), 1);
    a = v1 - v2;
    b1 = [data(~Lintersect,3) data(~Lintersect,1) zeros(sum(~Lintersect),1)] - v2;
    b2 = [data(~Lintersect,3)+1 data(~Lintersect,2) zeros(sum(~Lintersect),1)] - v2;

    denom = sqrt(sum(a.^2,2));
    d1s = sqrt(sum(cross(a,b1,2).^2,2))./denom;
    d2s = sqrt(sum(cross(a,b2,2).^2,2))./denom;
    Ld2 = d2s<d1s;
    % err_ = sum(d1s(~Ld2).^2) + sum(d2s(Ld2).^2);
    err_ = sum(abs(d1s(~Ld2))) + sum(abs(d2s(Ld2)));

    % plot(b1(~Ld2,1)+bx1, b1(~Ld2,2)+by1, 'r.')
    % plot(b2( Ld2,1)+bx1, b2( Ld2,2)+by1, 'r.')
else
    err_ = 0;
end

if nargout > 1
    bound_ = [bx0 by0 bx1 by1];
    if nargout > 2
        Lintersect_ = Lintersect;
    end
end
