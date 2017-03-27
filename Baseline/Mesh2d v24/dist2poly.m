function L = dist2poly(p,edgexy,lim)

% Find the minimum distance from the points in P to the polygon defined by
% the edges in EDGEXY. LIM is an optional argument that defines an upper
% bound on the distance for each point.

% Uses (something like?) a double sweep-line approach to reduce the number
% of edges that are required to be tested in order to determine the closest
% edge for each point. On average only size(EDGEXY)/4 comparisons need to
% be made for each point.

if nargin<3
   lim = [];
end
np = size(p,1);
ne = size(edgexy,1);
if isempty(lim)
   lim = inf*ones(np,1);
end

% Choose the direction with the biggest range as the "y-coordinate" for the
% test. This should ensure that the sorting is done along the best
% direction for long and skinny problems wrt either the x or y axes.
dxy = max(p)-min(p);
if dxy(1)>dxy(2)
    % Flip co-ords if x range is bigger
    p       = p(:,[2,1]);
    edgexy  = edgexy(:,[2,1,4,3]);
end

% Ensure edgexy(:,[1,2]) contains the lower y value
swap           = edgexy(:,4)<edgexy(:,2);
edgexy(swap,:) = edgexy(swap,[3,4,1,2]);

% Sort edges
[i,i]          = sort(edgexy(:,2));                                        % Sort edges by lower y value
edgexy_lower   = edgexy(i,:);
[i,i]          = sort(edgexy(:,4));                                        % Sort edges by upper y value
edgexy_upper   = edgexy(i,:);

% Mean edge y value
ymean = 0.5*( sum(sum(edgexy(:,[2,4]))) )/ne;

% Alloc output
L = zeros(np,1);

% Loop through points
tol = 1000.0*eps*max(dxy);
for k = 1:np

   x = p(k,1);
   y = p(k,2);
   d = lim(k);

   if y<ymean

      % Loop through edges bottom up
      for j = 1:ne
         y2 = edgexy_lower(j,4);
         if y2>=(y-d)
            y1 = edgexy_lower(j,2);
            if y1<=(y+d)

               x1 = edgexy_lower(j,1);
               x2 = edgexy_lower(j,3);

               if x1<x2
                  xmin = x1;
                  xmax = x2;
               else
                  xmin = x2;
                  xmax = x1;
               end

               if xmin<=(x+d) && xmax>=(x-d)
                  % Calculate the distance along the normal projection from [x,y] to the jth edge
                  x2mx1 = x2-x1;
                  y2my1 = y2-y1;

                  r = ((x-x1)*x2mx1+(y-y1)*y2my1)/(x2mx1^2+y2my1^2);
                  if r>1.0                                                 % Limit to wall endpoints
                     r = 1.0;
                  elseif r<0.0
                     r = 0.0;
                  end

                  dj = (x1+r*x2mx1-x)^2+(y1+r*y2my1-y)^2;
                  if (dj<d^2) && (dj>tol)
                     d = sqrt(dj);
                  end

               end

            else
               break
            end
         end
      end

   else

      % Loop through edges top down
      for j = ne:-1:1
         y1 = edgexy_upper(j,2);
         if y1<=(y+d)
            y2 = edgexy_upper(j,4);
            if y2>=(y-d)

               x1 = edgexy_upper(j,1);
               x2 = edgexy_upper(j,3);

               if x1<x2
                  xmin = x1;
                  xmax = x2;
               else
                  xmin = x2;
                  xmax = x1;
               end

               if xmin<=(x+d) && xmax>=(x-d)
                  % Calculate the distance along the normal projection from [x,y] to the jth edge
                  x2mx1 = x2-x1;
                  y2my1 = y2-y1;  

                  r = ((x-x1)*x2mx1+(y-y1)*y2my1)/(x2mx1^2+y2my1^2);
                  if r>1.0                                                 % Limit to wall endpoints
                     r = 1.0;
                  elseif r<0.0
                     r = 0.0;
                  end

                  dj = (x1+r*x2mx1-x)^2+(y1+r*y2my1-y)^2;
                  if (dj<d^2) && (dj>tol)
                     d = sqrt(dj);
                  end

               end

            else
               break
            end
         end
      end

   end

   L(k) = d;

end

end      % dist2poly()
