function facedemo(n)

% FACETEST: Example polygonal geometries for MESHFACES.
%
%  facedemo(n)
%
% Darren Engwirda - 2007

switch n
   case 1

      node = [0.0, 0.0; 1.0,0.0; 1.0,1.0; 0.0,1.0; 1.01,0.0; 1.01,1.0; 3.0,0.0; 3.0,1.0];
      edge = [1,2; 2,3; 3,4; 4,1; 2,5; 5,6; 6,3; 5,7; 7,8; 8,6];
      face{1} = [1,2,3,4];
      face{2} = [5,6,7,2];
      face{3} = [8,9,10,6];

      meshfaces(node,edge,face);

   case 2

      % Geometry
      dtheta = pi/36;
      theta = (-pi:dtheta:(pi-dtheta))';
      node1 = [cos(theta), sin(theta)];
      node2 = [-2.0,-2.0; 2.0,-2.0; 2.0,2.0; -2.0, 2.0];
      edge1 = [(1:size(node1,1))',[(2:size(node1,1))'; 1]];
      edge2 = [1,2; 2,3; 3,4; 4,1];

      edge = [edge1; edge2+size(node1,1)];
      node = [node1; node2];

      face{1} = 1:size(edge1,1);
      face{2} = 1:size(edge,1);

      meshfaces(node,edge,face);

   otherwise
      error('Invalid demo. N must be between 1-2');
end

end      % facetest()
