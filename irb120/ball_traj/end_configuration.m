function [q_itc] = end_configuration(ball_intersection, robot_model)
%function [q_itc] = end_configuration(ball_intersection, robot_model)
%   end_configuration solves the inverse kinematic problem and returns a
%   configuration for the robot that catches the ball. z axis of the end
%   effector is tangential to the ball's trajectory. x axis points
%   downwards and it is perpendicular to z. y is chosen such that 
%   y = cross(z, x)
%   However, this function might be slow because
%   1). a robot model has to be copied as a parameter
%   2). the inverse kinematics might not be easy to solve

z = -ball_intersection(4:6);
z = z/norm(z, 2);

if z(3) == 0
    x = [0; 0; -1];
else
    x = [1; 0; -z(1)/z(3)];
    x = x/norm(x, 2);
end

y = cross(z, x);

T = [[x;0], [y;0], [z;0], [ball_intersection(1:3);1]];

q_itc = robot_model.ikine(T);

end

% test
% q_itc = end_configuration(test_traj(:,idx), irb120)
