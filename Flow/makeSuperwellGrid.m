function makeSuperwellGrid(G, newData)
%   MAKESUPERWELLS Visual representation of partition grid that overlays
%   the wilmington wells in order to create superwells
    load('/Users/annikahuprikar/Desktop/Spring2021/HSURV PRISE Summer 2021/Flow/fromMatlab/mrstResults/data/state0.mat')
    data1=data;
    superwellGrid = figure('Name', '2D Superwell Grid')
    % plotCellData(G, data1.pressure / 1e6)
    % view(2)
    xlim([3.76e05, 4.01e05]); ylim([3.727e06, 3.752e06]); % 4 wells have y = 0, x = 1.66e+05
    hold on;
    scatter(newData(:, 1), newData(:, 2))
    title("Well Locations in XY Cartesian Plane");
    xlabel("X Position"); ylabel("Y Position");
    % xlim([3.8e05, 3.97e05]); ylim([3.731e06, 3.742e06]); % 4 wells have y = 0, x = 1.66e+05
    hold off;

    
    dx = 0.01e05 / 5; dy = 0.001e06 / 5;
    cur_x = 3.76e05; end_x = 4.01e05;
    cur_y = 3.727e06; end_y = 3.752e06; 
    
    %{
    dx = 0.01e05 / 5; dy = 0.001e06 / 5;
    cur_x = 3.8e05; end_x = 3.97e05;
    cur_y = 3.731e06; end_y = 3.742e06; 
    %}
    x_interval = cur_x : end_x; y_interval = cur_y : end_y;
    
    %{
    dx = 0.01e05;
    dy = 0.001e06;
    cur_x = 3.8e05; end_x = 3.95e05;
    cur_y = 3.732e06; end_y = 3.742e06;
    x_interval = cur_x : end_x; y_interval = cur_y : end_y;
    %}
    
    figure(superwellGrid);
    hold on;
    while cur_x <= end_x
        plot(cur_x*ones(size(y_interval)), y_interval, '-k');
        cur_x = cur_x + dx;
    end

    while cur_y <= end_y
        plot(x_interval, cur_y*ones(size(x_interval)), '-r');
        cur_y = cur_y + dy;
    end
    hold off;
end