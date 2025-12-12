function colorstyle = get_red_blue_color_map()
n_colors = 256;
light_blue = [0.2 0.3 1];
white = [1 1 1];
blue_to_white = interp1([0 1], [light_blue; white], linspace(0, 1, n_colors/2));
light_red = [1 0.2 0.2];
white_to_red = interp1([0 1], [white; light_red], linspace(0, 1, n_colors/2));
colorstyle = [blue_to_white; white_to_red];
