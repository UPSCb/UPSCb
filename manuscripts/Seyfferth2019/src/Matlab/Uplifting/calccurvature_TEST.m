% calccurvature_TEST:
figure('Position', [481 255 1305 397])
%% Flat Line
X =  1:10;
Y = (1:10)*0;
hook=10;
curvature = calccurvature(X,Y,hook);
subplot(1,6,1), plot(X,Y)
title({'Flat Line', ['Curvature: ' num2str(curvature)]})

%% Line with uniform slope
X = 1:10;
Y = 1:10;
hook = 10;
curvature = calccurvature(X,Y,hook);
subplot(1,6,2), plot(X,Y), title(curvature)
title({'Line with uniform slope', ['Curvature: ' num2str(curvature)]})


%% Line with changing slope (ie curved)
X = 1:10;
Y = (X.^2)/10;
hook = 10;
curvature = calccurvature(X,Y,hook);
subplot(1,6,3),plot(X,Y), title(curvature)
title({'Line with changing slope ', ['Curvature: ' num2str(curvature)]})

%% Line with changing slope (ie curved)
X = 1:10;
Y = (X.^4)/1000;
hook = 10;
curvature = calccurvature(X,Y,hook);
subplot(1,6,4),plot(X,Y), title(curvature)
title({'Line with changing slope ', ['Curvature: ' num2str(curvature)]})


%% Negative slope line
X = 1:10;
Y = -X.^4/1000;
hook = 10;
curvature = calccurvature(X,Y,hook);
subplot(1,6,5),plot(X,Y), title(curvature)
title({'Negative slope line ', ['Curvature: ' num2str(curvature)]})

%% Half circle, radius 2
th = linspace( pi/2, -pi/2, 10);
R = 2;
X = R*cos(th);
Y = R*sin(th)+R;
hook = 10;
curvature = calccurvature(X,Y,hook);
subplot(1,6,6), plot(X,Y), title(curvature)
title({'Half circle, radius 2 ', ['Curvature: ' num2str(curvature)]})

pltstmp('curvature TEST')