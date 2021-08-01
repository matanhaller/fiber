close all;
clc;

epsilonR = reshape([1, 2],1 ,1, 1, []); 
eta = 1.8;
N = 20;
K = 100;
R = 1;
% x = linspace(-3, 3, 40); y = linspace(-3, 3, 40); z = 0;
r = linspace(1e-3, 1, 40);
theta = linspace(0, 2*pi, 100);
z = reshape(linspace(-5,5,4),1,1,[], 1);
[R, T] = meshgrid(r, theta);
X = R.*cos(T); Y = R.*sin(T);
% [X, Y] = meshgrid(x, y);


int = intensity(epsilonR, eta, X, Y, z, N, K, 1e-3);
energy = reshape(storedEnergy(int, R, r(2)-r(1), theta(2)-theta(1)),size(epsilonR, 4),[]);
for i = 1:size(epsilonR, 4)
    figure; plot(reshape(z,1,[]), reshape(energy(i,:),1,[]));
    hold on
    xlabel('z');
    ylabel('E');
    title('Stored Energy in The Cylinder as a Function of distance from center');
    xlim([-5, 5]);
end
% finish!!
% fplot(@(x) 0.5051 / (x^2 + 0.3887));

% legend('Stored Energy', 'Fitted Lorentzian');
% saveas(gcf, 'C:\Users\yansp\Dropbox\Matan and Yaniv\Static\Graphs\StoredEnergy.jpg');
%%
% for i=1:4
%     er = epsilonR(i);
%     int = intensity(er, eta, X, Y, z, N, K, 1e-3);
%     
%     subplot(2, 2, i); hold on;
%     surf(X, Y, int, 'EdgeColor', 'none');
%     view(2);
%     colorbar;
%     title(sprintf('$\\varepsilon=%d$', er), 'FontSize', 14, 'Interpreter', 'latex');
%     xlabel('$x$', 'FontSize', 14, 'Interpreter', 'latex');
%     ylabel('$y$', 'FontSize', 14, 'Interpreter', 'latex');
% %              saveas(gcf, 'C:\Users\yansp\Dropbox\Matan and Yaniv\Static\Graphs\Intensity.jpg');
% end


