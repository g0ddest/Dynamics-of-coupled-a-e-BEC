%% Исходные данные
% Константы
hbar    = 1; %1.05457172647e-34;   % Постоянная Планка (нормированная)
i_m = sqrt(-1);
kg = 1e34/9.109;
mr = 1e11/5.292;
s = 1e17/2.419;

% Параметры задачи
chi     = hbar*7.6e-7*mr^1.5/s;
Gamma1  = 1.629e-23*mr^3/s;
Gamma2  = 304.4/s;
beta1   = 2.108e-7/s;
beta2   = 3.344e-6/s;
alpha   = 134.06/s;
omega   = 2*pi*100/s;         	% Trap frequency
N       = 5e5;                  	% Initial number of atoms

m       = 100;		  				% atomic mass???
a       = 5.7e-9*mr;                	% s-wave scattering length (длина рассеивания s-волны)
aHO     = sqrt(hbar/(2*m*omega));	% see Ref 23

alpha1  = alpha/omega;
beta11  = beta1/omega;
beta21  = beta2/omega;

Gamma11 = Gamma1/(omega*aHO^2);
Gamma21 = Gamma2/omega;

lambda_a    = 4*pi*hbar.^2*a/m;	% сила взаимодействия атом-атом
lambda_m    = 4*pi*hbar.^2*a/m;	% молекула-молекула
lambda_am   = 4*pi*hbar.^2*a/m;	% атом-молекула

lambda_a1   = lambda_a/(aHO.^2*hbar*omega);
lambda_m1   = lambda_m/(aHO.^2*hbar*omega);
lambda_am1  = lambda_am/(aHO.^2*hbar*omega);

epsilon     = 1;    % o_O WTF?
epsilon1    = epsilon/(hbar*omega);

% Параметры метода
h       = 1e-3;                         % Шаг по пространству
dt      = 1e-4;                         % Шаг по времени

xmax    = 15;                           % Максимальная координата по пространству
T       = 1e-3;                         % Конечный момент времени

Jmax    = double(floor(xmax / h));      % Точек сетки по пространству
kmax    = floor(T / dt);                % Точек сетки по времени

x       = linspace(0, xmax, Jmax);      % Множество значений x в узлах сетки
tplot  = linspace(0, T, kmax);

% Начальные условия
phim    = zeros(kmax, Jmax);            % \phi_m (слои по времени — в строках) -- равен нулю.
phia    = zeros(kmax, Jmax);            % \phi_a (слои по времени — в строках) -- заполнить гауссианом (возможно лоренцаном на будущее).

%phim(1, :) = sin((1:Jmax) / Jmax * pi);
phia(1, :) = gaussmf(x, [xmax/12 0]);%sin((1:Jmax) / Jmax * pi);

%% Расчёт

%plot(abs(phia(1, :)), 'r');
%plot(abs(phim(1, :)), 'b');

for k = 1:kmax-1                        % Шаг по времени

    hold on
    plot(abs(phia(k, :)), 'r');
    plot(abs(phim(k, :)), 'b');

    PHIAa   = sparse(Jmax, Jmax);       % Коэффициенты для нахождения phia
    PHIAb   = zeros(Jmax, 1);           % Свободные члены для нахождения phia

    PHIMa   = sparse(Jmax, Jmax);       % Коэффициенты для нахождения phim
    PHIMb   = zeros(Jmax, 1);           % Свободные члены для нахождения phim

    for j = 2:Jmax-1                    % Шаг по пространству для phia
        
        % O_o WTF???
        PHIAa(1, 1)       = 1;                  % Левое граничное условие
        PHIAb(1)          = 0;                  %   (x=0)
        PHIAa(Jmax, Jmax) = 1;                  % Левое граничное условие
        PHIAb(Jmax)       = exp(-x(Jmax)^2/4);  %   (x=?)
        
        % Большая квадратная скобка для phia
        Za  = 0.5 * x(j) ^ 2 + lambda_a1 * (...   
                1 + 32 * a ^ 1.5 / (3 * sqrt(pi)) * abs(phia(k, j)) / x(j)...
            ) * abs(phia(k, j)) ^ 2 / (x(j) ^ 2) +...
            lambda_am1 * abs(phim(k, j)) ^ 2 / (x(j) ^ 2) - i_m * alpha1 -...
            beta11 - i_m * Gamma11 * abs(phia(k, j)) ^ 2 / (x(j) ^ 2);

        assert(~isnan(Za), sprintf('Oops, I did it again! (%d, %f)', j, x(j)));
        
        PHIAa(j, j-1)   = -1 / (4 * h ^ 2);
        PHIAa(j, j+1)   = -1 / (4 * h ^ 2);
        PHIAa(j,j)      = Za / 2 - i_m / dt + 1 / (2 * h ^ 2) + alpha1 * phim(k, j) / x(j);
        PHIAb(j)        = (-i_m / dt + 0.5 / (h ^ 2) + 0.5 * Za) * phia(k, j) +...
            1 / (4 * h ^ 2) * (phia(k, j+1) + phia(k, j-1)) ;% * conj(phia(k, j));
         
        % O_o WTF???
        PHIMa(1, 1)       = 1;                  % Левое граничное условие
        PHIMb(1)          = 0;                  %   (x=0)
        PHIMa(Jmax, Jmax) = 1;                  % Левое граничное условие
        PHIMb(Jmax)        = exp(-x(Jmax)^2/4);  %   (x=?)
         
        % Большая квадратная скобка для phim
        Zb  = x(j) ^ 2 + epsilon1 + lambda_m1 * abs(phim(k, j)) ^ 2 / x(j) ^ 2 +...
            lambda_am1 * abs(phia(k, j)) ^ 2 / x(j) ^ 2 - beta21 - i_m * Gamma21;
        
        PHIMa(j, j-1)   = -1 / (8 * h ^ 2);
        PHIMa(j, j+1)   = -1 / (8 * h ^ 2);
        PHIMa(j, j)     = -i_m / dt + 1 / (4 * h ^ 2) + Zb / 2;
        PHIMb(j)        = -i_m / dt + (phim(k, j+1) - 2 * phim(k, j) + phim(k, j-1)) / (8 * h ^ 2) -...
            Zb * phim(k, j) / 2 - alpha1 / 2 * phia(k, j) ^ 2 / x(j);
    end

    A = (PHIAa \ PHIAb);
    phia(k+1, :) = A.';   % Решаем трёхдиагональную систему
    phim(k+1, :) = (PHIMa \ PHIMb).';   % Решаем трёхдиагональную систему


end

figure(3); clf;
mesh( tplot, x, abs(phia'));  


colormap([0 0 0]);