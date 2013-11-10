%% Исходные данные
% Константы
hbar    = 1; %1.05457172647e-34;   % Постоянная Планка (нормированная)

% Параметры задачи
chi     = hbar * 7.6 * 1.44/5.291772 * 2.418884 * 1e-9 ;
Gamma1  = 1.629 * 2.418884/5.291772/5.291772/5.291772 * 1e-7;
Gamma2  = 304.4 * 2.418884 * 1e-17;
beta1   = 2.108 * 2.418884 * 1e-10;
beta2   = 3.344 * 2.418884 * 1e-17;
alpha   = 134.06 * 2.418884 * 1e-11;
omega   = 2 * pi * 100  * 2.418884 * 1e-17;         	% Trap frequency
N       = 5e5;                  	% Initial number of atoms

m       = 85;		  				% atomic mass???
a       = 17e-9 /5.291772 *1e11 ;                	% s-wave scattering length (длина рассеивания s-волны)
aHO     = sqrt(hbar/(2*m*omega));	% see Ref 23

alpha1  = alpha / omega;
beta11  = beta1 / omega;
beta21  = beta2 / omega;

Gamma11 = Gamma1 / (omega * aHO ^ 2);
Gamma21 = Gamma2 / omega;

lambda_a    = 4 * pi * hbar .^ 2 * a / m;	% сила взаимодействия атом-атом
lambda_m    = 4 * pi * hbar .^ 2 * a / m;	% молекула-молекула
lambda_am   = 4 * pi * hbar .^ 2 * a / m;	% атом-молекула

lambda_a1   = lambda_a / (aHO .^ 2 * hbar * omega);
lambda_m1   = lambda_m / (aHO .^ 2 * hbar * omega);
lambda_am1  = lambda_am / (aHO .^ 2 * hbar * omega);

epsilon     = 1;    % o_O WTF?
epsilon1    = epsilon / (hbar * omega);

% Параметры метода
h       = 1e-3;                         % Шаг по пространству
dt      = 1e-4;                         % Шаг по времени

xmax    = 15;                           % Максимальная координата по пространству
T       = 1e-3;                         % Конечный момент времени

Jmax    = double(floor(xmax / h));      % Точек сетки по пространству
kmax    = floor(T / dt);                % Точек сетки по времени

x       = linspace(0, xmax, Jmax);      % Множество значений x в узлах сетки

% Начальные условия
phim    = zeros(kmax, Jmax);            % \phi_m (слои по времени — в строках) -- равен нулю.
phia    = zeros(kmax, Jmax);            % \phi_a (слои по времени — в строках) -- заполнить гауссианом (возможно лоренцаном на будущее).

%phim(1, :) = sin((1:Jmax) / Jmax * pi);
phia(1, :) = gaussmf(x, [xmax/12 0]);%sin((1:Jmax) / Jmax * pi);

%% Расчёт
for k = 1:kmax-1                        % Шаг по времени
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
            lambda_am1 * abs(phim(k, j)) ^ 2 / (x(j) ^ 2) - 1i * alpha1 -...
            beta11 - 1i * Gamma11 * abs(phia(k, j)) ^ 2 / (x(j) ^ 2);

        assert(~isnan(Za), sprintf('Oops, I did it again! (%d, %f)', j, x(j)));
        
        PHIAa(j, j-1)   = -1 / (4 * h ^ 2);
        PHIAa(j, j+1)   = -1 / (4 * h ^ 2);
        PHIAa(j,j)      = Za / 2 - 1i / dt + 1 / (2 * h ^ 2);
        PHIAb(j)        = (-1i / dt + 0.5 / (h ^ 2) + 0.5 * Za) * phia(k, j) +...
            1 / (4 * h ^ 2) * (phia(k, j+1) + phia(k, j-1)) +...
             alpha1 * phim(k, j) * phia(k, j)' / x(j);
         
        % O_o WTF???
        PHIMa(1, 1)       = 1;                  % Левое граничное условие
        PHIMb(1)          = 0;                  %   (x=0)
        PHIMa(Jmax, Jmax) = 1;                  % Левое граничное условие
        PHIMb(Jmax)        = exp(-x(Jmax)^2/4);  %   (x=?)
         
        % Большая квадратная скобка для phim
        Zb  = x(j) ^ 2 + epsilon1 + lambda_m1 * abs(phim(k, j)) ^ 2 / x(j) ^ 2 +...
            lambda_am1 * abs(phia(k, j)) ^ 2 / x(j) ^ 2 - beta21 - 1i * Gamma21;
        
        PHIMa(j, j-1)   = -1 / (8 * h ^ 2);
        PHIMa(j, j+1)   = -1 / (8 * h ^ 2);
        PHIMa(j, j)     = -1i / dt + 1 / (4 * h ^ 2) + Zb / 2;
        PHIMb(j)        = -1i / dt + (phim(k, j+1) - 2 * phim(k, j) + phim(k, j-1)) / (8 * h ^ 2) -...
            Zb * phim(k, j) / 2 - alpha1 / 2 * phia(k, j) ^ 2 / x(j);
    end
    phia(k+1, :) = (PHIAa \ PHIAb).';   % Решаем трёхдиагональную систему
    phim(k+1, :) = (PHIMa \ PHIMb).';   % Решаем трёхдиагональную систему
	
	hold on
    plot(abs(phia(k+1, :)), 'r');
    plot(abs(phim(k+1, :)), 'b');
    title('k=2');
end

figure(3); clf;
mesh(x, tplot,  phia);  % Plot P(x,t) vs. x and t