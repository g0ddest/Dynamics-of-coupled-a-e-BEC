%% �������� ������
% ���������
hbar    = 1.05457172647e-34;    % ���������� ������ (�������������)

% ��������� ������
chi     = hbar * 7.6e-7;
Gamma1  = 1.629e-23;
Gamma2  = 304.4;
beta1   = 2.108e7;
beta2   = 3.344e6;
alpha   = 134.06;
omega   = 2 * pi * 100;         	% Trap frequency
N       = 5e5;                  	% Initial number of atoms

m       = 85;		  				% atomic mass???
a       = 17e-9;                	% s-wave scattering length (����� ����������� s-�����)
aHO     = sqrt(hbar/(2*m*omega));	% see Ref 23

alpha1  = alpha / omega;
beta11  = beta1 / omega;
beta21  = beta2 / omega;

Gamma11 = Gamma1 / (omega * aHO ^ 2);
Gamma21 = Gamma2 / omega;

lambda_a    = 4 * pi * hbar .^ 2 * a / m;	% ���� �������������� ����-����
lambda_m    = 4 * pi * hbar .^ 2 * a / m;	% ��������-��������
lambda_am   = 4 * pi * hbar .^ 2 * a / m;	% ����-��������

lambda_a1   = lambda_a / (aHO .^ 2 * hbar * omega);
lambda_m1   = lambda_m / (aHO .^ 2 * hbar * omega);
lambda_am1  = lambda_am / (aHO .^ 2 * hbar * omega);

epsilon     = 1;    % o_O WTF?
epsilon1    = epsilon / (hbar * omega);

% ��������� ������
h       = 1e-3;                         % ��� �� ������������
dt      = 1e-4;                         % ��� �� �������

xmax    = 15;                           % ������������ ���������� �� ������������
T       = 1e-3;                         % �������� ������ �������

Jmax    = double(floor(xmax / h));      % ����� ����� �� ������������
kmax    = floor(T / dt);                % ����� ����� �� �������

x       = linspace(0, xmax, Jmax);      % ��������� �������� x � ����� �����

% ��������� �������
phim    = zeros(kmax, Jmax);            % \phi_m (���� �� ������� � � �������) -- ����� ����.
phia    = zeros(kmax, Jmax);            % \phi_a (���� �� ������� � � �������) -- ��������� ���������� (�������� ���������� �� �������).

%phim(1, :) = sin((1:Jmax) / Jmax * pi);
phia(1, :) = gaussmf(x, [xmax/12 0]);%sin((1:Jmax) / Jmax * pi);

%% ������
for k = 1:kmax-1                        % ��� �� �������
    PHIAa   = sparse(Jmax, Jmax);       % ������������ ��� ���������� phia
    PHIAb   = zeros(Jmax, 1);           % ��������� ����� ��� ���������� phia

    PHIMa   = sparse(Jmax, Jmax);       % ������������ ��� ���������� phim
    PHIMb   = zeros(Jmax, 1);           % ��������� ����� ��� ���������� phim

    for j = 2:Jmax-1                    % ��� �� ������������ ��� phia
        
        % O_o WTF???
        PHIAa(1, 1)       = 1;                  % ����� ��������� �������
        PHIAb(1)          = 0;                  %   (x=0)
        PHIAa(Jmax, Jmax) = 1;                  % ����� ��������� �������
        PHIAb(Jmax)       = exp(-x(Jmax)^2/4);  %   (x=?)
        
        % ������� ���������� ������ ��� phia
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
        PHIMa(1, 1)       = 1;                  % ����� ��������� �������
        PHIMb(1)          = 0;                  %   (x=0)
        PHIMa(Jmax, Jmax) = 1;                  % ����� ��������� �������
        PHIMb(Jmax)        = exp(-x(Jmax)^2/4);  %   (x=?)
         
        % ������� ���������� ������ ��� phim
        Zb  = x(j) ^ 2 + epsilon1 + lambda_m1 * abs(phim(k, j)) ^ 2 / x(j) ^ 2 +...
            lambda_am1 * abs(phia(k, j)) ^ 2 / x(j) ^ 2 - beta21 - 1i * Gamma21;
        
        PHIMa(j, j-1)   = -1 / (8 * h ^ 2);
        PHIMa(j, j+1)   = -1 / (8 * h ^ 2);
        PHIMa(j, j)     = -1i / dt + 1 / (4 * h ^ 2) + Zb / 2;
        PHIMb(j)        = -1i / dt + (phim(k, j+1) - 2 * phim(k, j) + phim(k, j-1)) / (8 * h ^ 2) -...
            Zb * phim(k, j) / 2 - alpha1 / 2 * phia(k, j) ^ 2 / x(j);
    end
    phia(k+1, :) = (PHIAa \ PHIAb).';   % ������ ��������������� �������
    phim(k+1, :) = (PHIMa \ PHIMb).';   % ������ ��������������� �������
end

hold on
plot(abs(phia(1, :)), 'r');
plot(abs(phim(1, :)), 'b');
title('k=1');

figure
hold on
plot(abs(phia(k+1, :)), 'r');
plot(abs(phim(2, :)), 'b');
title('k=2');
