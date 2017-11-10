%% Parameter Setting
n_k = 400;  % K step
n_mu = 400; % Mu step
e_min = -pi; % Mu limit
e_max = pi;
t = 1;
T = 0.08 * t;
%% Hamitonian Setting 
syms kx_sym ky_sym kz_sym;
f_sym(1) = t * ( cos( ky_sym ) + 2 * cos( ky_sym / 2 ) *  cos( sqrt(3) * kx_sym / 2 ) );
f_sym(2) = t * ( -sin( ky_sym ) + 2 * sin( ky_sym / 2 ) *  cos( sqrt(3) * kx_sym / 2 ) );
MF = matlabFunction([f_sym(1);f_sym(2)]);
%% No Necessity to Modify blow// Automatical Calculation
%% Calculation Limit
mu_min = e_min;
mu_max = e_max;
mu = linspace(e_min, e_max, n_mu);
kx = linspace(-pi, pi, n_k);
ky = linspace(-pi, pi, n_k);
%% Calculation on Derivation 
for i = 1:2
    fx_sym(i) = diff(f_sym(i),kx_sym);
    fxx_sym(i) = diff(fx_sym(i),kx_sym);
    fy_sym(i) = diff(f_sym(i),ky_sym);
    fyy_sym(i) = diff(fy_sym(i),ky_sym);
    fxy_sym(i) = diff(fx_sym(i),ky_sym);
end
MFx = matlabFunction([fx_sym(1);fx_sym(2)]);
MFy = matlabFunction([fy_sym(1);fy_sym(2)]);
MFxx = matlabFunction([fxx_sym(1);fxx_sym(2)]);
MFyy = matlabFunction([fyy_sym(1);fyy_sym(2)]);
MFxy = matlabFunction([fyy_sym(1);fyy_sym(2)]);
for i = 1:n_mu
    chi_spin(i) = 0;
    chi_LP(i) = 0;
    chi_omega_a(i) = 0;
    chi_omega_b(i) = 0;
    chi_g(i) = 0;
    chi_tildeg_a(i) = 0;
    chi_tildeg_b(i) = 0;
    chi_tildeg_c(i) = 0;
end
%%
tic
for i = 1:n_k
    for j = 1:n_k
        % Numerical Application 
        f = MF(kx(i),ky(j));
        fx = MFx(kx(i),ky(j));
        fy = MFy(kx(i),ky(j));
        fxx = MFxx(kx(i),ky(j));
        fyy = MFyy(kx(i),ky(j));
        fxy = MFxy(kx(i),ky(j));
        % Temporary Data
        [chi_spin_T, chi_LP_T, chi_omega_a_T, chi_omega_b_T, chi_g_T, chi_tildeg_a_T, chi_tildeg_b_T, chi_tildeg_c_T] = MS(f,fx,fy,fxx,fyy,fxy,mu,T); 
        % Summation
        chi_spin = chi_spin + chi_spin_T;
        chi_LP = chi_LP + chi_LP_T;
        chi_omega_a = chi_omega_a + chi_omega_a_T;
        chi_omega_b = chi_omega_b + chi_omega_b_T;
        chi_g = chi_g + chi_g_T;
        chi_tildeg_a = chi_tildeg_a + chi_tildeg_a_T;
        chi_tildeg_b = chi_tildeg_b + chi_tildeg_b_T;
        chi_tildeg_c = chi_tildeg_c + chi_tildeg_c_T;
    end
end
toc
%%
    c_spin = chi_spin / ( n_k * n_k );
    c_LP = chi_LP / ( n_k * n_k );
    c_omega_a = chi_omega_a / ( n_k * n_k );
    c_omega_b = chi_omega_b / ( n_k * n_k );
    c_g = chi_g / ( n_k * n_k );
    c_tildeg_a = chi_tildeg_a / ( n_k * n_k );
    c_tildeg_b = chi_tildeg_b / ( n_k * n_k );
    c_tildeg_c = chi_tildeg_c / ( n_k * n_k );
    c_orb = c_LP + c_omega_a + c_omega_b + c_g + c_tildeg_a + c_tildeg_c;
%% Plot
figure

subplot(2,2,1)
plot(mu,c_orb,'Color',[237,3,69]/255,'LineWidth',2);
box on
set(gca,'FontSize',10,'Fontname','Cambria Math');
set(gca, 'XMinorTick', 'on');
set(gca, 'YMinorTick', 'on');
set(gca, 'XGrid', 'on');
set(gca, 'YGrid', 'on');
xlabel('$\mu/t$','interpret','latex','FontSize',16);
ylabel('$\chi_{orb}/\chi_0$','interpret','latex','FontSize',16);
subplot(2,2,2)

plot(mu,c_LP,'Color',[3,195,131]/255,'LineWidth',2);
box on
set(gca,'FontSize',10,'Fontname','Cambria Math');
set(gca, 'XMinorTick', 'on');
set(gca, 'YMinorTick', 'on');
set(gca, 'XGrid', 'on');
set(gca, 'YGrid', 'on');
xlabel('$\mu/t$','interpret','latex','FontSize',16);
ylabel('$\chi_{LP}/\chi_0$','interpret','latex','FontSize',16);

subplot(2,2,3)
plot(mu,c_g,'Color',[251,191,69]/255,'LineWidth',2);
box on
set(gca,'FontSize',10,'Fontname','Cambria Math');
set(gca, 'XMinorTick', 'on');
set(gca, 'YMinorTick', 'on');
set(gca, 'XGrid', 'on');
set(gca, 'YGrid', 'on');
xlabel('$\mu/t$','interpret','latex','FontSize',16);
ylabel('$\chi_{g}/\chi_0$','interpret','latex','FontSize',16);

subplot(2,2,4)
plot(mu,c_omega_a + c_omega_b,'Color',[17,1,65]/255,'LineWidth',2);
box on
set(gca,'FontSize',10,'Fontname','Cambria Math');
set(gca, 'XMinorTick', 'on');
set(gca, 'YMinorTick', 'on');
set(gca, 'XGrid', 'on');
set(gca, 'YGrid', 'on');
xlabel('$\mu/t$','interpret','latex','FontSize',16);
ylabel('$\chi_{\Omega}/\chi_0$','interpret','latex','FontSize',16);
