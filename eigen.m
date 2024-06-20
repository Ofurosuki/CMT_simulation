% 固定パラメータの設定
Gamma_plus = 0; % 自然減衰率の和 (定数)
Omega_plus = 0; % 角周波数の和 (定数)
kappa = 1.0; % 結合係数

% 変数範囲の設定
Omega_minus_range = linspace(-4.0, 4.0, 100); % 角周波数の差の範囲
Gamma_minus_range = linspace(0, 4.0, 100); % 自然減衰率の差の範囲

% 角周波数の差と自然減衰率の差を1/2にしたもの
Omega_minus_half_range = Omega_minus_range / 2;
Gamma_minus_half_range = Gamma_minus_range / 2;

% 固有値の実部を格納するための配列
lambda1_real_part = zeros(length(Gamma_minus_half_range), length(Omega_minus_half_range));
lambda2_real_part = zeros(length(Gamma_minus_half_range), length(Omega_minus_half_range));

% 固有値の虚部を格納するための配列
lambda1_imag_part = zeros(length(Gamma_minus_half_range), length(Omega_minus_half_range));
lambda2_imag_part = zeros(length(Gamma_minus_half_range), length(Omega_minus_half_range));

% 固有値の計算
for i = 1:length(Gamma_minus_half_range)
    for j = 1:length(Omega_minus_half_range)
        Gamma_minus_half = Gamma_minus_half_range(i);
        Omega_minus_half = Omega_minus_half_range(j);
        
        % ルートの中身を計算
        inside_sqrt = (2*Gamma_minus_half)^2 - (2*Omega_minus_half)^2 - 4*kappa^2 - 4*1i*(2*Gamma_minus_half)*(2*Omega_minus_half);
        
        % 固有値の計算
        sqrt_term = sqrt(inside_sqrt);
        lambda1 = (Gamma_plus/2) - 1i*(Omega_plus/2) + (sqrt_term/2);
        lambda2 = (Gamma_plus/2) - 1i*(Omega_plus/2) - (sqrt_term/2);
        
        % 実部を保存
        lambda1_real_part(i, j) = real(lambda1); % lambda1の実部を使用
        lambda2_real_part(i, j) = real(lambda2); % lambda2の実部を使用

        % 虚部を保存
        lambda1_imag_part(i, j) = imag(lambda1); % lambda1の虚部を使用
        lambda2_imag_part(i, j) = imag(lambda2); % lambda2の虚部を使用
    end
end

% メッシュグリッドの作成
[Omega_minus_half_mesh, Gamma_minus_half_mesh] = meshgrid(Omega_minus_half_range, Gamma_minus_half_range);

% 実部のグラフ
figure;

% プロット lambda1の実部
surf1 = surf(Omega_minus_half_mesh, Gamma_minus_half_mesh, lambda1_real_part);
hold on;

% プロット lambda2の実部
surf2 = surf(Omega_minus_half_mesh, Gamma_minus_half_mesh, lambda2_real_part);

% カラーマップの設定
colormap('parula');

% グラフの装飾
xlabel('Omega diffirence (rad/s)');
ylabel('Gamma diffirence');
zlabel('Real part of Eigenvalue');
title('Real Part of Eigenvalues');
colorbar;

% 凡例の追加
legend([surf1, surf2], {'\lambda_1 Real Part', '\lambda_2 Real Part'}, 'Color', 'white');
hold off;

% 虚部のグラフ
figure;

% プロット lambda1の虚部
surf3 = surf(Omega_minus_half_mesh, Gamma_minus_half_mesh, lambda1_imag_part);
hold on;

% プロット lambda2の虚部
surf4 = surf(Omega_minus_half_mesh, Gamma_minus_half_mesh, lambda2_imag_part);

% カラーマップの設定
colormap('parula');

% グラフの装飾
xlabel('Omega diffirence (rad/s)');
ylabel('Gamma difference');
zlabel('Imaginary part of Eigenvalue');
title('Imaginary Part oigenvalues');
colorbar;

% 凡例の追加
legend([surf3, surf4], {'\lambda_1 Imaginary Part', '\lambda_2 Imaginary Part'}, 'Color', 'white');
hold off;
