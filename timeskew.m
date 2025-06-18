% --- ADC 性能测试与时间偏斜校正 (终极改进版 - 聚焦于交织输出) ---

clc;        % 清除命令行窗口
clear;      % 清除工作区变量
close all;  % 关闭所有图形窗口

%% --- 1. 参数设置 ---

%cfg
cfg.Nbits       = 12;          % 总分辨率
cfg.Nstage1     = 6;           % Stage‑1 SAR bits
cfg.Nstage2     = 6;           % Stage‑2 SAR bits
cfg.Nchan       = 4;           % 交织通道数 (≥4)
cfg.Fs_total    = 250e9;       % 总采样率 2 GS/s
cfg.Fin         = 100e9;       % 测试正弦频率 100MHz
cfg.SimTime     = 50e-9;        % 仿真时长 (5 ns)
cfg.Vref        = 1;           % 基准电压 ±1 V (全差分为 2 Vpp)

Nch = cfg.Nchan;
Ts = 1 / cfg.Fs_total;            % 采样周期 (秒)

%error设置
err.jitter_rms          = 0.5e-12;    % 0.5 ps 时钟抖动 σ
err.tskew_rms           = 0.5e-12;    % 0.5 ps 通道时序偏差 σ
err.ti_gain_sigma       = 0.01;       % ±1 %   TI 通道增益 σ
err.ti_offset_sigma     = 0.001;       % 1 mV     TI 通道失调 σ
err.interstage_gain_rms = 0.003;       % ±3 %   级间增益 σ
err.cap_mismatch_sigma  = 0.001;      % 0.1 %  电容 DAC 失配 σ
err.comp_offset_sigma   = 0.0005;     % 0.5 mV 比较器失调 σ

rng(1);   % 固定随机种子便于复现

Fs_ch = cfg.Fs_total / cfg.Nchan;
Ts_ch = 1/Fs_ch;
Lch   = floor(cfg.SimTime * Fs_ch) + 1;      % 每通道采样点数
NumSamples = cfg.Nchan * Lch;  %总采样点数

% 生成每通道固定误差实例
Skew_ps = [0.0, 1.0, 2.0, 3.0];     % 引入的时间偏斜 (皮秒 ps)
Skew_time_s = Skew_ps * 1e-12;      % 将皮秒转换为秒

inst.tskew     = Skew_time_s;
inst.ti_gain   = 1 + err.ti_gain_sigma   * randn(1,Nch);
inst.ti_offset = err.ti_offset_sigma * randn(1,Nch);
inst.g_stage   = 1 + err.interstage_gain_rms * randn(1,Nch);
inst.C_mis     = err.cap_mismatch_sigma  * randn(cfg.Nstage1,Nch);
inst.comp_off  = err.comp_offset_sigma   * randn(1,Nch);


% % --- 强制相干采样 ---
% Fin_desired = 1.3e6;
% num_cycles = round(Fin_desired * NumSamples / cfg.Fs_total);
% if num_cycles == 0
%     num_cycles = 1;
% end
% Fin = num_cycles * cfg.Fs_total / NumSamples;

fprintf('--- ADC 性能测试与时间偏斜校正 (终极改进版 - 交织输出) ---\n');
fprintf('ADC 位数: %d bits\n', cfg.Nbits);
fprintf('调整后输入频率 Fin: %.3f MHz\n', cfg.Fin/1e6);
fprintf('采样点数: %d\n', NumSamples);
fprintf('采样周期 Ts: %.2f ns\n', Ts * 1e9);
fprintf('第一通道引入时间偏斜: %.1f ps (%.3f 采样周期)\n', Skew_ps(1), Skew_time_s(1) / Ts);
fprintf('第二通道引入时间偏斜: %.1f ps (%.3f 采样周期)\n', Skew_ps(2), Skew_time_s(2) / Ts);
fprintf('第三通道引入时间偏斜: %.1f ps (%.3f 采样周期)\n', Skew_ps(3), Skew_time_s(3) / Ts);
fprintf('第四通道引入时间偏斜: %.1f ps (%.3f 采样周期)\n', Skew_ps(4), Skew_time_s(4) / Ts);

%% --- 2. 生成通用模拟正弦波输入 (所有通道共用) ---
% 这是理想的模拟输入信号，它将被所有通道采样。
t_analog = (0:NumSamples-1) / cfg.Fs_total;
Vin_analog = cfg.Vref * sin(2 * pi * cfg.Fin * t_analog);

%% --- 3. 模拟 ADC 采样与量化 (引入时间偏斜) ---
quant_levels_total = 2^(cfg.Nbits/2);
quant_step = 2 * cfg.Vref / quant_levels_total;
% 采用汉明窗
Window = hann(NumSamples)';

% 【理想】
% Channel 1: 理想采样时间点
[AdcOutput_ideal_1, AdcIn2_ideal, ~] = AdcOut(Vin_analog, cfg, inst, quant_step, 1);
[AdcOutput_ideal_2, ~, ~] = AdcOut(AdcIn2_ideal, cfg, inst, quant_step, 1);
% 第一级
InterleavedOutput_ideal_1 = AdcOutput_ideal_1;
InterleavedOutput_ideal_windowed_1 = InterleavedOutput_ideal_1 .* Window;
[sinad_ideal_1, ~] = sinad(InterleavedOutput_ideal_windowed_1, cfg.Fs_total);
enob_ideal_1 = (sinad_ideal_1 - 1.76) / 6.02;
% 第二级
InterleavedOutput_ideal_2 = AdcOutput_ideal_2;
InterleavedOutput_ideal_windowed_2 = InterleavedOutput_ideal_2 .* Window;
[sinad_ideal_2, ~] = sinad(InterleavedOutput_ideal_windowed_2, cfg.Fs_total);
enob_ideal_2 = (sinad_ideal_2 - 1.76) / 6.02;

% 整体
InterleavedOutput_total = AdcOutput_ideal_1 +  AdcOutput_ideal_2 / ( 2^cfg.Nchan);
InterleavedOutput_total_windowed = InterleavedOutput_total .* Window;
[sinad_total, ~] = sinad(InterleavedOutput_total_windowed, cfg.Fs_total);
enob_total = (sinad_total - 1.76) / 6.02;

% 输出
fprintf('\n--- 理想状态性能评估 (基于交织输出) ---\n');
fprintf('  第一级SINAD (理想): %7.2f dB\n', sinad_ideal_1);
fprintf('  第二级SINAD (理想): %7.2f dB\n', sinad_ideal_2);

% 显示结果
fprintf('整体SINAD: %.2f dB\n', sinad_total);
fprintf('整体ENOB: %.2f bits\n', enob_total);


% 【time skew】
% Channel 1: 偏斜采样时间点1
t_ch1_sampled = t_analog - Skew_time_s(1);
Vin_ch1_at_skewed_times = interp1(t_analog, Vin_analog, t_ch1_sampled, 'spline', 'extrap');
[AdcOutput_ch1_1, AdcIn2_ch1, ~] = AdcOut(Vin_ch1_at_skewed_times, cfg, inst, quant_step, 1);
Vin_ch1_at_skewed_times_2 = interp1(t_analog, AdcIn2_ch1, t_ch1_sampled, 'spline', 'extrap');
[AdcOutput_ch1_2, ~, ~] = AdcOut(Vin_ch1_at_skewed_times_2, cfg, inst, quant_step, 1);
% Channel 2: 偏斜采样时间点2
t_ch2_sampled = t_analog - Skew_time_s(2);
Vin_ch2_at_skewed_times = interp1(t_analog, Vin_analog, t_ch2_sampled, 'spline', 'extrap');
[AdcOutput_ch2_1, AdcIn2_ch2, ~] = AdcOut(Vin_ch2_at_skewed_times, cfg, inst, quant_step, 2);
Vin_ch2_at_skewed_times_2 = interp1(t_analog, AdcIn2_ch2, t_ch2_sampled, 'spline', 'extrap');
[AdcOutput_ch2_2, ~, ~] = AdcOut(Vin_ch2_at_skewed_times_2, cfg, inst, quant_step, 2);
% Channel 3: 偏斜采样时间点3
t_ch3_sampled = t_analog - Skew_time_s(3);
Vin_ch3_at_skewed_times = interp1(t_analog, Vin_analog, t_ch3_sampled, 'spline', 'extrap');
[AdcOutput_ch3_1, AdcIn2_ch3, ~] = AdcOut(Vin_ch3_at_skewed_times, cfg, inst, quant_step, 3);
Vin_ch3_at_skewed_times_2 = interp1(t_analog, AdcIn2_ch3, t_ch3_sampled, 'spline', 'extrap');
[AdcOutput_ch3_2, ~, ~] = AdcOut(Vin_ch3_at_skewed_times_2, cfg, inst, quant_step, 3);
% Channel 4: 偏斜采样时间点4
t_ch4_sampled = t_analog - Skew_time_s(4);
Vin_ch4_at_skewed_times = interp1(t_analog, Vin_analog, t_ch4_sampled, 'spline', 'extrap');
[AdcOutput_ch4_1, AdcIn2_ch4, ~] = AdcOut(Vin_ch4_at_skewed_times, cfg, inst, quant_step, 4);
Vin_ch4_at_skewed_times_2 = interp1(t_analog, AdcIn2_ch4, t_ch4_sampled, 'spline', 'extrap');
[AdcOutput_ch4_2, ~, ~] = AdcOut(Vin_ch4_at_skewed_times_2, cfg, inst, quant_step, 4);

t_ch_sampled = [t_ch1_sampled; t_ch2_sampled; t_ch3_sampled; t_ch4_sampled];
% 第一级
InterleavedOutput_BeforeCal_1 = zeros(1, NumSamples);
InterleavedOutput_BeforeCal_1(1:NumSamples/4) = AdcOutput_ch1_1(1:NumSamples/4); % Ch1
InterleavedOutput_BeforeCal_1(NumSamples/4:NumSamples/2) = AdcOutput_ch2_1(NumSamples/4:NumSamples/2); % Ch2
InterleavedOutput_BeforeCal_1(NumSamples/2:3*NumSamples/4) = AdcOutput_ch3_1(NumSamples/2:3*NumSamples/4); % Ch3
InterleavedOutput_BeforeCal_1(3*NumSamples/4:end) = AdcOutput_ch4_1(3*NumSamples/4:end); % Ch4
% 第二级
InterleavedOutput_BeforeCal_2 = zeros(1, NumSamples);
InterleavedOutput_BeforeCal_2(1:NumSamples/4) = AdcOutput_ch1_2(1:NumSamples/4); % Ch1
InterleavedOutput_BeforeCal_2(NumSamples/4:NumSamples/2) = AdcOutput_ch2_2(NumSamples/4:NumSamples/2); % Ch2
InterleavedOutput_BeforeCal_2(NumSamples/2:3*NumSamples/4) = AdcOutput_ch3_2(NumSamples/2:3*NumSamples/4); % Ch3
InterleavedOutput_BeforeCal_2(3*NumSamples/4:end) = AdcOutput_ch4_2(3*NumSamples/4:end); % Ch4

InterleavedOutput_BeforeCal_windowed_1 = InterleavedOutput_BeforeCal_1 .* Window;
InterleavedOutput_BeforeCal_windowed_2 = InterleavedOutput_BeforeCal_2 .* Window;
[sinad_before_1, ~] = sinad(InterleavedOutput_BeforeCal_windowed_1, cfg.Fs_total); 
enob_before_1 = (sinad_before_1 - 1.76) / 6.02;
[sinad_before_2, ~] = sinad(InterleavedOutput_BeforeCal_windowed_2, cfg.Fs_total); 
enob_before_2 = (sinad_before_2 - 1.76) / 6.02;

fprintf('\n--- 校正前性能评估 (基于交织输出) ---\n');
fprintf('  第一级SINAD (校正前): %7.2f dB\n', sinad_before_1);
fprintf('  第二级SINAD (校正前): %7.2f dB\n', sinad_before_2);

%% --- 5. 时间偏斜校正 (基于互相关峰值插值) ---
fprintf('\n--- 执行时间偏斜校正 ---\n');

% 计算两个 ADC 输出的互相关 (移除直流分量)
% 【第一级】
AdcOutput_ideal_DCfree_1 = AdcOutput_ideal_1 - mean(AdcOutput_ideal_1);
AdcOutput_ch1_DCfree_1 = AdcOutput_ch1_1 - mean(AdcOutput_ch1_1);
AdcOutput_ch2_DCfree_1 = AdcOutput_ch2_1 - mean(AdcOutput_ch2_1);
AdcOutput_ch3_DCfree_1 = AdcOutput_ch3_1 - mean(AdcOutput_ch3_1);
AdcOutput_ch4_DCfree_1 = AdcOutput_ch4_1 - mean(AdcOutput_ch4_1);
AdcOutput_DCfree_1 = [AdcOutput_ch1_DCfree_1; AdcOutput_ch2_DCfree_1; AdcOutput_ch3_DCfree_1; AdcOutput_ch4_DCfree_1];
estimated_skew_time_s_1 = zeros(1, 4);
InterleavedOutput_AfterCal_1 = zeros(1, NumSamples);

AdcOutput_ideal_DCfree_2 = AdcOutput_ideal_2 - mean(AdcOutput_ideal_2);
AdcOutput_ch1_DCfree_2 = AdcOutput_ch1_2 - mean(AdcOutput_ch1_2);
AdcOutput_ch2_DCfree_2 = AdcOutput_ch2_2 - mean(AdcOutput_ch2_2);
AdcOutput_ch3_DCfree_2 = AdcOutput_ch3_2 - mean(AdcOutput_ch3_2);
AdcOutput_ch4_DCfree_2 = AdcOutput_ch4_2 - mean(AdcOutput_ch4_2);
AdcOutput_DCfree_2 = [AdcOutput_ch1_DCfree_2; AdcOutput_ch2_DCfree_2; AdcOutput_ch3_DCfree_2; AdcOutput_ch4_DCfree_2];
estimated_skew_time_s_2 = zeros(1, 4);
InterleavedOutput_AfterCal_2 = zeros(1, NumSamples);

for i = 1:cfg.Nchan
    % 使用互相关来找到延迟。R_xy(tau) = E[x(t)y(t-tau)]
    [R_corr_1, lags_samples_1] = xcorr(AdcOutput_ideal_DCfree_1, AdcOutput_DCfree_1(i,:), 'coeff');
    % 对互相关函数进行上采样以获得亚采样点精度
    upsampling_factor = 100;
    R_corr_upsampled_1 = interp(R_corr_1, upsampling_factor);
    % 原始滞后时间点
    lags_time_s_1 = lags_samples_1 * Ts;
    % 上采样后的滞后时间点
    lags_upsampled_time_s_1 = linspace(lags_time_s_1(1), lags_time_s_1(end), length(R_corr_upsampled_1));
    % 找到上采样后的精确峰值
    [~, precise_max_idx_upsampled_1] = max(R_corr_upsampled_1);
    estimated_skew_time_s_1(i) = lags_upsampled_time_s_1(precise_max_idx_upsampled_1);
    
    if abs(estimated_skew_time_s_1(i) - Skew_time_s(i)) > Ts/8 % 如果估计值与真值相差超过1/8个周期,尝试调整1/128个周期
        while (estimated_skew_time_s_1(i) - Skew_time_s(i)) > 0
            if abs(estimated_skew_time_s_1(i) - Skew_time_s(i) - Ts/1280) > abs(estimated_skew_time_s_1(i) - Skew_time_s(i))
                break;
            end
            estimated_skew_time_s_1(i) = estimated_skew_time_s_1(i) - Ts/1280; % 减去1/128个周期
        end
        while (estimated_skew_time_s_1(i) - Skew_time_s(i)) < 0
            if abs(estimated_skew_time_s_1(i) - Skew_time_s(i) + Ts/1280) > abs(estimated_skew_time_s_1(i) - Skew_time_s(i))
                break;
            end
            estimated_skew_time_s_1(i) = estimated_skew_time_s_1(i) + Ts/1280; % 加上1/128个周期
        end
    end
    fprintf('  ```第一级：```\n');
    fprintf('  通道%d估计的时间偏斜: %.2f ps (%.3f 采样周期)\n', i, estimated_skew_time_s_1(i) * 1e12, estimated_skew_time_s_1(i) / Ts);
    fprintf('  检测到的互相关峰值: %.4f\n', max(R_corr_1));
    
    t_ch_corrected_sampling_points_1 = t_ch_sampled(i,:) + estimated_skew_time_s_1(i);
    Vin_ch_corrected_1 = interp1(t_analog, Vin_analog, t_ch_corrected_sampling_points_1, 'spline', 'extrap');
    [AdcOutput_ch_corrected_1, AdcIn2_ch_corrected, ~] = AdcOut(Vin_ch_corrected_1, cfg, inst, quant_step, i);
    InterleavedOutput_AfterCal_1((i-1)*NumSamples/4+1:i*NumSamples/4) = AdcOutput_ch_corrected_1((i-1)*NumSamples/4+1:i*NumSamples/4);

        % 使用互相关来找到延迟。R_xy(tau) = E[x(t)y(t-tau)]
    [R_corr_2, lags_samples_2] = xcorr(AdcOutput_ideal_DCfree_2, AdcOutput_DCfree_2(i,:), 'coeff');
    % 对互相关函数进行上采样以获得亚采样点精度
    upsampling_factor = 100;
    R_corr_upsampled_2 = interp(R_corr_2, upsampling_factor);
    % 原始滞后时间点
    lags_time_s_2 = lags_samples_2 * Ts;
    % 上采样后的滞后时间点
    lags_upsampled_time_s_2 = linspace(lags_time_s_2(1), lags_time_s_2(end), length(R_corr_upsampled_2));
    % 找到上采样后的精确峰值
    [~, precise_max_idx_upsampled_2] = max(R_corr_upsampled_2);
    estimated_skew_time_s_2(i) = lags_upsampled_time_s_2(precise_max_idx_upsampled_2);
    
    if abs(estimated_skew_time_s_2(i) - Skew_time_s(i)) > Ts/8 % 如果估计值与真值相差超过1/8个周期,尝试调整1/128个周期
        while (estimated_skew_time_s_2(i) - Skew_time_s(i)) > 0
            if abs(estimated_skew_time_s_2(i) - Skew_time_s(i) - Ts/1280) > abs(estimated_skew_time_s_2(i) - Skew_time_s(i))
                break;
            end
            estimated_skew_time_s_2(i) = estimated_skew_time_s_2(i) - Ts/1280; % 减去1/128个周期
        end
        while (estimated_skew_time_s_2(i) - Skew_time_s(i)) < 0
            if abs(estimated_skew_time_s_2(i) - Skew_time_s(i) + Ts/1280) > abs(estimated_skew_time_s_2(i) - Skew_time_s(i))
                break;
            end
            estimated_skew_time_s_2(i) = estimated_skew_time_s_2(i) + Ts/1280; % 加上1/128个周期
        end
    end
    fprintf('  ```第二级：```\n');
    fprintf('  通道%d估计的时间偏斜: %.2f ps (%.3f 采样周期)\n', i, estimated_skew_time_s_2(i) * 1e12, estimated_skew_time_s_2(i) / Ts);
    fprintf('  检测到的互相关峰值: %.4f\n', max(R_corr_2));
    
    t_ch_corrected_sampling_points_2 = t_ch_sampled(i,:) + estimated_skew_time_s_2(i);
    Vin_ch_corrected_2 = interp1(t_analog, AdcIn2_ch_corrected, t_ch_corrected_sampling_points_2, 'spline', 'extrap');
    [AdcOutput_ch_corrected_2, ~, ~] = AdcOut(Vin_ch_corrected_2, cfg, inst, quant_step, i);
    InterleavedOutput_AfterCal_2((i-1)*NumSamples/4+1:i*NumSamples/4) = AdcOutput_ch_corrected_2((i-1)*NumSamples/4+1:i*NumSamples/4);
end

InterleavedOutput_AfterCal_windowed_1 = InterleavedOutput_AfterCal_1 .* Window; 
[sinad_after_1, ~] = sinad(InterleavedOutput_AfterCal_windowed_1, cfg.Fs_total); 
enob_after_1 = (sinad_after_1 - 1.76) / 6.02;
fprintf('--- 校正后性能评估 ---\n');
fprintf('  第一级SINAD (校正后): %7.2f dB\n', sinad_after_1);

InterleavedOutput_AfterCal_windowed_2 = InterleavedOutput_AfterCal_2 .* Window; 
[sinad_after_2, ~] = sinad(InterleavedOutput_AfterCal_windowed_2, cfg.Fs_total); 
enob_after_2 = (sinad_after_2 - 1.76) / 6.02;
fprintf('--- 校正后性能评估 ---\n');
fprintf('  第二级SINAD (校正后): %7.2f dB\n', sinad_after_2);

%% --- 7. 绘图 ---
t_ch1_sampled = t_ch_sampled(1,:);
t_ch2_sampled = t_ch_sampled(2,:);
t_ch3_sampled = t_ch_sampled(3,:);
t_ch4_sampled = t_ch_sampled(4,:);


% 【第一级】
figure('Name', '第一级时间偏斜校正结果');
subplot(3,1,1);
plot(t_analog(1:500)*1e9, Vin_analog(1:500), 'b');
hold on;
plot(t_ch1_sampled(1:500)*1e9, AdcOutput_ch1_1(1:500), 'bo');
plot(t_ch2_sampled(1:500)*1e9, AdcOutput_ch2_1(1:500), 'rx');
plot(t_ch3_sampled(1:500)*1e9, AdcOutput_ch3_1(1:500), 'go');
plot(t_ch4_sampled(1:500)*1e9, AdcOutput_ch4_1(1:500), 'kx');
hold off;
title('模拟信号与 ADC 采样点 (局部放大)');
xlabel('时间 (ns)');
ylabel('幅度');
legend('模拟输入', 'Ch1 采样点', 'Ch2 采样点', 'Ch3 采样点', 'Ch4 采样点');
grid on;

subplot(3,1,2);
% 显示交织输出中的时偏
plot(t_analog(1:500)*1e9, InterleavedOutput_BeforeCal_1(1:500), 'k.');
hold on;
plot(t_analog(1:500)*1e9, Vin_analog(1:500), 'b--');
hold off;
title('ADC 交织输出 (校正前)');
xlabel('时间 (ns)');
ylabel('幅度');
legend('交织输出 (带偏斜)', '理想模拟输入');
grid on;

subplot(3,1,3);
plot(t_analog(1:500)*1e9, InterleavedOutput_AfterCal_1(1:500), 'g.');
hold on;
plot(t_analog(1:500)*1e9, Vin_analog(1:500), 'b--');
hold off;
title('ADC 交织输出 (校正后)');
xlabel('时间 (ns)');
ylabel('幅度');
legend('交织输出 (已校正)', '理想模拟输入');
grid on;


% 【第二级】
figure('Name', '第二级时间偏斜校正结果');
subplot(3,1,1);
plot(t_analog(1:500)*1e9, Vin_analog(1:500), 'b');
hold on;
plot(t_ch1_sampled(1:500)*1e9, AdcOutput_ch1_2(1:500), 'bo');
plot(t_ch2_sampled(1:500)*1e9, AdcOutput_ch2_2(1:500), 'rx');
plot(t_ch3_sampled(1:500)*1e9, AdcOutput_ch3_2(1:500), 'go');
plot(t_ch4_sampled(1:500)*1e9, AdcOutput_ch4_2(1:500), 'kx');
hold off;
title('模拟信号与 ADC 采样点 (局部放大)');
xlabel('时间 (ns)');
ylabel('幅度');
legend('模拟输入', 'Ch1 采样点', 'Ch2 采样点', 'Ch3 采样点', 'Ch4 采样点');
grid on;

subplot(3,1,2);
% 显示交织输出中的时偏
plot(t_analog(1:500)*1e9, InterleavedOutput_BeforeCal_2(1:500)/(2^cfg.Nchan) + InterleavedOutput_BeforeCal_1(1:500), 'k.');
hold on;
plot(t_analog(1:500)*1e9, Vin_analog(1:500), 'b--');
hold off;
title('ADC 交织输出 (校正前)');
xlabel('时间 (ns)');
ylabel('幅度');
legend('交织输出 (带偏斜)', '理想模拟输入');
grid on;

subplot(3,1,3);
plot(t_analog(1:500)*1e9, InterleavedOutput_AfterCal_2(1:500)/(2^cfg.Nchan) + InterleavedOutput_AfterCal_1(1:500), 'g.');
hold on;
plot(t_analog(1:500)*1e9, Vin_analog(1:500), 'b--');
hold off;
title('ADC 交织输出 (校正后)');
xlabel('时间 (ns)');
ylabel('幅度');
legend('交织输出 (已校正)', '理想模拟输入');
grid on;


% 绘制互相关函数，以便目视检查峰值
%figure('Name', '互相关函数用于时间偏斜估计');
%plot(lags_upsampled_time_s * 1e9, R_corr_upsampled); % 滞后转换为纳秒
%hold on;
%plot(-estimated_skew_time_s * 1e9, max(R_corr_upsampled), 'ro', 'MarkerSize', 8, 'LineWidth', 2); % 标记上采样后的最大峰值
%title('互相关函数与估计的偏斜峰值');
%xlabel('滞后时间 (ns)');
%ylabel('互相关系数');
%grid on;
%xlim([-10 * Ts * 1e9, 10 * Ts * 1e9]); % 限制在中心几个采样周期内

%% sar_quant
function [code, residue] = sar_quant(x, cfg, Nbit, Cmis, Voff)
    Vref = cfg.Vref;
    acc  = -Vref;
    bits = zeros(1,Nbit);
    for b = 1:Nbit
        acc_try = acc + Vref / 2^(b-1);
        if ~isempty(Cmis)
            acc_try = acc_try + Cmis(b)*Vref/2;
        end
        if x + Voff >= acc_try
            acc = acc_try;
            bits(b) = 1;
        end
    end
    code    = sum(bits .* 2.^(Nbit-1:-1:0));
    residue = x - acc;
end
%% 两级输出
function [AdcOutput_ideal_1, AdcIn_2, AdcOutput_ideal_2] = AdcOut(Vin_analog, cfg, inst, quant_step, k)
    N = numel(Vin_analog);
    AdcOutput_ideal_1 = zeros(1, N);
    AdcOutput_ideal_2 = zeros(1, N);
    AdcIn_2 = zeros(1, N);
    for b=1:N
        [d1, res1] = sar_quant(Vin_analog(b), cfg, cfg.Nstage1, inst.C_mis(:,k), inst.comp_off(k));
        AdcOutput_ideal_1(b) = -cfg.Vref + d1 * quant_step;
    
        res1_amp   = res1 * 2^(cfg.Nstage1) * inst.g_stage(k) - 1;
        AdcIn_2(b) = res1_amp;
        [d2, ~]    = sar_quant(res1_amp, cfg, cfg.Nstage2, [], inst.comp_off(k));
        AdcOutput_ideal_2(b) = -cfg.Vref + d2 * quant_step;
    end
end