
clear; clc; close all;

%cfg
cfg.Nbits       = 12;          % 总分辨率
cfg.Nstage1     = 6;           % Stage‑1 SAR bits
cfg.Nstage2     = 6;           % Stage‑2 SAR bits
cfg.Nchan       = 4;           % 交织通道数 (≥4)
cfg.Fs_total    = 2e9;         % 总采样率 2 GS/s
cfg.Fin         = 100e6;       % 测试正弦频率
cfg.SimTime     = 5e-6;        % 仿真时长 (5 µs)
cfg.Vref        = 1;           % 基准电压 ±1 V (全差分为 2 Vpp)

%error设置
err.jitter_rms          = 0.5e-12;    % 0.5 ps 时钟抖动 σ
err.tskew_rms           = 0.5e-12;    % 0.5 ps 通道时序偏差 σ
err.ti_gain_sigma       = 0.01;       % ±1 %   TI 通道增益 σ
err.ti_offset_sigma     = 0.001;       % 1 mV     TI 通道失调 σ
err.interstage_gain_rms = 0.003;       % ±3 %   级间增益 σ
err.cap_mismatch_sigma  = 0.001;      % 0.1 %  电容 DAC 失配 σ
err.comp_offset_sigma   = 0.0005;     % 0.5 mV 比较器失调 σ

rng(1);   % 固定随机种子便于复现



[out_dig, meta] = TI_PipeSAR_core(cfg, err);


[enob, sinad, sfdr, thd] = dynamic_metrics(out_dig, cfg.Nbits, ...
                                           cfg.Fs_total, cfg.Fin);

fprintf('\n=== Dynamic metrics (all errors on) ===\n');
fprintf('   SINAD : %6.2f  dB\n', sinad);
fprintf('   ENOB  : %6.2f  bits\n', enob);
fprintf('   SFDR  : %6.2f  dB\n', sfdr);
fprintf('   THD   : %6.2f  dB\n', thd);




figure;  plot_fft_dbfs(out_dig, cfg.Fs_total);
title('TI‑Pipelined‑SAR ADC Output Spectrum');


%顶层
function [dig_out, meta] = TI_PipeSAR_core(cfg, err)


    Nch   = cfg.Nchan;
    Fs_ch = cfg.Fs_total / Nch;
    Ts_ch = 1/Fs_ch;
    Lch   = floor(cfg.SimTime * Fs_ch) + 1;      % 每通道采样点数
    
    % 生成每通道固定误差实例
    inst.tskew     = err.tskew_rms           * randn(1,Nch);
    inst.ti_gain   = 1 + err.ti_gain_sigma   * randn(1,Nch);
    inst.ti_offset =     err.ti_offset_sigma * randn(1,Nch);
    inst.g_stage   = 1 + err.interstage_gain_rms * randn(1,Nch);
    inst.C_mis     = err.cap_mismatch_sigma  * randn(cfg.Nstage1,Nch);
    inst.comp_off  = err.comp_offset_sigma   * randn(1,Nch);

    % 逐通道采样 & 量化
    ch_out = cell(1,Nch);
    for k = 1:Nch
        % 名义采样时刻 + 通道 skew + 抖动
        t_nom = ( (0:Lch-1)*Ts_ch ) + (k-1)*Ts_ch/Nch;
        t_real = t_nom + inst.tskew(k) ...
                      + err.jitter_rms * randn(size(t_nom));
        xs = sin(2*pi*cfg.Fin*t_real);
        ch_out{k} = pipeline_channel(xs, cfg, inst, k);
    end
    
    % 交织
    dig_out = re_interleave(ch_out);
    meta.inst = inst;
end


function dig = pipeline_channel(xs, cfg, inst, k)
    N = numel(xs);
    dig = zeros(1,N,'uint16');
    for n = 1:N
        
        [d1, res1] = sar_quant(xs(n), cfg, cfg.Nstage1, inst.C_mis(:,k), inst.comp_off(k));
        
        res1_amp   = res1 * 2^(cfg.Nstage1) * inst.g_stage(k) - 1;
        
        [d2, ~]    = sar_quant(res1_amp, cfg, cfg.Nstage2, [], inst.comp_off(k));
        
        dig(n) = bitshift(uint16(d1), cfg.Nstage2) + uint16(d2);
    end
    % 加 TI 通道增益/偏移失配
    dig = uint16( double(dig) .* inst.ti_gain(k) ...
                + inst.ti_offset(k) * 2^(cfg.Nbits-1) );
end


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
    code    = uint16( sum(bits .* 2.^(Nbit-1:-1:0)) );
    residue = x - acc;
end


function y = re_interleave(C)
    Nch = numel(C);   L = numel(C{1});
    y   = zeros(1, Nch*L, 'uint16');
    for k = 1:Nch
        y(k:Nch:end) = C{k};
    end
end


function [enob, sinad, sfdr, thd] = dynamic_metrics(code, Nbit, Fs, Fin)
    x  = double(code) - 2^(Nbit-1);
    N  = length(x);
    w  = win_blackmanharris(N);
    X  = fft(x .* w');
    X  = X(1:N/2);
    X2 = abs(X).^2;
    X2 = X2 / max(X2);                       % 归一化至 dBFS
    fbin = round(Fin/Fs * N);
    DC = sum(X2(1:10))
    Psig = sum(X2(fbin - 100 : fbin + 100))
                          
    % THD (2‑5 次谐波)
    thd_pow = 0;
    for h = 2:5
        idx = h*fbin;
        if idx <= length(X2) - 50; thd_pow = thd_pow + sum(X2(idx-50:idx+50)); end
    end
    thd = 10*log10(thd_pow/Psig);
    Pnoise = sum(X2) - DC - Psig -thd_pow     % 去掉 DC 与基波
    sinad = 10*log10(Psig/Pnoise);
    enob  = (sinad - 1.76)/6.02;
    sfdr  = 10*log10(Psig / max(X2(fbin:end)));
end



function w = win_blackmanharris(N)
    n  = (0:N-1)';
    a0 = 0.35875;
    a1 = 0.48829;
    a2 = 0.14128;
    a3 = 0.01168;
    w  = a0 ...
       - a1*cos(2*pi*n/(N-1)) ...
       + a2*cos(4*pi*n/(N-1)) ...
       - a3*cos(6*pi*n/(N-1));
end


function plot_fft_dbfs(code, Fs)
    x  = double(code);
    N  = length(x);
    w  = win_blackmanharris(N);
    X  = fft( (x-mean(x)) .* w');
    P  = 20*log10(abs(X(1:N/2))/max(abs(X)));
    f  = (0:N/2-1)/N * Fs / 1e6; % MHz
    plot(f, P); grid on;
    xlabel('Frequency (MHz)'); ylabel('Amplitude (dBFS)');
end
