clear; clc; close all;

%cfg
cfg.Nbits       = 12;          % 总分辨率
cfg.Nstage1     = 6;           % Stage‑1 SAR bits
cfg.Nstage2     = 6;           % Stage‑2 SAR bits
cfg.Nchan       = 4;           % 交织通道数 (≥4)
cfg.Fs_total    = 2e9;         % 总采样率 2 GS/s 静态分析无意义
cfg.Fin         = 100e6;       % 测试正弦频率 静态分析无意义
cfg.SimTime     = 5e-6;        % 仿真时长 (5 µs) 静态分析无意义
cfg.Vref        = 1;           % 基准电压 ±1 V (全差分为 2 Vpp)

%error设置
err.jitter_rms          = 0;    % 时钟抖动 σ 由于是静态分析设置为0
err.tskew_rms           = 0;    % 通道时序偏差 σ 由于是静态分析设置为0
err.ti_gain_sigma       = 0;       % ±1 %   TI 通道增益 σ
err.ti_offset_sigma     = 0;       % 1 mV     TI 通道失调 σ
err.interstage_gain_rms = 0;       % ±3 %   级间增益 σ
err.cap_mismatch_sigma  = 0;      % 0.1 %  电容 DAC 失配 σ
err.comp_offset_sigma   = 0;     % 0.5 mV 比较器失调 σ

rng(1);   % 固定随机种子便于复现

%这里生成用于得到转移曲线和直方图测试结果的斜坡信号
x = -0.99999:0.00001:1;
y = TI_PipeSAR_core(x, cfg, err);
plot(x, y);

[dnl, inl]=static_dnl_inl(y, 12);
DNL=max(abs(dnl))
INL=max(abs(inl))

%顶层
function [dig_out, meta] = TI_PipeSAR_core(x, cfg, err)


    Nch   = cfg.Nchan;
    
    
    % 生成每通道固定误差实例
    %inst.tskew     = err.tskew_rms           * randn(1,Nch);
    inst.ti_gain   = 1 + err.ti_gain_sigma   * randn(1,Nch);
    inst.ti_offset =     err.ti_offset_sigma * randn(1,Nch);
    inst.g_stage   = 1 + err.interstage_gain_rms * randn(1,Nch);
    inst.C_mis     = err.cap_mismatch_sigma  * randn(cfg.Nstage1,Nch);
    inst.comp_off  = err.comp_offset_sigma   * randn(1,Nch);

    % 逐通道采样 & 量化
    ch_out = cell(1,Nch);
    for k = 1:Nch
        
        xs = x(k:Nch:end);
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

function [dnl, inl] = static_dnl_inl(code, Nbit)
    cd   = double(code);
    nb   = 2^Nbit;
    h = histcounts(cd, 0:nb);
    ex   = mean(h);
    dnl  = (h - ex)/ex;
    inl  = cumsum(dnl);
end