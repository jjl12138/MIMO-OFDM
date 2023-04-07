clear;
clc;

%% 参数设置
% 
ofdm.Nb      = 100;                 % ofdm符号数
ofdm.Nt      = 2;                   % 发射天线数
ofdm.Nr      = 2;                   % 接收天线数
ofdm.K       = 64;                   % 子载波数
ofdm.G       = 1/4;                 % 保护间隔
ofdm.Mod     = 16;                   % 调制阶数
ofdm.PSpace  = 1;                   % 导频间隔

chan.SNR_dB  = 0:2:30;                  % 信噪比
chan.L       = 6;                   % 信道抽头数


ofdm.PPos    = 1:(ofdm.PSpace+1):ofdm.K;    % OFDM 导频位置
ofdm.PL      = length(ofdm.PPos);           % 导频个数
ofdm.DPos    = setxor(1:ofdm.K,ofdm.PPos);  % 数据位置
ofdm.DL      = length(ofdm.DPos);           % 数据个数
ofdm.BER     = zeros(length(chan.SNR_dB),ofdm.Nb); % 误码率
chan.NMSE    = zeros(length(chan.SNR_dB),ofdm.Nb); % NMSE                             
chan.sigma   = sqrt(10.^(-0.1.*chan.SNR_dB)); % 噪声功率

% QAM调制的归一化因子
temp         = 0:ofdm.Mod-1;           
temp         = qammod(temp,ofdm.Mod);  
temp         = abs(temp).^2;           
temp         = mean(temp);             
ofdm.ModNorm = 1/sqrt(temp);           

%% 仿真
for nb = 1:ofdm.Nb
    for ind = 1:length(chan.SNR_dB)
        snri = chan.SNR_dB(ind);

        %% 生成数据
        ofdm.d      = randi(ofdm.Mod,ofdm.DL,ofdm.Nt)-1;
    
        ofdm.dMod   = zeros(ofdm.K,ofdm.Nt);    
       
        ofdm.dMod(ofdm.DPos,:) = ofdm.ModNorm*qammod(ofdm.d,ofdm.Mod);
    
        %% 导频插入
        
        for nt = 1:ofdm.Nt
            ofdm.dMod(ofdm.PPos,nt) = exp(-sqrt(-1)*2*pi*(nt-1)*chan.L*(1:ofdm.PL).'/ofdm.PL);
        end
     
        ofdm.pow = var(ofdm.dMod(:))+abs(mean(ofdm.dMod(:)))^2;
    
        %% IFFT
       
        ofdm.ifft = sqrt(ofdm.K)*ifft(ofdm.dMod,ofdm.K);
    
        %% 循环前缀
    
        ofdm.ifftG = [ofdm.ifft(ofdm.K*(1-ofdm.G)+1:ofdm.K,:);ofdm.ifft];
        %% 信道
    
        chan.Coeff = 1/sqrt(2)*1/sqrt(chan.L)*(randn(ofdm.Nt,ofdm.Nr,chan.L)+sqrt(-1)*randn(ofdm.Nt,ofdm.Nr,chan.L));
        
        %% 信号通过信道并叠加噪声
    
        ofdm.Y = zeros(ofdm.K*(1+ofdm.G),ofdm.Nr);
    
        for nt=1:ofdm.Nt
            for nr=1:ofdm.Nr
                ofdm.Y(:,nr) = ofdm.Y(:,nr) + filter(squeeze(chan.Coeff(nt,nr,:)),1,ofdm.ifftG(:,nt));
            end
        end
    
        % adding noise
        ofdm.Y = ofdm.Y + chan.sigma(ind)*1/sqrt(2)*(         randn(ofdm.K*(1+ofdm.G),ofdm.Nr)+...
            sqrt(-1)*randn(ofdm.K*(1+ofdm.G),ofdm.Nr)     );
    
        %% 去除循环前缀

        ofdm.fftG = ofdm.Y(ofdm.K*ofdm.G+1:ofdm.K*(1+ofdm.G),:);
    
        %% FFT 
        
        ofdm.fft  = 1/sqrt(ofdm.K)*fft(ofdm.fftG,ofdm.K);
    
        %% 信道估计
    
        F = dftmtx(ofdm.K);
        F = F(:,1:chan.L);
        
        chan.CoeffEst = zeros(ofdm.Nt,ofdm.Nr,chan.L);

        for nr = 1 : ofdm.Nr
            chan.A = zeros(ofdm.PL,chan.L*ofdm.Nt);
            for nt = 1 : ofdm.Nt
                chan.A(:,(1:chan.L)+(nt-1)*chan.L) = diag(ofdm.dMod(ofdm.PPos,nt))*F(ofdm.PPos,:);
            end
%             ChanEst = pinv(chan.A)*ofdm.fft(ofdm.PPos,nr);
            ChanEst = inv(chan.A' * chan.A + chan.sigma(ind)*eye(chan.L*ofdm.Nt))*chan.A'*ofdm.fft(ofdm.PPos,nr);
            for nt = 1 : ofdm.Nt
                chan.CoeffEst(nt,nr,:) = ChanEst((1:chan.L)+(nt-1)*chan.L);
            end
        end
    
        %% MSE (Mean Square error calculation)
        chan.NMSE(ind,nb) = (norm(chan.Coeff - chan.CoeffEst,'fro')/norm(chan.Coeff,'fro'))^2;
    
    %     disp(['MSE of channel estimation (theory)     is : ',num2str(chan.MSE_Theory)])
    %     disp(['MSE of channel estimation (simulation) is : ',num2str(chan.MSE_Simulation)])
    
        %% 解调

        chan.CoeffEstFreq = fft(chan.CoeffEst,ofdm.K,3);
        chan.CoeffEstFreq = permute(chan.CoeffEstFreq,[3 1 2]);

        ofdm.H = chan.CoeffEstFreq;
  
        ofdm.dDemod = zeros(ofdm.DL,ofdm.Nt);

        for dl = 1 : ofdm.DL
            H_temp = reshape(chan.CoeffEstFreq(ofdm.DPos(dl),:,:),ofdm.Nt,[]).';
            W_MMSE = pinv(H_temp' * H_temp + chan.sigma(ind)^2*eye(ofdm.Nt))*H_temp';
            ofdm.dDemod(dl,:) = W_MMSE * squeeze(ofdm.fft(ofdm.DPos(dl),:)).';            
        end

        % demodulation
        
        ofdm.dEst = qamdemod(1/ofdm.ModNorm * ofdm.dDemod,ofdm.Mod);
        % BER calculation
        [~,ofdm.BER(ind,nb)]  = biterr(ofdm.d(:),ofdm.dEst(:),log2(ofdm.Mod));

    end
end


% 2. MMSE 均衡器
% 1. 多天线的信号检测，把接收信号表示成矩阵-向量形式
% disp(['BER is = ',num2str(BER)]);
% disp(['MSE of channel estimation (simulation) is : ',num2str(chan_MSE_Simulation)]);

mean_NMSE = mean(chan.NMSE,2);
mean_BER = mean(ofdm.BER,2);

subplot(1,2,1);
semilogy(chan.SNR_dB,mean_NMSE,'rx-');
xlabel('SNR [dB]');
ylabel(['NMSE']);


subplot(1,2,2);
semilogy(chan.SNR_dB,mean_BER,'cd-');
xlabel('SNR [dB]');
ylabel(['BER']);






