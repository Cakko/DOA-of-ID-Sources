function [CRB] = func_CRB_MultiIDS(xm,theta,angularSpread,SNR,L,wavelength)
%FUNC_CRB_MULTIIDS Summary of this function goes here
%   Detailed explanation goes here
% Written by Jiaqi LI
% Date: 2024/05/07
% 本函数采用Slepian-Bangs Formula计算多非相干分布式信源的CRB
lxm = length(xm);
Q = length(theta);

%% 接收信号协方差模型
% 有如下假设
Sigma_s2 = 10^(SNR/10);
Sigma_n2 = 1;
rho = Sigma_s2 * 1; % 假设每条散射径的增益为1
eta = [theta;angularSpread];

% 生成由eta决定的期望协方差矩阵R和归一化信号协方差PHI(PHI是存储了Q个信息矩阵:MXMXQ)
[PHI,R,invR] = genR(xm,lxm,Q,eta,wavelength,rho,Sigma_n2);

%% 
Ut = zeros(lxm^2,Q);
Ua = zeros(lxm^2,Q);
Vs = zeros(lxm^2,Q);
R05 = pinv(R^0.5);
for qq = 1:Q
    eta_sel = eta(:,qq); % 选择第qq个信源计算Deta
    eta1 = eta_sel(1);  eta2 = eta_sel(2);
    PHIq = PHI(:,:,qq); 
    % 首先为了简便，给出C矩阵
    [C1,C2,C3] = genC(eta2,xm,lxm,wavelength); 
    % 这里给出PHI对eta的各偏导项, Dt, Ds
    [Dt,Ds] = genDPHI(eta1,rho,PHIq,C1,C2,C3);

    Ut(:,qq) = vec(R05*Dt*R05);
    Ua(:,qq) = vec(R05*Ds*R05);
    Vs(:,qq) = vec(R05*PHIq*R05);
end
Vn = vec(invR);
U = [Ut,Ua];
V = [Vs,Vn];
PITV = eye(lxm^2) - V*pinv(V'*V)*V'; % V的正交投影

CRB = real(diag(pinv(U'*PITV*U)))/L;
CRB = sqrt((CRB(1 : 2 * Q))) * 180 / pi;

%%%--------------------------MAIN FUNCTION END--------------------------%%%
%%%--------------------------NEXT: NESTED FUNCTION--------------------------%%%
function [rsl_vec] = vec(A)
    % vectorization
    rsl_vec = reshape(A,[],1);
end

function [PHI,R,invR] = genR(xm,lxm,Q,eta,wavelength,rho,Sigma_n2)
    % phi0, PHI是归一化的信号协方差矩阵
    % PHI存储Q个MXM维的信息
    PHI = zeros(lxm,lxm,Q);
    Rs = 0;
    for qqq = 1:Q
        psi0 = eye(lxm);
        for kk = 1:lxm
            for ll = 1:kk-1
                psi0(kk,ll) = exp(1i*2*pi/wavelength*(xm(kk)-xm(ll))*sin(eta(1,qqq)))...
                    .* exp(-0.5 * (2*pi/wavelength*(xm(kk)-xm(ll))*eta(2,qqq)*cos(eta(1,qqq))).^2); % Gaussian Dist.
                psi0(ll,kk) = conj(psi0(kk,ll));
            end
        end
        Rs = Rs + rho * psi0;
        PHI(:,:,qqq) = psi0;
    end
    % R是ECM
    R = Rs + Sigma_n2*eye(lxm);
    invR = inv(R);
end

function [C1,C2,C3] = genC(asb,xm,lxm,wavelength)
    C1 = zeros(lxm);    
    C2 = zeros(lxm);
    C3 = zeros(lxm);
    for kkk = 1:lxm
        for lll = 1:kkk-1
            C1(kkk,lll) = 2*pi/wavelength*(xm(kkk)-xm(lll));
            C2(kkk,lll) = (2*pi/wavelength*(xm(kkk)-xm(lll))*(asb))^2;
            C3(kkk,lll) = (2*pi/wavelength*(xm(kkk)-xm(lll)))^2*(asb);
            C1(lll,kkk) = -C1(kkk,lll);
            C2(lll,kkk) = C2(kkk,lll);
            C3(lll,kkk) = C3(kkk,lll);       
        end
    end
end

function [Dt,Ds] = genDPHI(sita,rho,PHI,C1,C2,C3)
    % 计算对PHI对于每组eta的导数
    Dt = rho * (1i*cos(sita).*C1+0.5*sin(2*sita).*C2) .* PHI;
    Ds = rho * (-(cos(sita))^2.*C3) .* PHI;
end


end

