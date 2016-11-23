

clc; clear variables; close all;

n1=6; n2= 6; n3=6; n4=6; n5=6; % Total no of signals in data sets
m1=14;m2=18;m3=20; m4=22; m5=24;% Total no of sensors in data sets

var_1 = [1 1 1 1 1 1];  % Variance of sources in data sets
var_2 = [1 1 1 1 1 1];
var_3 = [1 1 1 1 1 1];
var_4 = [1 1 1 1 1 1];
var_5 = [1 1 1 1 1 1];

M = 400;

% Cross Correlation matrices satisfying the assumption in [6],[12] and [13]
% R12 = diag([0.8,  0.8, 0,   0    0   0 ].*sqrt(var_1/2).*sqrt(var_2/2));
% R23 = diag([0.9,  0.9, 0  , 0    0   0 ].*sqrt(var_2/2).*sqrt(var_3/2));
% R31 = diag([0.8,  0.8, 0,   0    0   0 ].*sqrt(var_3/2).*sqrt(var_1/2));
% R14 = diag([0.85, 0.8, 0,   0,   0   0 ].*sqrt(var_1/2).*sqrt(var_4/2));
% R24 = diag([0.9,  0.8, 0,   0    0   0 ].*sqrt(var_2/2).*sqrt(var_4/2));
% R34 = diag([0.9,  0.8, 0,   0    0   0 ].*sqrt(var_3/2).*sqrt(var_4/2));
% R15 = diag([0.85, 0.8, 0,   0    0   0 ].*sqrt(var_1/2).*sqrt(var_5/2));
% R25 = diag([0.9,  0.75, 0,   0    0   0 ].*sqrt(var_2/2).*sqrt(var_5/2));
% R35 = diag([0.9,  0.8, 0,   0    0   0 ].*sqrt(var_3/2).*sqrt(var_5/2));
% R45 = diag([0.85, 0.8, 0,   0    0   0 ].*sqrt(var_4/2).*sqrt(var_5/2));

% Cross Correlation matrices not satisfying the assumption in [6],[12] and [13]
R12 = diag([0.8,  0.8, 0.7,   0    0.7   0 ].*sqrt(var_1/2).*sqrt(var_2/2));
R23 = diag([0.9,  0.9, 0.5  , 0    0   0.6 ].*sqrt(var_2/2).*sqrt(var_3/2));
R31 = diag([0.8,  0.8, 0.7,   0    0   0.5 ].*sqrt(var_3/2).*sqrt(var_1/2));
R14 = diag([0.85, 0.8, 0.7,   0.7,   0   0 ].*sqrt(var_1/2).*sqrt(var_4/2));
R24 = diag([0.9,  0.8, 0.6,   0.6    0   0.6 ].*sqrt(var_2/2).*sqrt(var_4/2));
R34 = diag([0.9,  0.8, 0,   0    0.5   0.6 ].*sqrt(var_3/2).*sqrt(var_4/2));
R15 = diag([0.85, 0.8, 0,   0.7    0.5   0 ].*sqrt(var_1/2).*sqrt(var_5/2));
R25 = diag([0.9,  0.75, 0,   0    0.5   0.6 ].*sqrt(var_2/2).*sqrt(var_5/2));
R35 = diag([0.9,  0.8, 0,   0.5    0   0.5 ].*sqrt(var_3/2).*sqrt(var_5/2));
R45 = diag([0.85, 0.8, 0,   0.6    0.6   0.6 ].*sqrt(var_4/2).*sqrt(var_5/2));

d = 2;

B = 500; % Number of Bootstrap iterations
alpha = 0.05;

num_iterations = 5*1e2;
SNR_Vector = [ -10 -7.5 -5 -4 -3 -2.5 -1 0 5 10 15 ];

m=[m1,m2,m3,m4,m5];

for sample = 1:length(SNR_Vector)
    
    var_n = var_1(1)/(10^(SNR_Vector(sample)/10))
    var_n_1 = var_n; var_n_2 = var_n; var_n_3 = var_n; var_n_4 = var_n; var_n_5 = var_n; % Variance of noise
    
    for iter=1:num_iterations
        
        display(['iteration=' num2str(iter)]);
        
        % Source Generation
        % Generating Source Matrices according to the variance and correlation coefficient specified above
        
        idx1= n1; idx2=n2+n1; idx3=n1+n2+n3; idx4=n1+n2+n3+n4; idx5=n1+n2+n3+n4+n5;
        
        E = zeros(idx5);
        
        E(1:idx1,1:idx1) = diag(var_1/2);
        E(idx1+1:idx2,idx1+1:idx2) = diag(var_2/2);
        E(idx2+1:idx3,idx2+1:idx3) = diag(var_3/2);
        E(idx3+1:idx4,idx3+1:idx4) = diag(var_4/2);
        E(idx4+1:idx5,idx4+1:idx5) = diag(var_5/2);
        
        E(1:idx1,idx1+1:idx2) = R12; E(idx1+1:idx2,1:idx1) = R12';
        E(idx1+1:idx2,idx2+1:idx3) = R23; E(idx2+1:idx3,idx1+1:idx2) = R23';
        E(idx2+1:idx3,1:idx1) = R31; E(1:idx1,idx2+1:idx3) = R31';
        E(1:idx1,idx3+1:idx4) = R14; E(idx1+1:idx2,idx3+1:idx4) = R24; E(idx2+1:idx3,idx3+1:idx4) = R34;
        E(idx3+1:idx4,1:idx1) = R14'; E(idx3+1:idx4,idx1+1:idx2) = R24'; E(idx3+1:idx4,idx2+1:idx3) = R34';
        E(1:idx1,idx4+1:idx5) = R15; E(idx1+1:idx2,idx4+1:idx5) = R25; E(idx2+1:idx3,idx4+1:idx5) = R35; E(idx3+1:idx4,idx4+1:idx5) = R45;
        E(idx4+1:idx5,1:idx1) = R15'; E(idx4+1:idx5,idx1+1:idx2) = R25'; E(idx4+1:idx5,idx2+1:idx3) = R35'; E(idx4+1:idx5,idx3+1:idx4) = R45';
        
        mu = zeros(idx5,1);
        F = mvnrnd(mu,E,M);
        
        S1_r = F(:,1:idx1)'; S2_r = F(:,idx1+1:idx2)';  S3_r= F(:,idx2+1:idx3)'; S4_r= F(:,idx3+1:idx4)'; S5_r= F(:,idx4+1:idx5)';
        
        F = mvnrnd(mu,E,M);
        S1_c = F(:,1:idx1)'; S2_c = F(:,idx1+1:idx2)';  S3_c= F(:,idx2+1:idx3)'; S4_c= F(:,idx3+1:idx4)'; S5_c= F(:,idx4+1:idx5)';
        
        S1 = (S1_r + 1i*S1_c); S2 = (S2_r + 1i*S2_c); S3 = (S3_r + 1i*S3_c); S4 = (S4_r + 1i*S4_c); S5 = (S5_r + 1i*S5_c);
        
        % Mixing Matrices Generation
        
        A1 = orth(rand(m1,n1)+ 1i*rand(m1,n1));
        A2 = orth(rand(m2,n2)+ 1i*rand(m2,n2));
        A3 = orth(rand(m3,n3)+ 1i*rand(m3,n3));
        A4 = orth(rand(m4,n4)+ 1i*rand(m4,n4));
        A5 = orth(rand(m5,n5)+ 1i*rand(m5,n5));
        
        % Noise Matrices Generation
        
        N1 = sqrt(var_n_1/2)*(randn(m1,M)+ 1i*randn(m1,M));
        N2 = sqrt(var_n_2/2)*(randn(m2,M)+ 1i*randn(m2,M));
        N3 = sqrt(var_n_3/2)*(randn(m3,M)+ 1i*randn(m3,M));
        N4 = sqrt(var_n_4/2)*(randn(m4,M)+ 1i*randn(m4,M));
        N5 = sqrt(var_n_5/2)*(randn(m5,M)+ 1i*randn(m5,M));
        
        %  Colored Noise
        A= [sqrt(8) sqrt(1) sqrt(1) ]; % Filter coefficient vector
        V_N1 = filter(1,A,N1);
        V_N2 = filter(1,A,N2);
        V_N3 =  filter(1,A,N3);
        V_N4 =  filter(1,A,N4);
        V_N5 =  filter(1,A,N5);
        
        % Data sets generation
        
        X1 = A1*S1 + V_N1;
        X2 = A2*S2 + V_N2;
        X3 = A3*S3 + V_N3;
        X4 = A4*S4 + V_N4;
        X5 = A5*S5 + V_N5;
        
        % CCA and ESTIMATION OF CORRELATED SIGNALS
        
        [~,~,Vx1] = svd(X1,'econ');
        [~,~,Vx2] = svd(X2,'econ');
        [~,~,Vx3] = svd(X3,'econ');
        [~,~,Vx4] = svd(X4,'econ');
        [~,~,Vx5] = svd(X5,'econ');
        
        K = svd(Vx1'*Vx2*Vx2'*Vx3*Vx3'*Vx4*Vx4'*Vx5*Vx5'*Vx3*Vx3'*Vx1*Vx1'*Vx4*Vx4'*Vx2*Vx2'*Vx5*Vx5'*Vx1);
        
        % Wu method
        
        [~,~, Vx2x3x4x5] = svd([X2;X3;X4;X5],'econ');
        [~,~, Vx3x4x5] = svd([X3;X4;X5],'econ');
        [~,~, Vx4x5] = svd([X4;X5],'econ');
        
        % Wu's Likelihood - Order 1,2,3,4,5
        Kx1_x2x3x4x5 = svd(Vx1'*Vx2x3x4x5);
        Kx2_x3x4x5 = svd(Vx2'*Vx3x4x5);
        Kx3_x4x5 = svd(Vx3'*Vx4x5);
        Kx4_x5 = svd(Vx4'*Vx5);
        
        p_arrays=5; m = [m1,m2,m3,m4,m5];
        K_Wu = zeros(p_arrays-1,max(m));
        K_Wu(1,1:m1) = Kx1_x2x3x4x5'; K_Wu(2,1:m2) = Kx2_x3x4x5'; K_Wu(3,1:m3) = Kx3_x4x5'; K_Wu(4,1:m4) = Kx4_x5';
        C_S_Wu(iter) = Wu_Correlated_Signals_Unequal_Sensors(M,m,p_arrays,K_Wu,a,beta);
        
        % Bootstrap Operation
        for b=1:B
            
            [Y1,I] = datasample(X1,M,2); % Bootstrap data matrix
            Y2 = X2(:,I);
            Y3 = X3(:,I);
            Y4 = X4(:,I);
            Y5 = X5(:,I);
            
            
            % Centering around the mean
            Y1 = Y1 - repmat(mean(Y1,2),1,M);
            Y2 = Y2 - repmat(mean(Y2,2),1,M);
            Y3 = Y3 - repmat(mean(Y3,2),1,M);
            Y4 = Y4 - repmat(mean(Y4,2),1,M);
            Y5 = Y5 - repmat(mean(Y5,2),1,M);
            
            [~,~,Vy1] = svd(Y1,'econ');
            [~,~,Vy2] = svd(Y2,'econ');
            [~,~,Vy3] = svd(Y3,'econ');
            [~,~,Vy4] = svd(Y4,'econ');
            [~,~,Vy5] = svd(Y5,'econ');
            K_star = svd(Vy1'*Vy2*Vy2'*Vy3*Vy3'*Vy4*Vy4'*Vy5*Vy5'*Vy3*Vy3'*Vy1*Vy1'*Vy4*Vy4'*Vy2*Vy2'*Vy5*Vy5'*Vy1);
            
            K_star_matrix(:,b) = K_star;
            
        end
        
        % Hypothesis Testing using Bootstrap
 
        k_cap_type_3_wo_bc = hypothesis_testing_bootstrap(3,alpha,K,K_star_matrix,B);
        k_cap_vector_3_wo_bc(iter) = k_cap_type_3_wo_bc
        
        % Max-Min Multiple data sets [12]
        
        [~,~, Vx2x3x4x5] = svd([X2;X3;X4;X5],'econ');
        [~,~, Vx3x4x5] = svd([X3;X4;X5],'econ');
        [~,~, Vx4x5] = svd([X4;X5],'econ');
        P_fa = 0.001;
        rmax = min([floor(2*M/(3*5)),m1]);
        [No_CS_Max_Min_5, r_d] = Max_Min_Five_Datasets_Complex_Fast(M,P_fa, Vx1,Vx2x3x4x5,Vx2, Vx3x4x5, Vx3, Vx4x5, Vx4, Vx5, rmax);
        No_CS_Vector_Max_Min_5(iter) = No_CS_Max_Min_5
        
        % Detector [13]
        intuitive_rr_eusipco(iter) = intuitive_rank_red_eusipco(X1,X2,X3,X4,X5,M,P_fa)
        
        warning('off','all')
        
    end
    
    Prob_detection_bt_3_wo_bc(sample) = length(find(k_cap_vector_3_wo_bc == d))/num_iterations;
    Prob_detection_Wu(sample) = length(find(C_S_Wu == d))/num_iterations;  
    Prob_detection_Max_Min_5(sample) = length(find(No_CS_Vector_Max_Min_5 == d))/num_iterations;
    Prob_detection_Intuitive_rr_eusipco_5(sample) = length(find(intuitive_rr_eusipco == d))/num_iterations;
     
    Mean_bt_3_wo_bc(sample) = mean(k_cap_vector_3_wo_bc)
    Var_bt_3_wo_bc(sample) = var(k_cap_vector_3_wo_bc);
    
    Mean_Wu(sample) = mean(C_S_Wu)
    Var_Wu(sample) = var(C_S_Wu);
    Mean_Max_Min_5(sample) = mean(No_CS_Vector_Max_Min_5)
    Var_Max_Min_5(sample) = var(No_CS_Vector_Max_Min_5);
    Mean_Intuitive_rr_eusipco_5(sample) = mean(intuitive_rr_eusipco)
    Var_Intuitive_rr_eusipco_5(sample) = var(intuitive_rr_eusipco);
    
end


figure();
h3 = plot(SNR_Vector,Mean_bt_3_wo_bc,'ro--','markersize',12,'Linewidth',2); hold on;
h4 = plot(SNR_Vector,Mean_Wu,'ks-.','markersize',12,'Linewidth',2);
h5 = plot(SNR_Vector,Mean_Max_Min_5,'b*--','markersize',12,'Linewidth',2);
h6 = plot(SNR_Vector,Mean_Intuitive_rr_eusipco_5,'mx-.','markersize',12,'Linewidth',2);
h7 = plot(SNR_Vector,2*ones(size(SNR_Vector)),'g-','Linewidth',2);
h = [h3,h4,h5,h6,h7];
a = xlabel('SNR (dB)','fontsize',20,'FontName','Times New Roman');
b = ylabel('Mean value of $\hat{d}$','fontsize',20,'FontName','Times New Roman');
c = legend('Proposed detector','ITC MDL [6]','Max-Min detector [12]','Detector based on [13]','True number','FontName','Times New Roman','fontsize',20);
set(a,'interpreter','latex');
set(b,'interpreter','latex');
set(c,'interpreter','latex');
% set(gca,'FontSize',20);
set(gca,'FontName','Times New Roman','FontSize',20);

