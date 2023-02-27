% This program is to construct a random LDPC codes by using the bit filling method[1]. 
 % 
 % Ref:
 % [1] J. Campellot, D. S. Modhat and S. Rajagopalant. "Designing LDPC
 % Codes Using Bit-Filling". ICC 2001, pp. 55–59, 2001.
 %
 % The paper[1] can be downloaded from Web site of IEEE Explore.
 %
 
 %   Copyright (C) 2007, Guangrong Fan. MCL. BIT.
 %   $Revision: 1.0 $  $Date: 2007/07/30 21:12:41 $
 
 clear;
 clc;
 
 %============ Parameters related to parity check matrix ============%
 M = 100;        % The number of check matrix rows
 N =200;         % Number of parity-check matrix columns VN
 cols_w = 4;      % The column weight
 girth =6;       % The demand girth length
 right_eye = 1;       % (flag = 1), or the unit matrix on the left (flag = else)
 %============ Initialization ============%
 H = sparse(M,N);
 n = 0;
 U1 = []; 
 ck_deg(1:1:M) = 0;
 
 max_rows = ceil(cols_w*N/M);    % The maximum row weight

 %======================= Bit Filling =======================%
 while ( (n == 0) | ((i==cols_w) & (~isempty(F0))) ) 
    A = 1:1:M;
   
    %For c belongs to U1, set ck_deg(c) to increase by 1, let H(c,n) = 1
    for j = 1:1:length(U1)
        H(U1(j),n) = 1;
        ck_deg(U1(j)) = ck_deg(U1(j)) + 1;
    end
    l = find(ck_deg == max_rows); 
    A(l) = [];     %%Remove the node whose degree of check node is greater than max_rows                                                           
    %%%%%%%%%%%%%%%%%%%%
    i = 0;
    U1 = [];
    U = [];
    F0 = A;
    while ( (i < cols_w) & (~isempty(F0)) )
            for j = 1:1:length(U);
                l = (find(F0==U(j)));
                F0(l) = [];   %Calculate F0 = A\U
            end
            if ( ~isempty(F0) )  
                [min_y,min_i] = min_v(ck_deg(F0));
                min_v_num = length(min_i);
                if ( min_v_num > 1)                  
                  new_cn = F0(min_i(unidrnd(min_v_num)));
                else  
                  new_cn = F0(min_i);
                end;
                U1 = [U1,new_cn];  
                U_tmp = srh_cn_set(new_cn,H,girth);
                U = [U1,U_tmp];
                U = unique(U);
                i = i + 1;                             
            end
         
    end  %%while
    if ( (i == cols_w) & (~isempty(F0)))
      n = n+1; 
    end;
 end %%while
 %%%%%%%%%%%%%%%%%%%%%%%%%

 %================ Get the parity check matrix H ================%
 if ( n > N )
   H = sparse(H(:,1:1:N));
 else
   H = sparse(H(:,1:1:n));    
 end;
 
H=full(H);
x=H;
[Hm,Hn]=size(x);
x=x(randperm(Hm),:);
x=x(:,randperm(Hn));
[outputH,outputG]=GF2_Gauss_systematic(x,right_eye);%  
 
 save('H.mat','outputH');
 
 save('G.mat','outputG');
