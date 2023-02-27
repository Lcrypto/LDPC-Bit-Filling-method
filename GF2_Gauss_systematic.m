function [outputH,outputG]=GF2_Gauss_systematic(inputx,inputy);
%%%: Input check matrix x, perform Gaussian elimination on binary domain to generate H,
% then generate generation matrix G corresponding to H, output H, G ----% %%
%%%%%----Reference basis: x (Gaussian elimination) -> H=[I|P] -> G=[P'|I]G ------%%%%%
H = inputx;
flag = inputy;%flag, to determine whether the last generated G is the unit matrix on the right (flag = 1), or the unit matrix on the left (flag = else)
H1 = H;
[m,n] = size(H);

%%%---- Check each line of H one by one, H(i,:), where i=1:m, if H(i,i)==0, then ----%%%%
%%%----Find the first non-zero element in the row, record the column of the non-zero element, and compare the column of H with ----%%%%
%%%----The i-th column of H is exchanged; if H(i,i)==1, there is no need for column exchange. Then check after the exchange ----%%%%
%%%----or H(:,i) that does not need to be exchanged, from line i+1 to line m below the H(i,i) element ----%%%%
%%%----Whether there are non-zero elements, if yes, add the i-th line of H to this line, if not, then check the value of H ----%%%%
%%%----The next line, cycle. The result is to convert H into a bottom left with all 0, H(i,i) all 1 ----%%%%
%%%----matrix, recorded as Hmid, and then superimpose Hmid from row m to row j, j=m-1:-1:1, -----%%%%
%%%---- After the loop, the left half is the matrix [I|P] -----%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%

for i=1:m; % Check each row of H one by one, H(i,:)
   if H1(i,i)==0;% If H(i,i)==0 branch
      j=min(find(H1(i,:)));%  If H(i,i)==0 branch
      H1(:,[i j])=H1(:,[j,i]);% Exchange the column of H with the i-th column of H, the H matrix is used as the matrix to retain the output, and H1 is used for operation
      H(:,[i j])=H(:,[j,i]);
      for k=i+1:m;% Whether there are any non-zero elements from the i+1th line to the m line below the H(i,i) element
          if H1(k,i)==1;% If yes, superimpose the i-th row of H to this row
             H1(k,:)=H1(i,:)+H1(k,:);
             H1(k,:)=mod(H1(k,:),2);
          end;% None, check the next line of H
      end;
   else;% If the branch of H(i,i)==1, there is no need to exchange lines, only need to check whether there are non-zero elements in the i+1th line to the m line below the H(i,i) element
      for k=i+1:m;
          if H1(k,i)==1;% If yes, superimpose the i-th row of H to this row
             H1(k,:)=H1(i,:)+H1(k,:);
             H1(k,:)=mod(H1(k,:),2);
          end;
      end;
   end;% None, check the next line of H
end;

for i=m:-1:2;% Superimpose from line m to line j, j=m-1:-1:1
  for k=i-1:-1:1;
    if H1(k,i)==1;
       H1(k,:)=H1(i,:)+H1(k,:);
       H1(k,:)=mod(H1(k,:),2);
    end;
  end;% After the cycle, the matrix H1 = [I|P] with the left half as the unit matrix and the H corresponding to H1 are obtained
end;

%%%%-----------------Step 3: Get matrix G----------------%%%%%
if flag ==1;% The right side of % G is the identity matrix
   PP = H1(:,m+1:n);% Take out the m+1 to nth columns of H1
   G = [PP.' diag(ones(1,n-m))];%  H1=[I|P] -> G=[P'|I]
   outputH = H;
   outputG = G;
else;% The left side of G is the identity matrix
   PP = H1(:,m+1:n);
   GT = [diag(ones(1,n-m));PP];
   outputH = [H(:,m+1:n) H(:,1:m)];
   outputG = GT.';
end;