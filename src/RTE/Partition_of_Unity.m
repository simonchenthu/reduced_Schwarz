function POU = Partition_of_Unity(t,s,dx,Nv)

[v0,~] = legendre_quad(Nv-1); v0 = v0';
x = t(1):dx:s(end); 
[xx,~] = meshgrid(x,v0); 
Mx = length(t);

L = @(x,alpha) exp(-1./(1-min(x/alpha,1))); %  Prototype function

POU = zeros([size(xx),Mx]);

for j = 1:Mx
    
    xx_dis = abs(xx-(s(j)+t(j))/2);
    
    alpha_x = (s(j)-t(j))/2;
    
    POU(:,:,j) = L(xx_dis,alpha_x);
    
    %         figure(99)
    %         mesh(xx,yy,POU(:,:,j,k));
end

weight = sum(POU,3);
POU = POU./weight;
POU(isnan(POU)|isinf(POU)) = 0;

POU(:,1,1) = 1; POU(:,end,end) = 1;

if Nv == 1
    POU = reshape(POU,[],Mx);
end

% for j = 1:Mx
%     figure(99)
%     if Nv >1
%         mesh(xx,vv,POU(:,:,j));
%     else
%         plot(x,POU(:,:,j));
%     end
%     pause;
% end


end