v=[5+6i; 7+9i; 3+2i; 8+4i];
[u,W]=approx(v);
disp(v-W*u)
function [u,W]=approx(v)
    N=length(v);
    u=zeros(2,1);
    u(1)=(max(abs(v))+min(abs(v)))/2;
    u(2)=(max(abs(v))-min(abs(v)))/2;
    W=zeros(N,2);
    for i=1:N
        a1 = abs(v(i))^2 + max(abs(v)) * min(abs(v));
        b1 = abs(v(i)) * (max(abs(v)) + min(abs(v)));
        theta_i = angle(v(i)) - acos(min(max(a1 / b1, -1), 1));  % Clamped inside [-1,1] to avoid domain error
        W(i,1) = exp(1j * theta_i);

        a2 = abs(v(i))^2 - max(abs(v)) * min(abs(v));
        b2 = abs(v(i)) * (max(abs(v)) - min(abs(v)));
        theta_i = angle(v(i)) + acos(min(max(a2 / b2, -1), 1));  % Clamped inside [-1,1] to avoid domain error
        W(i,2) = exp(1j * theta_i);
    end 
end