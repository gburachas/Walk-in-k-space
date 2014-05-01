
function x =  findroot(f,eps)

x0 = 0.5; % initial value
x = x0;

while abs(f(x))>eps
    if f(x)>0
        x=x/2;
    elseif f(x)<0
        x = (x+1)/2;
    else
        return x;
    end
end

return x
