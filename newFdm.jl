using LinearAlgebra;
using Plots;
# Let's create our range
N = 100;
Δx = 1/N;
x = 0:Δx:1-Δx;
T = zeros(N,N);
for i in 1:N
    for j in 1:N
        if i == j
            T[i, j] = -2
        elseif abs(i - j) == 1
            T[i, j] = 1
        end
    end
end
H = T*(-1/(2*Δx^2));
val, vec = eigen(H);

# We must normalize our eigenvectors
y1 = vec[:,1] ;
A_1 = 1/sqrt(sum(abs2,y2));
y2 = vec[:,2] ;
plot(x,[y1,y2]);
savefig("myplot.png") ;