% check differentiation...

format long

x2 = [0.50, 0.60, 0.83];
x3 = [0.65, 0.72, 0.21];

x1 = [0.44, 0.11,-0.91];

dx = 1e-3*[0.05, 0.023, -0.07];

U0 = dot(x3-x2,x2-x1) / norm(x3-x2) / norm(x2-x1)

gradU = -(x3-x2) / norm(x3-x2) / norm(x2-x1) - dot(x3-x2,x2-x1) / norm(x3-x2) * (x1-x2) / norm(x2-x1)^3;

x1 = x1+dx;
U1 = dot(x3-x2,x2-x1) / norm(x3-x2) / norm(x2-x1)

U0 + dot(gradU,dx)