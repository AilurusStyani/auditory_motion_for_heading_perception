framei = 0:0.01:100;
frequency = 1/5;
y = square(frequency*framei, 20) /2 + 0.5;
z = square(frequency*framei + pi/2, 20) /2 + 0.5;
h = square(frequency*framei + pi, 20) /2 + 0.5;
q = square(frequency*framei + pi*3/2, 20) /2 + 0.5;
plot(framei,y)
hold on
plot(framei,z)
hold on
plot(framei,h)
hold on
plot(framei,q)

framei = 0:0.01:100;
frequency = 1/5;
for i = 1:8
    y = square(frequency*framei + 2*pi/8*(i-1), 10) /2 + 0.5;
    plot(framei,y)
hold on
end

% rng('Shuffle')
i = 1:8;
Shuffle(i)
