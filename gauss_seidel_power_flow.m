%% Gauss-Seidel Power Flow Calculation
clc; clear;

%% Input Data
max_iter = 100;       % Maximum iterations
error_iter = 1e-8;    % Iteration error
baseMVA = 100;        % Base MVA
bus_data = [           % Bus data
    % Bus No.  Type (1=Slack, 2=PV, 3=PQ)  P Load (MW)  Q Load (MVAR)  Initial Voltage  Ground Admittance
    1       1       0       0       1.01    0;
    2       3       0.6     0.35    1       0;
    3       3       0.7     0.42    1       18;
    4       3       0.8     0.5     1       15;
    5       2       0.65    0.36    1       0
];
bus_data(:, 6) = bus_data(:, 6) / baseMVA; % Ground admittance normalization
gen_data = [           % Generator data
    1       0       0;
    5       1.9     0;
];
line_data = [           % Line data
    % From Bus  To Bus  R (pu)   X (pu)   Charging (pu)  Tap Ratio
    1         2       0.0108  0.0649   0              1;
    1         4       0.0235  0.0941   0              1;
    2         5       0.0118  0.0471   0              1;
    3         5       0.0147  0.0588   0              1;
    4         5       0.0118  0.0529   0              1;
    2         3       0       0.04     0              0.975
];
capacitor_data = [           % Capacitor data
    % Bus  Rating (MVAR)
    3     18;
    4     15
];

%% Calculate Admittance Matrix
n_bus = size(bus_data, 1);
n_gen = size(gen_data, 1);
n_line = size(line_data, 1);
n_capacitor = size(capacitor_data, 1);

PQ = find(bus_data(:, 2) == 3);
PV = find(bus_data(:, 2) == 2);
V0 = find(bus_data(:, 2) == 1);

G_B = zeros(n_bus, n_bus); % Admittance matrix
for i = 1:n_line
    from_bus = line_data(i, 1);
    to_bus = line_data(i, 2);
    Y = 1 / complex(line_data(i, 3), line_data(i, 4));
    G_B(from_bus, to_bus) = G_B(from_bus, to_bus) - Y / line_data(i, 6);
    G_B(to_bus, from_bus) = G_B(to_bus, from_bus) - Y / line_data(i, 6);
    G_B(from_bus, from_bus) = G_B(from_bus, from_bus) + Y / (line_data(i, 6)^2) + complex(0, line_data(i, 5));
    G_B(to_bus, to_bus) = G_B(to_bus, to_bus) + Y + complex(0, line_data(i, 5));
end
for i = 1:n_bus
    G_B(i, i) = G_B(i, i) + 1j * bus_data(i, 6);
end

%% Initialize Variables
P = zeros(n_bus, 1); Q = zeros(n_bus, 1); V = bus_data(:, 5); theta = zeros(n_bus, 1);
V_k = V .* exp(1j * theta);
Q_PV_k = zeros(n_bus, 1);
V_k_1 = V_k;
iteration = 0;

%% Gauss-Seidel Iteration
while iteration <= max_iter
    for i = 1:n_bus
        a = find(i == PV);
        b = find(i == gen_data(:, 1));
        if i ~= V0
            S = -(bus_data(i, 3) + 1j * bus_data(i, 4));
            if ~isempty(a)
                Q_PV_k(i) = imag(V_k(i) * (G_B(i, :) * V_k_1)');
                S = gen_data(b, 2) - bus_data(i, 3) + 1j * Q_PV_k(i);
            end
            V_k_1(i) = ((S / V_k(i))' - G_B(i, :) * V_k_1 + G_B(i, i) * V_k(i)) / G_B(i, i);
        end
        if ~isempty(a)
            V_k_1(i) = bus_data(i, 5) * V_k_1(i) / abs(V_k_1(i));
        end
    end
    c1 = max(abs(abs(V_k_1) - abs(V_k)));
    c2 = max(abs(angle(V_k_1) - angle(V_k)));
    if c1 <= error_iter && c2 <= error_iter
        break;
    end
    V_k = V_k_1;
    iteration = iteration + 1;
end

%% Line Flows
V = abs(V_k);
theta = angle(V_k);
P_b = zeros(n_line, 2);
Q_b = zeros(n_line, 2);
for i = 1:n_line
    from_bus = line_data(i, 1);
    to_bus = line_data(i, 2);
    Y = -1 / complex(line_data(i, 3), line_data(i, 4));
    Y_f = Y / (line_data(i, 6)^2);
    Y_ft = Y / line_data(i, 6);
    P_b(i, 1) = -1 * V(from_bus) * V(from_bus) * real(Y_f) + V(from_bus) * V(to_bus) * (real(Y_ft) * cos(theta(from_bus) - theta(to_bus))...
        + imag(Y_ft) * sin(theta(from_bus) - theta(to_bus)));
    P_b(i, 2) = -1 * V(to_bus) * V(to_bus) * real(Y) + V(to_bus) * V(from_bus) * (real(Y_ft) * cos(theta(to_bus) - theta(from_bus))...
        + imag(Y_ft) * sin(theta(to_bus) - theta(from_bus)));
    Q_b(i, 1) = V(from_bus) * V(from_bus) * imag(Y_f) + V(from_bus) * V(to_bus) * (real(Y_ft) * sin(theta(from_bus) - theta(to_bus))...
        - imag(Y_ft) * cos(theta(from_bus) - theta(to_bus))) - V(from_bus) * V(from_bus) * line_data(i, 5);
    Q_b(i, 2) = V(to_bus) * V(to_bus) * imag(Y) + V(to_bus) * V(from_bus) * (real(Y_ft) * sin(theta(to_bus) - theta(from_bus))...
        - imag(Y_ft) * cos(theta(to_bus) - theta(from_bus))) - V(to_bus) * V(to_bus) * line_data(i, 5);
end

%% Calculate Generator Power Outputs
results.gen = zeros(n_gen, 2);
for i = 1:n_gen
    generator_bus = gen_data(i, 1);
    connected_lines_from = find(line_data(:, 1) == generator_bus);
    connected_lines_to = find(line_data(:, 2) == generator_bus);
    results.gen(i, 1) = bus_data(generator_bus, 3) + sum(P_b(connected_lines_from, 1)) + sum(P_b(connected_lines_to, 2));
    results.gen(i, 2) = bus_data(generator_bus, 4) + sum(Q_b(connected_lines_from, 1)) + sum(Q_b(connected_lines_to, 2));
end

%% Prepare Results
results.bus = [V, theta];
results.gen = results.gen * baseMVA;
results.branch = [P_b, Q_b] * baseMVA;

%% Display Results
fprintf('Iterations: %d\n', iteration);
fprintf('***********************************Bus Data*******************************************\n');
fprintf(' Bus No.      Voltage Magnitude (pu)     Voltage Angle (degrees)\n');
for i = 1:size(bus_data, 1)
    fprintf('%6d%25.4f%26.4f\n', i, results.bus(i, 1), results.bus(i, 2) / pi * 180);
end
fprintf('\n***********************************Generator Power Data***********************************\n');
fprintf(' Bus No.      Active Power Output (MW)     Reactive Power Output (MVAR)\n');
for i = 1:size(gen_data, 1)
    fprintf('%6d%26.3f%26.3f\n', gen_data(i, 1), results.gen(i, 1), results.gen(i, 2));
end
fprintf('\n***********************************Line Power Flow***********************************\n');
fprintf(' From Bus      To Bus      From Bus Active Power Injection (MW)     From Bus Reactive Power Injection (MVAR)     To Bus Active Power Injection (MW)     To Bus Reactive Power Injection (MVAR)     Line Loss Active Power (MW)     Line Loss Reactive Power (MVAR)\n');
for i = 1:n_line
    fprintf('%8d%13d%30.3f%34.3f%37.3f%44.3f%30.3f%33.3f\n', line_data(i, 1), line_data(i, 2), ...
        results.branch(i, 1), results.branch(i, 3), -results.branch(i, 2), -results.branch(i, 4), ...
        results.branch(i, 1) + results.branch(i, 2), results.branch(i, 3) + results.branch(i, 4));
end
fprintf('\n***********************************Bus Admittance Matrix***********************************\n');
for i = 1:n_bus
    for j = 1:n_bus
        fprintf('%10.4f%+10.4fi', real(G_B(i, j)), imag(G_B(i, j)));
    end
    fprintf('\n');
end
fprintf('\n');
