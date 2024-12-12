% Parameters
N = 500;          % Length of the taper (number of points)
NW = 4;           % Time-bandwidth product
K = 2*NW - 1;     % Number of tapers (2*NW - 1 is typical)

% Generate DPSS tapers
[tapers, eigenvalues] = dpss(N, NW, K);

% Plot the tapers
figure;
for i = 1:K
    subplot(K, 1, i); 
    plot(tapers(:, i), 'LineWidth', 1.5);
    title(['Taper ', num2str(i), ', Eigenvalue: ', num2str(eigenvalues(i))]);
    xlabel('Time Index'); ylabel('Amplitude');
    grid on;
end

% Explanation of messy taper
figure;
plot(1:K, eigenvalues, '-o', 'LineWidth', 1.5);
xlabel('Taper Index'); ylabel('Eigenvalue');
title('Energy Concentration (Eigenvalues) for DPSS Tapers');
grid on;

% Parameters
N = 500;          % Length of the sequence
NW = 4;           % Time-bandwidth product
W = NW / N;       % Half-bandwidth (normalized)

% Create the Toeplitz matrix
n = 0:(N-1);
C = zeros(N, N);  % Initialize the matrix
for i = 1:N
    for j = 1:N
        % Sinc kernel defining the DPSS problem
        C(i, j) = sin(2 * pi * W * (i - j)) / (pi * (i - j));
        if i == j
            C(i, j) = 2 * W; % Correct diagonal element
        end
    end
end

% Solve the eigenvalue problem
[eigVectors, eigValuesMatrix] = eig(C);

% Extract eigenvalues (diagonal elements of eigValuesMatrix)
eigenvalues = diag(eigValuesMatrix);

% Sort eigenvalues and corresponding eigenvectors
[eigenvalues, idx] = sort(eigenvalues, 'descend');
eigVectors = eigVectors(:, idx);

% Display eigenvalues
disp('Eigenvalues:');
disp(eigenvalues(1:10)); % Display the first 10 eigenvalues

% Plot tapers
figure;
for i = 1:10
    subplot(10, 1, i);
    plot(eigVectors(:, i), 'LineWidth', 1.5);
    title(['Taper ', num2str(i), ', Eigenvalue: ', num2str(eigenvalues(i))]);
    xlabel('Time Index'); ylabel('Amplitude');
    grid on;
end
