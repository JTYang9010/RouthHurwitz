%% Routh-Hurwitz Memory Comparison
% Compare core storage of Full Routh, Rolling Routh, roots(), and eig().
% IMPORTANT: These are analytical storage estimates, NOT measured peak memory.
clear; clc; close all;

%% 1. Define polynomial orders
n = [10 12 20 50 100 200 500 1000];
m = ceil((n + 1)/2);
bytesPerDouble = 8;

%% 2. Core array storage models
% Full Routh table (n+1) by m
fullBytes = bytesPerDouble .* (n+1) .* m;
% Rolling algorithm's three main row buffers (prev2, prev1, current)
rollingBytes = bytesPerDouble .* 3 .* m;
% Dense companion matrix is a reasonable O(n^2) storage reference for roots.
% MATLAB roots() internals may use additional/different storage.
companionBytes = bytesPerDouble .* n.^2;
% Dense state matrix supplied to eig (INPUT size, not eig's workspace)
stateMatrixBytes = bytesPerDouble .* n.^2;

%% 3. Tabulate results
MiB = 1024^2;
T = table(n(:), fullBytes(:)/MiB, rollingBytes(:)/MiB, ...
    companionBytes(:)/MiB, stateMatrixBytes(:)/MiB, ...
    fullBytes(:)./rollingBytes(:), ...
    'VariableNames', {'Order','FullRouth_MiB','Rolling3Rows_MiB', ...
    'CompanionMatrix_MiB','StateMatrix_MiB','RouthStorageRatio'});
disp(T);

%% 4. Visualization: all four storage proxies
figure('Name','Memory Storage Comparison');
tiledlayout(1,2,'TileSpacing','compact');
nexttile;
plot(n, fullBytes/1024,'o-','LineWidth',1.6); hold on;
plot(n, rollingBytes/1024,'s-','LineWidth',1.6);
plot(n, companionBytes/1024,'^-','LineWidth',1.6);
plot(n, stateMatrixBytes/1024,'--','LineWidth',1.2);
grid on;
xlabel('Polynomial order n'); ylabel('Core array storage (KiB)');
title('Linear axes');
legend('Full Routh table','Rolling (3 rows)', ...
    'Companion matrix proxy for roots','Input A for eig', ...
    'Location','northwest');
nexttile;
loglog(n, fullBytes/1024,'o-','LineWidth',1.6); hold on;
loglog(n, rollingBytes/1024,'s-','LineWidth',1.6);
loglog(n, companionBytes/1024,'^-','LineWidth',1.6);
loglog(n, stateMatrixBytes/1024,'--','LineWidth',1.2);
grid on;
xlabel('Polynomial order n'); ylabel('Core array storage (KiB)');
title('Log-log axes');
legend('Full Routh O(n^2)','Rolling O(n)', ...
    'Companion matrix O(n^2)','State matrix O(n^2)', ...
    'Location','northwest');
sgtitle('Analytical Storage Estimates (not measured peak memory)');

%% 5. Memory reduction from Full to Rolling Routh
figure('Name','Routh Memory Reduction');
plot(n, fullBytes./rollingBytes, 'o-','LineWidth',1.7);
grid on;
xlabel('Polynomial order n'); ylabel('Full / Rolling core storage');
title('Routh storage reduction factor');

%% 6. Focus on the 12th-order experiment
idx12 = find(n == 12, 1);
fprintf('\n12th-order core arrays:\n');
fprintf('Full Routh table      : %.0f bytes\n', fullBytes(idx12));
fprintf('Rolling 3 row buffers : %.0f bytes\n', rollingBytes(idx12));
fprintf('Companion matrix      : %.0f bytes (proxy for roots)\n', companionBytes(idx12));
fprintf('Dense state matrix A  : %.0f bytes (input to eig)\n', stateMatrixBytes(idx12));
fprintf('Full/Rolling ratio    : %.2fx\n', fullBytes(idx12)/rollingBytes(idx12));

%% 7. Interpretation
% Full Routh: O(n^2) core storage.
% Rolling Routh: O(n) core storage.
% roots(): A dense companion-matrix formulation uses O(n^2) storage,
% but the plot is NOT MATLAB roots()'s actual memory use.
% eig(): Dense A alone takes O(n^2); internal workspace is unmeasured.
% Input construction, temporary arrays, and peak process memory are excluded.
