% 1. DEFINE PARAMETERS
n = 1024;
k = 512;
m = n - k;

fprintf('Generating LDPC (1024, 512) matrices...\n');

% 2. CREATE SYSTEMATIC H matrix [P | I]
% Using a random P and an Identity matrix I ensures H is full rank.
P = randi([0 1], m, k);
I = eye(m);
H = [P I];

% 3. CREATE SYSTEMATIC G matrix [I | P']
% In GF(2), if H = [P | I], then G = [I | P^T]
G = [eye(k) P'];

% 4. VERIFY G * H' = 0 (Optional but recommended)
check = mod(G * H', 2);
if any(check(:))
    error('Verification failed: G*H'' is not zero.');
else
    fprintf('Verification successful: G*H'' = 0.\n');
end

% 5. SAVE TO MAT FILE
% We save H as sparse (efficient) and G as full (standard for generator matrices)
H = sparse(H);
G = full(G); 
save('ldpc_1024_512.mat', 'G', 'H');
fprintf('File "ldpc_1024_512.mat" created successfully.\n\n');

% ---------------------------------------------------------
% 6. YOUR DROP-IN CODE STARTS HERE
% ---------------------------------------------------------
S = load('ldpc_1024_512.mat');
G = S.G; 
H = S.H;
[k, n] = size(G);
code_label = 'LDPC_1024_512';

fprintf('Loaded: %s (k=%d, n=%d)\n', code_label, k, n);
