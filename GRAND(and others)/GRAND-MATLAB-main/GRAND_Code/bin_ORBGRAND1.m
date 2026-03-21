% 1-line ORGBRAND (soft detection) 

% From 
% K. R. Duffy, W. An, and M. Medard, "Ordered reliability bits guessing 
% random additive noise decoding,” IEEE Trans. Signal Process., vol. 70, 
% pp. 4528-4542, 2022.
% Uses the Landslide algorithm from that paper and the light-weight
% 1-line implementation introduced in 
% K. Galligan, M. Médard, K. R. Duffy, "Block turbo decoding with ORBGRAND"
% arXiv:2207.11149, 2022.

% 1-line ORBGRAND differs from Basic ORBGRAND by considering a dynamically
% determined intercept in its statistical model. 

% ORBGRAND with 1-line fitting provides an observable, but
% limited, improvement over the basic version because their only difference is that the
% 1-line version starts from the quantized value of 𝐿1 instead of the origin in the basic
% version.


% Inputs:
%   n               - code length
%   H               - Parity check matrix or CRC function
%   max_query       - Maximum number of code-book queries to abandonment
%   y_soft          - Channel soft information
%
% Outputs:
%   y_decoded       - Decoded ML codeword
%   err_vec         - Putative noise
%   n_guesses       - Number of guesses performed
%   abandoned       - 1 if abandoned, 0 if found a codeword

function [y_decoded,err_vec,n_guesses,abandoned] = bin_ORBGRAND1(H,max_query,y_soft)

    % Hard demodulate
    y_demod = (y_soft<0);
    n=length(y_demod);

    n_guesses = 1;
    err_vec = zeros(1,n);

    %Abandon default values
    y_decoded = -1*ones(size(y_demod));

    % First query is demodulated string     
    Hy = mod(H*y_demod',2);
    if Hy==zeros(size(Hy))
        y_decoded = y_demod;
        abandoned = 0;
        return;
    end
    
    % If the demodulated string is not in the codebook, decode.
    % Bit reliability
    reliability = abs(y_soft);
    [L, ind_order] = sort(reliability,'ascend');

    % Inverse sort order
    inv_perm = zeros(1,length(ind_order));
    for ii=1:length(ind_order)
           inv_perm(ind_order(ii))=ii;    
    end

    % A simple and effective method to pick c is to fit a line
    % through the points (1, L1) and (n/2, Ln/2), where n is the code length.

    % Slope
    beta = (L(round(n/2))-L(1))/(round(n/2)-1);
    % Intercept
    c = max(round(L(1)/beta-1),0);

    % This gives the best estimate of (1, L1), the least reliable and
    % most important bit, and accurately approximates the remaining
    % values if they follow a line-like distribution

    % This is the H columns reordered to put in ML order
    test_H = H(:,ind_order);

    %Noise sequences are generated in order of their total weight wt
    % For each wt we iterate over all pairs of non-negative integers
    % (wH, wL) such that wT = cwH + wL, and the landslide
    % algorithm generates all noise sequences for each pair

    % Total starting weight
    wt=c+1; %minimum possible weight

    while n_guesses<max_query && wt<=c*n+n*(n+1)/2
        % Hamming weight
        w=max(1,ceil((1+2*(n+c)-sqrt((1+2*(n+c))^2-8*wt))/2));
        while w<=n
            % Logistic weight
            W = wt-c*w;
            if W<w*(w+1)/2
                break; % invalid pair
            else
                % Make error vectors
                % Internally converts W and n to W' and n'.
                noise_locations = landslide(W,w,n); 
                % For each error vector
                for jj=1:size(noise_locations,1)
                   n_guesses = n_guesses +1;
                    err_vec =zeros(1,n);
                    err_vec(noise_locations(jj,:))=1;
                    if (Hy == mod(test_H*err_vec',2))
                        err_vec = err_vec(inv_perm);
                        y_decoded = mod(y_demod-err_vec,2);
                        abandoned = 0;
                        return;
                    end
                end
            end
            % Increment Hamming weight
            w=w+1;
        end
        wt=wt+1;
    end
    % If we max out on queries or total weight
    abandoned = 1;
    err_vec = zeros(size(y_demod));
end

