clear; clc; close all;

%% Parameters
n        = 256;
Bmag     = 5;
msbBits  = 3;
LWmax    = 104;
nFrames  = 2000;
EbN0dB   = 7;

sorterNames = ["bitonic", "odd_even", "bubble", "insertion", "merge"];

stats = struct();

for s = 1:numel(sorterNames)
    name = sorterNames(s);
    stats.(name).tieProbFrames = [];
    stats.(name).cmpCntFrames  = [];
    stats.(name).topOverlap    = zeros(nFrames,1);
end

%% Monte Carlo
for f = 1:nFrames

    % Replace this block with your real LLR generation if desired
    x = randi([0 1],1,n);
    bpsk = 1 - 2*x;

    EbN0 = 10^(EbN0dB/10);
    sigma = sqrt(1/(2*EbN0));

    y = bpsk + sigma*randn(1,n);
    llr = 2*y/(sigma^2);

    magTrue = abs(llr);

    % 5-bit quantized magnitude: 0...31
    magQ = round(magTrue);
    magQ = min(max(magQ,0),2^Bmag-1);

    [~, refIdx] = sort(magTrue,'ascend');

    for s = 1:numel(sorterNames)
        name = sorterNames(s);

        switch name
            case "bitonic"
                [idx, tieCnt, cmpCnt] = bitonic_sort_pa(magQ, msbBits, Bmag);

            case "odd_even"
                [idx, tieCnt, cmpCnt] = odd_even_sort_pa(magQ, msbBits, Bmag);

            case "bubble"
                [idx, tieCnt, cmpCnt] = bubble_sort_pa(magQ, msbBits, Bmag);

            case "insertion"
                [idx, tieCnt, cmpCnt] = insertion_sort_pa(magQ, msbBits, Bmag);

            case "merge"
                [idx, tieCnt, cmpCnt] = merge_sort_pa(magQ, msbBits, Bmag);
        end

        tieProb = tieCnt ./ max(cmpCnt,1);

        stats.(name).tieProbFrames = pad_and_append(stats.(name).tieProbFrames, tieProb);
        stats.(name).cmpCntFrames  = pad_and_append(stats.(name).cmpCntFrames, cmpCnt);

        K = min(LWmax,n);
        stats.(name).topOverlap(f) = numel(intersect(idx(1:K), refIdx(1:K))) / K;
    end
end

%% Print summary
fprintf("\n==== Pure MATLAB reduced-precision sorter framework ====\n");
fprintf("n = %d, Bmag = %d, MSB bits = %d, Eb/N0 = %.1f dB, frames = %d\n", ...
    n, Bmag, msbBits, EbN0dB, nFrames);

for s = 1:numel(sorterNames)
    name = sorterNames(s);

    meanTie = mean(stats.(name).tieProbFrames,1,'omitnan');
    meanOverlap = mean(stats.(name).topOverlap);

    fprintf("\n%s sorter\n", upper(name));
    fprintf("Mean top-LWmax overlap = %.2f %%\n", 100*meanOverlap);
    fprintf("Mean tie probability per stage/pass/level:\n");
    disp(meanTie);
end

%% Plot tie probability
figure;
hold on; grid on;

for s = 1:numel(sorterNames)
    name = sorterNames(s);
    meanTie = mean(stats.(name).tieProbFrames,1,'omitnan');
    plot(1:numel(meanTie), meanTie, '-o', 'DisplayName', name);
end

xlabel('Stage / pass / iteration / merge level');
ylabel('MSB tie probability');
title(sprintf('Tie probability, n=%d, MSB%d, Eb/N0=%.1f dB', n, msbBits, EbN0dB));
legend('Location','best');

%% Plot top-LWmax overlap
figure;
bar(categorical(sorterNames), ...
    arrayfun(@(s) mean(stats.(sorterNames(s)).topOverlap), 1:numel(sorterNames)));
grid on;
ylabel(sprintf('Top-%d overlap with true-valued ordering', LWmax));
title('Ordering preservation in decoder-relevant region');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [swap, isTie] = pa_compare(a, b, msbBits, Bmag)
    shift = Bmag - msbBits;

    aMSB = floor(double(a) / 2^shift);
    bMSB = floor(double(b) / 2^shift);

    isTie = (aMSB == bMSB);
    swap  = (aMSB > bMSB); % ascending order
end

function [vals, idx, isTie] = compare_exchange(vals, idx, i, j, msbBits, Bmag)
    [swap, isTie] = pa_compare(vals(i), vals(j), msbBits, Bmag);

    if swap
        tmp = vals(i); vals(i) = vals(j); vals(j) = tmp;
        tmp = idx(i);  idx(i)  = idx(j);  idx(j)  = tmp;
    end
end

function A = pad_and_append(A, row)
    row = row(:).';

    if isempty(A)
        A = row;
        return;
    end

    oldCols = size(A,2);
    newCols = numel(row);

    if newCols > oldCols
        A(:,oldCols+1:newCols) = NaN;
    elseif newCols < oldCols
        row(newCols+1:oldCols) = NaN;
    end

    A = [A; row];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bubble sorter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [idx, tieCnt, cmpCnt] = bubble_sort_pa(mag, msbBits, Bmag)

    n = numel(mag);
    vals = mag(:).';
    idx  = 1:n;

    tieCnt = zeros(1,n-1);
    cmpCnt = zeros(1,n-1);

    for pass = 1:n-1
        for i = 1:n-pass
            [vals, idx, isTie] = compare_exchange(vals, idx, i, i+1, msbBits, Bmag);

            tieCnt(pass) = tieCnt(pass) + isTie;
            cmpCnt(pass) = cmpCnt(pass) + 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Odd-even transposition sorter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [idx, tieCnt, cmpCnt] = odd_even_sort_pa(mag, msbBits, Bmag)

    n = numel(mag);
    vals = mag(:).';
    idx  = 1:n;

    tieCnt = zeros(1,n);
    cmpCnt = zeros(1,n);

    for phase = 1:n
        if mod(phase,2) == 1
            startIdx = 1;
        else
            startIdx = 2;
        end

        for i = startIdx:2:n-1
            [vals, idx, isTie] = compare_exchange(vals, idx, i, i+1, msbBits, Bmag);

            tieCnt(phase) = tieCnt(phase) + isTie;
            cmpCnt(phase) = cmpCnt(phase) + 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Insertion sorter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [idx, tieCnt, cmpCnt] = insertion_sort_pa(mag, msbBits, Bmag)

    n = numel(mag);
    vals = mag(:).';
    idx  = 1:n;

    tieCnt = zeros(1,n-1);
    cmpCnt = zeros(1,n-1);

    for i = 2:n
        j = i;

        while j > 1
            [swap, isTie] = pa_compare(vals(j-1), vals(j), msbBits, Bmag);

            tieCnt(i-1) = tieCnt(i-1) + isTie;
            cmpCnt(i-1) = cmpCnt(i-1) + 1;

            if swap
                tmp = vals(j-1); vals(j-1) = vals(j); vals(j) = tmp;
                tmp = idx(j-1);  idx(j-1)  = idx(j);  idx(j)  = tmp;
                j = j - 1;
            else
                break;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Merge sorter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [idx, tieCnt, cmpCnt] = merge_sort_pa(mag, msbBits, Bmag)

    n = numel(mag);
    vals = mag(:).';
    idx  = 1:n;

    nLevels = ceil(log2(n));
    tieCnt = zeros(1,nLevels);
    cmpCnt = zeros(1,nLevels);

    width = 1;
    level = 1;

    while width < n
        left = 1;

        while left <= n
            mid   = min(left + width - 1, n);
            right = min(left + 2*width - 1, n);

            if mid < right
                [vals(left:right), idx(left:right), t, c] = merge_blocks_pa( ...
                    vals(left:mid), idx(left:mid), ...
                    vals(mid+1:right), idx(mid+1:right), ...
                    msbBits, Bmag);

                tieCnt(level) = tieCnt(level) + t;
                cmpCnt(level) = cmpCnt(level) + c;
            end

            left = left + 2*width;
        end

        width = 2*width;
        level = level + 1;
    end
end

function [outVals, outIdx, tieCnt, cmpCnt] = merge_blocks_pa(valsA, idxA, valsB, idxB, msbBits, Bmag)

    i = 1;
    j = 1;
    k = 1;

    outVals = zeros(1,numel(valsA)+numel(valsB));
    outIdx  = zeros(1,numel(idxA)+numel(idxB));

    tieCnt = 0;
    cmpCnt = 0;

    while i <= numel(valsA) && j <= numel(valsB)
        [swap, isTie] = pa_compare(valsA(i), valsB(j), msbBits, Bmag);

        tieCnt = tieCnt + isTie;
        cmpCnt = cmpCnt + 1;

        if swap
            outVals(k) = valsB(j);
            outIdx(k)  = idxB(j);
            j = j + 1;
        else
            outVals(k) = valsA(i);
            outIdx(k)  = idxA(i);
            i = i + 1;
        end

        k = k + 1;
    end

    while i <= numel(valsA)
        outVals(k) = valsA(i);
        outIdx(k)  = idxA(i);
        i = i + 1;
        k = k + 1;
    end

    while j <= numel(valsB)
        outVals(k) = valsB(j);
        outIdx(k)  = idxB(j);
        j = j + 1;
        k = k + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bitonic sorter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [idx, tieCnt, cmpCnt] = bitonic_sort_pa(mag, msbBits, Bmag)

    n = numel(mag);

    if abs(log2(n) - round(log2(n))) > 0
        error('Bitonic sorter requires n to be a power of 2.');
    end

    vals = mag(:).';
    idx  = 1:n;

    nStages = log2(n);
    tieCnt = zeros(1,nStages);
    cmpCnt = zeros(1,nStages);

    for stage = 1:nStages
        k = 2^stage;

        for j = k/2:-1:1
            if bitand(j, j-1) ~= 0
                continue;
            end

            for i = 1:n
                ixj = bitxor(i-1,j) + 1;

                if ixj > i
                    ascending = mod(floor((i-1)/k),2) == 0;

                    [swap, isTie] = pa_compare(vals(i), vals(ixj), msbBits, Bmag);

                    tieCnt(stage) = tieCnt(stage) + isTie;
                    cmpCnt(stage) = cmpCnt(stage) + 1;

                    doSwap = (ascending && swap) || (~ascending && ~swap && vals(i) ~= vals(ixj));

                    if doSwap
                        tmp = vals(i); vals(i) = vals(ixj); vals(ixj) = tmp;
                        tmp = idx(i);  idx(i)  = idx(ixj); idx(ixj)  = tmp;
                    end
                end
            end
        end
    end

    % Final safety: sort according to reduced MSB key only.
    % This keeps output deterministic if the bitonic network direction logic
    % is not desired for non-hardware experiments.
    shift = Bmag - msbBits;
    key = floor(double(vals)/2^shift);
    [~, p] = sort(key,'ascend','stable');
    idx = idx(p);
end