% Sample plotted output from simulations with [128,116] codes and 
% GRAND (hard detection), ORBGRAND (soft detection) MATLAB implementation, ORBGRAND (soft detection) 
% C implementation and ORBGRAND1 (soft detection). 

clear codes

n=128;
k=116;
n_CODES = 0;

% ORBGRAND (soft detection) C implementation
DECODER = 'ORBGRAND_C';

code.class = 'CAPOLAR';
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = DECODER;
code.LT = '--d';
code.color = 'b';
codes(n_CODES).code = code; 

code.class = 'PAC';
conv_code = [1 1 1 0 1 1 0 1];
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(bin2dec(num2str(conv_code))) '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = DECODER;
code.LT = '--s';
code.color = [0.1250, 0.6940, 0.9290];
codes(n_CODES).code = code; 

code.class = 'CRC';
poly='0x8f3';
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' poly '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = [DECODER ' ' poly];
code.LT = '--o';
code.color = 'r';
codes(n_CODES).code = code; 

code.class = 'RLC';
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = DECODER;
code.LT = '--x';
code.color = 'm';
codes(n_CODES).code = code; 


% ORBGRAND (soft detection) MATLAB implementation
DECODER = 'ORBGRAND_M'; 

code.class = 'CAPOLAR';
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = DECODER;
code.LT = '-d';
code.color = 'b';
codes(n_CODES).code = code; 

code.class = 'PAC';
conv_code = [1 1 1 0 1 1 0 1];
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(bin2dec(num2str(conv_code))) '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = DECODER;
code.LT = '-s';
code.color = [0.1250, 0.6940, 0.9290];
codes(n_CODES).code = code; 

code.class = 'CRC';
poly='0x8f3';
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' poly '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = [DECODER ' ' poly];
code.LT = '-o';
code.color = 'r';
codes(n_CODES).code = code; 

code.class = 'RLC';
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = DECODER;
code.LT = '-x';
code.color = 'm';
codes(n_CODES).code = code; 


% ORBGRAND HW 
DECODER = 'ORBGRAND-BASELINE';

code.class = 'CAPOLAR';
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = DECODER;
code.LT = '-.d';
code.color = 'b';
codes(n_CODES).code = code; 

code.class = 'PAC';
conv_code = [1 1 1 0 1 1 0 1];
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(bin2dec(num2str(conv_code))) '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = DECODER;
code.LT = '-.s';
code.color = [0.1250, 0.6940, 0.9290];
codes(n_CODES).code = code; 

code.class = 'CRC';
poly='0x8f3';
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' poly '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = [DECODER ' ' poly];
code.LT = '-.o';
code.color = 'r';
codes(n_CODES).code = code; 

code.class = 'RLC';
n_CODES=n_CODES+1;
filename = ['../RESULTS/' DECODER '_' code.class '_' num2str(n) '_' num2str(k) '_1.mat'];
load(filename,'code');
code.decoder = DECODER;
code.LT = '-.x';
code.color = 'm';
codes(n_CODES).code = code; 

make_fig(codes,1)

function make_fig(codes,fig_no)

    FONT=16;
    MS=8;
    LW=2;
    num_codes = length(codes);
    
    figure(fig_no)
    clf
    % BLER
    subplot(1,3,1)
    hold on
    for ii=1:num_codes
        code_info = [codes(ii).code.decoder ', ' codes(ii).code.class ' [' num2str(codes(ii).code.n) ',' num2str(codes(ii).code.k) '], R=' num2str(codes(ii).code.k/codes(ii).code.n,'%.2f')];
        plot(codes(ii).code.ebn0,codes(ii).code.BLER,codes(ii).code.LT,'displayname',code_info,'color',codes(ii).code.color,'LineWidth',LW, 'MarkerSize',MS)
    end
    hold off
    legend('show','Location','NorthEast');
    ylabel('BLER')
    xlabel('Eb/N0')
    grid 'on'
    xl = xlim;
    xlim([xl(1)-1 xl(2)+1])
    set(gca, 'YScale', 'log')
    yl=ylim; 
    ylim([10^floor(log10(yl(1))) 10^0])
    set(gca,'FontSize',FONT)
    
    % BER
    subplot(1,3,2)
    hold on
    for ii=1:num_codes
        code_info = [codes(ii).code.decoder ', ' codes(ii).code.class ' [' num2str(codes(ii).code.n) ',' num2str(codes(ii).code.k) '], R=' num2str(codes(ii).code.k/codes(ii).code.n,'%.2f')];
        plot(codes(ii).code.ebn0,codes(ii).code.BER,codes(ii).code.LT,'displayname',code_info,'color',codes(ii).code.color,'LineWidth',LW, 'MarkerSize',MS)
    end
    hold off
    legend('show','Location','NorthEast');
    ylabel('BER')
    xlabel('Eb/N0')
    grid 'on'
    xl = xlim;
    xlim([xl(1)-1 xl(2)+1])
    set(gca, 'YScale', 'log')
    yl=ylim; 
    ylim([10^floor(log10(yl(1))) 10^0])
    set(gca,'FontSize',FONT)

    % Average number of codebook queries
    subplot(1,3,3)
    hold on
    for ii=1:num_codes
        code_info = [codes(ii).code.decoder ', ' codes(ii).code.class ' [' num2str(codes(ii).code.n) ',' num2str(codes(ii).code.k) '], R=' num2str(codes(ii).code.k/codes(ii).code.n,'%.2f')];
        plot(codes(ii).code.ebn0,codes(ii).code.EG,codes(ii).code.LT,'displayname',code_info,'color',codes(ii).code.color,'LineWidth',LW, 'MarkerSize',MS)
    end
    hold off
    legend('show','Location','NorthEast');
    ylabel('Average number of code-book queries')
    xlabel('Eb/N0')
    grid 'on'
    xl = xlim;
    xlim([xl(1)-1 xl(2)+1])
    set(gca, 'YScale', 'log')
    set(gca,'FontSize',FONT)

 
end

