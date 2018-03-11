%% Before running, set up the testbench cell name on line 8
% your decoder outputs should be labeled as Y0,Y1,...,Y7
% your decoder inputs shouled be labeled as A,B,C
%
% (c) 2018 Oscar Castaneda, Olalekan Afuye, Charles Jeon & Christoph Studer

% set up the name of your testbench cell
tb_name = 'decoder_mark2_sim';

% set up cds_srr function
addpath('/opt/cadence/MMSIM151/tools.lnx86/spectre/matlab/64bit');

% directory that contains the simulation outputs
directory = sprintf('%s/simulation/%s/spectre/schematic/psf', getenv('HOME'), tb_name);

% set up basic parameters
Vdd = 1.2; % define vdd

% define period
period_a = 800; % bit2 (MSB)
period_b = 400; % bit1
period_c = 200; % bit0 (LSB)

% get input signals
a = cds_srr(directory, 'tran-tran', 'A', 0);
b = cds_srr(directory, 'tran-tran', 'B', 0);
c = cds_srr(directory, 'tran-tran', 'C', 0);

% convert time into ps
t_ps = a.time*1e12;

% extract voltages of signals
a = a.V;
b = b.V;
c = c.V;

% get output signals and put them together in a table where the i-th
% column corresponds to the 'Y(i-1)' output
y_mtx = [];
for i=1:8
    signal_name = ['Y',int2str(i-1)];
    y = cds_srr(directory, 'tran-tran', signal_name, 0);
    y_mtx = [y_mtx y.V];
end

exp_y_mtx = zeros(size(y_mtx));
mydecoder_output = zeros(8,8);
exp_decoder_output = zeros(8,8);

% we sample at the midpoint for each output for the second period
t_ps_sample = 2*period_a + 5000 + 100 + (0:7)*200;

%% decoder output
err_flag = 0;
for i=1:8
    % find t_ps closest to t_ps_sample
    [~,t_ps_idx] = min(abs(t_ps-t_ps_sample(i)));
    
    % measure the outputs and declare 1 if it is greater than Vdd/2
    mydecoder_output(i,:) = y_mtx(t_ps_idx,:) > (Vdd/2);
    
    % create expected output waveform
    a_bits = (a > Vdd/2);
    b_bits = (b > Vdd/2);
    c_bits = (c > Vdd/2);
    vec_bits = [a_bits b_bits c_bits];
    exp_dec = bi2de(vec_bits,'left-msb');
    exp_y_mtx(:,i) = Vdd*(exp_dec == (i-1));
    
    % expected decoder output is given by:
    exp_decoder_output(i,exp_dec(t_ps_idx)+1) = 1;
    
    if (exp_decoder_output(i,:) ~= mydecoder_output(i,:))

    disp(['Expected output for input '...
        'A=' num2str(vec_bits(t_ps_idx,1)) ...
        ' B=' num2str(vec_bits(t_ps_idx,2)) ...
        ' C=' num2str(vec_bits(t_ps_idx,3)) ...
        ' is y' num2str(8-i) ...
        ' but measured output is y' num2str(exp_dec(t_ps_idx))...
        ]) 
    err_flag  = err_flag + 1;
    end
end

if err_flag == 0
    disp('Your decoder circuit has no errors :)')
end
    %
%% plots
figure(1)
set(gcf,'units','pixel');
set(gcf,'position',[0,200,800,600]);
% we have 8 outputs so total time period elapsed is 2*period_a = 1600ps
subplot(2,1,1);
plot(t_ps,c,'k',t_ps,b,'b',t_ps,a,'r','linewidth',4)
legend('c','b','a')
grid on
% note that our simulations start after 5ns so we add it to our window
xlim(2*period_a+5000 + [0 2*period_a])
title('inputs')
xlabel('time [ps]')
ylabel('input to decoder [V]')
set(gca,'fontsize',14)

% output from your circuit
subplot(2,1,2);
plot(t_ps,y_mtx,'linewidth',2)
legend('y0','y1','y2','y3','y4','y5','y6','y7')
grid on
% note that our simulations start after 5ns so we add it to our window
xlim(2*period_a + 5000 + [0 2*period_a])
title('outputs')
xlabel('time [ps]')
ylabel('output [V]')
set(gca,'fontsize',14)

% print to SVG (vector graphics file)
print -dsvg Lab2_decoder_IO.svg

figure(2)
set(gcf,'units','pixel');
set(gcf,'position',[800,200,800,600]);
for i=1:8
    subplot(2,4,i)
    plot(t_ps,y_mtx(:,i),'k',t_ps,exp_y_mtx(:,i),'r','linewidth',3);
    grid on
    legend('act.','exp.','location','south')
    xlim(2*period_a + 5000 + [0 2*period_a])
    ylim([-1 1.5])
    xlabel('time [ps]')
    ylabel('output [V]')
    title(['y' num2str(i-1)])
    set(gca,'fontsize',10)
end

% print to SVG (vector graphics file)
print -dsvg Lab2_decoder_expresp.svg