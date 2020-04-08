clear
clc


files = {'cylinder_399K_deg4_4cpu.csv', 'cylinder_399K_deg4_8cpu.csv','cylinder_399K_deg4_12cpu.csv','cylinder_399K_deg4_24cpu.csv',...
         'cylinder_399K_deg4_48cpu.csv','cylinder_399K_deg4_96cpu.csv' };

time = zeros(size(files));
flops = time;
for i=1:size(files,2)
    filename = files{i};
    df = readtable(filename,'ReadVariableNames',true);
    time(i) = sum(df{:,4});
    flops(i) = sum(df{:,8});       
end

MFlops_perTime = 10^(-6)*flops./time;
Haswell_numRank_deg_4 = [4 8 12 24 48 96];

prb_title = 'Linear 3D Elasticity (degree 4) \nu=0.3 , E=1e6, 399K elements, 75M dofs';
hardware_sky = 'Skylake (Xeon Gold 6126)';
semilogx(Haswell_numRank_deg_4,MFlops_perTime, 'r-o');
legend(hardware_sky,'Location','northwest');
title(prb_title);
ylabel('MFops')
xlabel('Number of MPI ranks')
L = MFlops_perTime(:)/2400*100;
text(Haswell_numRank_deg_4,MFlops_perTime, num2str(L(1:length(Haswell_numRank_deg_4))), 'VerticalAlignment','bottom');
% L = strsplit(sprintf('(%.1f, %.1f)\n', [Haswell_numRank_deg_4(:) MFlops_perTime(:)]'), '\n');
% text(Haswell_numRank_deg_4, MFlops_perTime, L(1:length(Haswell_numRank_deg_4)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom');
% ylim([50 160])
% xlim([0 110])
 xticks([4 8 12 24 48 96])
 xticklabels({'4', '8','12', '24', '48', '96'})
grid on