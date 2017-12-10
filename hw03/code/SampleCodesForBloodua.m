% Sample codes for the calculation of blood ua at different SO2 level with the OMLC compiled molar extinction coef. of Hb and HbO2
%	http://omlc.org/spectra/hemoglobin/summary.html
%	
%															Edited by Meng-Lin Li, Ph. D., 11/19/2015
%															Dept. of Electrical Engineering, 
%															National Tsing Hua University	

% blood
% ua = (2.303) e (x g/liter)/(64,500 g Hb/mole) 
% ,where x is the number of grams per liter. A typical value of x for whole blood is x=150 g Hb/liter
x = 150;   % in g Hb/liter
MW_Hb = 64500;  % g Hb/mole, molecular weight of Hb
load e_HbO2_Hb; % molar extinction coef. table, (lambda, HbO2, Hb)

figure
semilogy(e_HbO2_Hb(:,1), e_HbO2_Hb(:,2), 'r-', e_HbO2_Hb(:,1), e_HbO2_Hb(:,3), 'b--', 'linewidth',2);
xlabel('Wavelenth(nm)');
ylabel('Molar extinction coefficient (1/cm/M)');
legend('HbO2', 'Hb');
legend('boxoff')

% SO2 = 95%
SO2 = 0.95;
ua_SO2_95 = 2.303*e_HbO2_Hb(:,2)*x*SO2/MW_Hb + 2.303*e_HbO2_Hb(:,3)*x*(1-SO2)/MW_Hb; % assume molecular weight of Hb and HbO2 are approximately the same

% SO2 = 70%
SO2 = 0.7;
ua_SO2_70 = 2.303*e_HbO2_Hb(:,2)*x*SO2/MW_Hb + 2.303*e_HbO2_Hb(:,3)*x*(1-SO2)/MW_Hb;

figure
plot(e_HbO2_Hb(:,1), ua_SO2_95, 'r-', e_HbO2_Hb(:,1), ua_SO2_70, 'b-', 'linewidth',2);
xlabel('Wavelenth(nm)');
ylabel('Absorption coefficient (1/cm)')
legend('Blood with SO_2=95%', 'Blood with SO_2=70%');
legend('boxoff')
axis tight



% SO2 = 100%
SO2 = 1;
ua_SO2_100 = 2.303*e_HbO2_Hb(:,2)*x*SO2/MW_Hb + 2.303*e_HbO2_Hb(:,3)*x*(1-SO2)/MW_Hb;

% SO2 = 0%
SO2 = 0;
ua_SO2_0 = 2.303*e_HbO2_Hb(:,2)*x*SO2/MW_Hb + 2.303*e_HbO2_Hb(:,3)*x*(1-SO2)/MW_Hb;

figure
plot(e_HbO2_Hb(:,1), ua_SO2_100, 'r-', e_HbO2_Hb(:,1), ua_SO2_0, 'b-', 'linewidth',2);
xlabel('Wavelenth(nm)');
ylabel('Absorption coefficient (1/cm)')
legend('Blood with SO_2=100%', 'Blood with SO_2=0%');
legend('boxoff')
axis tight






