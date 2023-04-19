x_ax = linspace(dx/2,R-dx/2,R/dx);
Th_ax = DX_mat(:,2).*(180/pi);
Sr_ax = DX_mat(:,2).*transpose((dx/2):dx:(R-dx/2));
Sc_ax = DX_mat(:,3);
s_ax = DX_mat(:,4);
phi_ax = DX_mat(:,5).*(180/pi);
sc_ax = DX_mat(:,5).*transpose((dx/2):dx:(R-dx/2));
H_ax = DX_mat(:,6).*60;
h_ax = DX_mat(:,7).*60;
d_ax = DX_mat(:,8);

xline1 = 0; xline2 = 0;
for ii = 2:length(FMLM)
    if FMLM(ii,2) == 1 && FMLM(ii-1,2) < 1
        xline1 = FMLM(ii,1) - dx/2;
    end
    if FMLM(ii,2) == 2 && FMLM(ii-1,2) < 2
        xline2 = FMLM(ii,1) - dx/2;
    end
end

set(groot, 'defaultLineLineWidth',1);
figure('Name', strcat('parameters (all) - ',num2str(rho_0),' - ',num2str(scen))) %('DefaultAxesColorOrder', [0 0 1; 1 0 0])
grid on;
%yticks(0:5);
plot(x_ax,Sr_ax,'k-', x_ax,Sc_ax,'r-', x_ax,s_ax,'k--',x_ax,sc_ax,'r--');
ylim([0, 4]);
if xline1 > 0, xline(xline1,'--',{'Fixed','feeder'},'LabelOrientation','horizontal'); end
if xline2 > 0, xline(xline2,'--',{'DR','feeder'},'LabelOrientation','horizontal'); end
%title(['Spacing and Headway | r = ' num2str(r_opt) ' km' ' | R = ' num2str(R) ' km']); 
xlabel('Distance from center x (km)');
ylabel('Spacing [km]'); %left y-axis
hold on
yyaxis right
ylim([0, 10]);%ylim([0, max([max(H_ax),max(h_ax)])+1]);
plot(x_ax,H_ax,'b-.',x_ax,h_ax,'b:');
ylabel('Headway [min]'); %right y-axis
legend({'S_r(x)','S_c(x)','s(x)','s_c(x)','H(x)','h(x)'},'Location','northoutside','NumColumns',6)
%saveas(gcf,strcat('struct - ',num2str(rho_0),' - ',num2str(f),' - ',num2str(v_p),'.png'));
hold off

% figure('Name', strcat('parameters - ',num2str(rho_0),' - ',num2str(f)))
% subplot(3,1,1);
% plot(x_ax, Th_ax,'k-', x_ax, phi_ax, 'r--')
% ylim([0, max([45,max(Th_ax)+10])]);
% ylabel('Angular spacing [°]'); %left y-axis
% legend({'Theta_r(x)','phi(x)'},'Location','northwest','NumColumns',2)
% subplot(3,1,2);
% plot(x_ax, Sc_ax,'r-', x_ax, s_ax, 'k--')
% ylim([0, max(Sc_ax)+1]);
% ylabel('Linear spacing [km]'); %left y-axis
% legend({'S_c(x)','s(x)'},'Location','northwest','NumColumns',2)
% subplot(3,1,3);
% plot(x_ax, H_ax,'b-.', x_ax, h_ax, 'b:')
% ylim([0, (max(H_ax)+1)]);
% xlabel('Distance from center x (km)');
% ylabel('Headway [min]'); %left y-axis
% legend({'H(x)','h(x)'},'Location','northwest','NumColumns',2)
% saveas(gcf,strcat('param - ',num2str(rho_0),' - ',num2str(f),' - ',num2str(v_p),'.png'));

% 
% figure('Name', strcat('pie - ',num2str(rho_0),' - ',num2str(f)))
% explode = [1,1,1,1,1,1];
% pg = pie([Z_L,Z_V,Z_M,Z_A,Z_W,Z_T],explode);
% %title(['Total cost Z = ' num2str(round(Z/1000)) ' k' char(8364)]);
% labels = {'Z_L','Z_V','Z_M','Z_A','Z_W','Z_T'};
% legend(labels,'Location','southeastoutside','Orientation','horizontal','NumColumns',1);
% saveas(gcf,strcat('pie - ',num2str(rho_0),' - ',num2str(f),' - ',num2str(v_p),'.png'));

% figure
% plot(x_ax,rho_0.*exp(-gamma*x_ax));
% title('demand plot');
% hold on
% plot(x_ax,32.*ones(1,R/dx),'--');
% xlabel('Distance from center x (km)');
% ylabel('trips/km^2-h'); %left y-axis