% Transfer plot

load('impulsieTransit.mat')

figure
plot3(-mu,0,0,'^k')
hold on
plot3(1-mu,0,0,'^k')
plot3(LP(1,1),0,0,'xk')
grid on
clf
plot3(-mu,0,0,'^k')
hold on
plot3(1-mu,0,0,'vk')
plot3(LP(1,1),0,0,'xk')
plot3(rr_halo1(:,1),rr_halo1(:,2),rr_halo1(:,3),'-b','LineWidth',2)
plot3(manif(best_rrManif_s_p_indx).rr_s_p(:,1),manif(best_rrManif_s_p_indx).rr_s_p(:,2),...
    manif(best_rrManif_s_p_indx).rr_s_p(:,3),'-c','LineWidth',2);
plot3(rr_impTrans(:,1),rr_impTrans(:,2),rr_impTrans(:,3),'-r','LineWidth',2);

legend('Earth','Moon','L1','Halo','Stable manif.','Impulsive Transfer');

grid on; axis equal;
xlabel('x'); ylabel('y'); zlabel('z');