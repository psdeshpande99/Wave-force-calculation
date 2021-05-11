figure('Name','lever arm','NumberTitle','off','Units','normalized','Position',[0 0 1 1]);

x = [1 2 3 4 5];
hig = [h1 h2 h3 h4 h5]; 
w1 = 0.5; 
bar(hig,w1,'FaceColor',[0.2 0.2 0.5])
hck = [h1c h2c h3c h4c h5c];
w2 = .25;
hold on
bar(hck,w2,'FaceColor',[0 0.7 0.7])
hold off
grid on
ylabel('Lever Arm(m)')
legend({'Igawaki','Chakrabart'},'Location','northwest')
ax = gca;
ax.XTick = [1 2 3 4 5]; 
ax.XTickLabels = {'No current','Uniform(+)','Uniform(-)','Varying(+)','Varying(-)'};
ax.XTickLabelRotation = 45;

figure('Name','Fmax','NumberTitle','off','Units','normalized','Position',[0 0 1 1]);

x = [1 2 3 4 5];
hig = [Ftmax_nocurrent, Ftmax_uniform_current_supporting,Ftmax_unform_oppo,Ftmax_varying_current_supporting,Ftmax_varying_current_opposing]; 
w1 = 0.5; 
bar(hig,w1,'FaceColor',[0.2 0.2 0.5])
hck = [Ftmax_nocurrentc, Ftmax_uniform_current_supportingc,Ftmax_unform_oppoc,Ftmax_varying_current_supportingc,Ftmax_varying_current_opposingc];
w2 = .25;
hold on
bar(hck,w2,'FaceColor',[0 0.7 0.7])
hold off
grid on
ylabel('Fmax(N)')
legend({'Igawaki','Chakrabart'},'Location','northwest')
ax = gca;
ax.XTick = [1 2 3 4 5]; 
ax.XTickLabels = {'No current','Uniform(+)','Uniform(-)','Varying(+)','Varying(-)'};
ax.XTickLabelRotation = 45;



figure('Name','Mmax','NumberTitle','off','Units','normalized','Position',[0 0 1 1]);

x = [1 2 3 4 5];
hig = [Mtmax_nocurrent, Mtmax_uniform_current_supporting,Mtmax_unform_oppo,Mtmax_varying_current_supporting,Mtmax_varying_current_opposing]; 
w1 = 0.5; 
bar(hig,w1,'FaceColor',[0.2 0.2 0.5])
hck = [Mtmax_nocurrentc, Mtmax_uniform_current_supportingc,Mtmax_unform_oppoc,Mtmax_varying_current_supportingc,Mtmax_varying_current_opposingc];
w2 = .25;
hold on
bar(hck,w2,'FaceColor',[0 0.7 0.7])
hold off
grid on
ylabel('Mmax(Nm)')
legend({'Igawaki','Chakrabart'},'Location','northwest')
ax = gca;
ax.XTick = [1 2 3 4 5]; 
ax.XTickLabels = {'No current','Uniform(+)','Uniform(-)','Varying(+)','Varying(-)'};
ax.XTickLabelRotation = 45;