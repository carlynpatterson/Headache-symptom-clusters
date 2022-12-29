% Headache-associated symptoms multiple correspondence analysis across all headache participants

load PfizerHAdata

% select participants 6 - 17 years with a pedmidas score and headache
% severity score
HA = data(data.age>=6 & data.age<18 & ~isnan(data.p_pedmidas_score) & (data.p_sev_overall=='mild' | data.p_sev_overall=='mod' | data.p_sev_overall=='sev'),:);
clear data

% calculate age at the time of filling out the form in days
HA.ageDays = datenum(HA.visit_dt)-datenum(HA.dob);

% create a variable for the number of prescription preventive medications
% tried
HA.num_prescr_prevent = sum(table2array(HA(:,[318:321 323:332 335:343])),2);

HA.allodynia = sum(table2array(HA(:,817:820)),2); % clinician entered data
HA.allodynia(HA.allodynia>0) = 1;

HA.pulsate = sum(table2array(HA(:,[123 124 133])),2);
HA.pulsate(HA.pulsate>1) = 1;
HA.pressure = sum(table2array(HA(:,[126:128 131])),2);
HA.pressure(HA.pressure>1) = 1;
HA.neuralgia = sum(table2array(HA(:,[125 129 130 132])),2);
HA.neuralgia(HA.neuralgia>1) = 1;

HA.severity_grade = ones(height(HA),1);
HA.severity_grade(HA.p_sev_overall=='mod') = 2;
HA.severity_grade(HA.p_sev_overall=='sev') = 3;

HA.pedmidas_grade = zeros(height(HA),1);
HA.pedmidas_grade(HA.p_pedmidas_score<10) = 0;
HA.pedmidas_grade(HA.p_pedmidas_score>10 & HA.p_pedmidas_score<=30) = 1;
HA.pedmidas_grade(HA.p_pedmidas_score>30 & HA.p_pedmidas_score<=50) = 2;
HA.pedmidas_grade(HA.p_pedmidas_score>50) = 3;

HA.pedmidas_grade2 = categorical(HA.pedmidas_grade);
HA.pedmidas_grade2(HA.p_pedmidas_score<10) = 'none';
HA.pedmidas_grade2(HA.p_pedmidas_score>10 & HA.p_pedmidas_score<=30) = 'mild';
HA.pedmidas_grade2(HA.p_pedmidas_score>30 & HA.p_pedmidas_score<=50) = 'mod';
HA.pedmidas_grade2(HA.p_pedmidas_score>50) = 'sev';

HA.freq_bad = NaN*ones(height(HA),1);
HA.freq_bad (HA.p_fre_bad=='never') = 1;
HA.freq_bad (HA.p_fre_bad=='1mo') = 2;
HA.freq_bad (HA.p_fre_bad=='1to3mo') = 3;
HA.freq_bad (HA.p_fre_bad=='1wk') = 4;
HA.freq_bad (HA.p_fre_bad=='2to3wk') = 5;
HA.freq_bad (HA.p_fre_bad=='3wk') = 6;
HA.freq_bad (HA.p_fre_bad=='daily') = 7;
HA.freq_bad (HA.p_fre_bad=='always') = 8;

HA.ageYr = round(HA.age);

%% Calculate MCA for all entries

% MCA (only associated sx)
haASx=table2array(HA(:,[245 247 250:255 257 259:260 208 209]));
var_HAaSx=char('nausea','vomiting','light sensitivity','smell sensitivity','sound sensitivity','light headed',...
    'room spinning','balance problems','ear ringing','neck pain','difficulty thinking','blurry vision',...
    'double vision');

% reconfigure data so they can be used in the MCA function, mcorran. All variables need to be converted to binary

binary_hx=cell(size(haASx));
binary_struct=NaN*ones(size(haASx,2),1);
for x=1:size(haASx,2)
    temp=haASx(:,x);
    outcome=unique(temp);
        for y=1:size(haASx,1)
            binary_struct(x,:)=2;
                switch temp(y)
                    case outcome(1)
                        binary_hx{y,x}=[1 0];
                    case outcome(2)
                        binary_hx{y,x}=[0 1];
                end
        end
end

% concatonate each subjects binary outcomes
binary_Hx=NaN*ones(size(binary_hx,1),size(var_HAaSx,1)*2);
temp=[];
for x=1:size(binary_hx,1)
    for y=1:size(binary_hx,2)
        temp=cat(2,temp,cell2mat(binary_hx(x,y)));
    end
    binary_Hx(x,:)=temp;
    temp=[];
end

% Calculate MCA
[~,~,~,~,~,MCA_corrHAaSx,sx_scores] = mcorran3(binary_Hx,binary_struct','var_names',var_HAaSx);

% Calculate MCA scores
MCA_no=3;
MCA_score_HAaSx=NaN*ones(size(binary_Hx,1),MCA_no);
for x=1:size(binary_Hx,1)
    for y=1:MCA_no
        temp1=binary_Hx(x,:);
        temp2=MCA_corrHAaSx(:,y);
        r=temp1*temp2;
        MCA_score_HAaSx(x,y)=r;
    end
end

HA.MCA1_HAaSx = MCA_score_HAaSx(:,1);
HA.MCA2_HAaSx = MCA_score_HAaSx(:,2).*-1;

%% plot symptom presence MCA map with dots representing the relative frequency of a symptom
figure
hold on
set(gca,'TickDir','out'); set(gca,'Box','off');
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XLim = [-1 2]; ax.YLim = [-1 2];


for x = 1:size(var_HAaSx,1)
    HAa = haASx(:,x);
    HAa_prct = length(HAa(HAa==1))./length(HAa);
    plot(sx_scores(1,x),sx_scores(2,x).*-1,'ok','MarkerSize',round(100*HAa_prct))
    text(sx_scores(1,x)+0.01,-1*sx_scores(2,x)+0.01,var_HAaSx(x,:))
end

%% Headache characteristics

figure

subplot(2,2,1)
histogram(HA.p_sev_overall,'Normalization','probability')
set(gca,'TickDir','out'); set(gca,'Box','off');
title('usual headache severity')

subplot(2,2,2)
histogram(HA.p_pedmidas_score,'Normalization','probability')
set(gca,'TickDir','out'); set(gca,'Box','off');
title('pedmidas')

subplot(2,2,3)
histogram(HA.p_fre_bad,'Normalization','probability')
set(gca,'TickDir','out'); set(gca,'Box','off');
title('frequency of bad headaches')

subplot(2,2,4)
histogram(HA.pulsate,'Normalization','probability','FaceColor','k')
hold on
histogram(HA.pressure,'Normalization','probability','FaceColor','b')
histogram(HA.neuralgia,'Normalization','probability','FaceColor','r')
set(gca,'TickDir','out'); set(gca,'Box','off');
title('headache quality')


%% Compare MCA and headache burden metrics

figure
subplot(2,3,1)
hold on
mild1 = HA.MCA1_HAaSx(HA.severity_grade==1);
mod1 = HA.MCA1_HAaSx(HA.severity_grade==2);
sev1 = HA.MCA1_HAaSx(HA.severity_grade==3);
jtr = (rand(size(mild1))-0.5).*0.1;
plot(ones(size(mild1))+jtr,mild1,'.k','MarkerSize',8)
errorbar(1,prctile(mild1,50),abs(diff(prctile(mild1,[25 50]))),abs(diff(prctile(mild1,[50 75]))),'+r','MarkerSize',12,'MarkerEdgeColor','r')
jtr = (rand(size(mod1))-0.5).*0.1;
plot(2*ones(size(mod1))+jtr,mod1,'.k','MarkerSize',8)
errorbar(2,prctile(mod1,50),abs(diff(prctile(mod1,[25 50]))),abs(diff(prctile(mod1,[50 75]))),'+r','MarkerSize',12,'MarkerEdgeColor','r')
jtr = (rand(size(sev1))-0.5).*0.1;
plot(3*ones(size(sev1))+jtr,sev1,'.k','MarkerSize',8)
errorbar(3,prctile(sev1,50),abs(diff(prctile(sev1,[25 50]))),abs(diff(prctile(sev1,[50 75]))),'+r','MarkerSize',12,'MarkerEdgeColor','r')
[rho,pval] = corr([HA.severity_grade HA.MCA1_HAaSx],'Type','Spearman');
title(sprintf('severity grade, Spearmans Rho = %1.2f, p = %1.1e',[rho(1,2) pval(1,2)]))
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XTick = 1:3; ax.XTickLabel = {'mild','moderate','severe'}; ylabel('MCA1 factor loading'); ax.XLim = [0 4];

% bootstrap analysis of headache severity vs. MCA rho value

bootstat = bootstrp(6662,@(a,b)corr([a b],'Type','Spearman'),HA.severity_grade,HA.MCA1_HAaSx);
fprintf('HA severity vs. MCA1, Rho = %1.2f; 95CI %1.2f - %1.2f \n',prctile(bootstat(:,2),[50 2.5 97.5]))


bootstat = bootstrp(6662,@(a,b)corr([a b],'Type','Spearman'),HA.severity_grade,HA.MCA2_HAaSx);
fprintf('HA severity vs. MCA2, Rho = %1.2f; 95CI %1.2f - %1.2f \n',prctile(bootstat(:,2),[50 2.5 97.5]))

subplot(2,3,2)
hold on
none2 = HA.MCA1_HAaSx(HA.pedmidas_grade==0);
mild2 = HA.MCA1_HAaSx(HA.pedmidas_grade==1);
mod2 = HA.MCA1_HAaSx(HA.pedmidas_grade==2);
sev2 = HA.MCA1_HAaSx(HA.pedmidas_grade==3);
jtr = (rand(size(none2))-0.5).*0.1;
plot(0*ones(size(none2))+jtr,none2,'.k','MarkerSize',8)
errorbar(0,prctile(none2,50),abs(diff(prctile(none2,[25 50]))),abs(diff(prctile(none2,[50 75]))),'+r','MarkerSize',12,'MarkerEdgeColor','r')
jtr = (rand(size(mild2))-0.5).*0.1;
plot(ones(size(mild2))+jtr,mild2,'.k','MarkerSize',8)
errorbar(1,prctile(mild2,50),abs(diff(prctile(mild2,[25 50]))),abs(diff(prctile(mild2,[50 75]))),'+r','MarkerSize',12,'MarkerEdgeColor','r')
jtr = (rand(size(mod2))-0.5).*0.1;
plot(2*ones(size(mod2))+jtr,mod2,'.k','MarkerSize',8)
errorbar(2,prctile(mod2,50),abs(diff(prctile(mod2,[25 50]))),abs(diff(prctile(mod2,[50 75]))),'+r','MarkerSize',12,'MarkerEdgeColor','r')
jtr = (rand(size(sev2))-0.5).*0.1;
plot(3*ones(size(sev2))+jtr,sev2,'.k','MarkerSize',8)
errorbar(3,prctile(sev2,50),abs(diff(prctile(sev2,[25 50]))),abs(diff(prctile(sev2,[50 75]))),'+r','MarkerSize',12,'MarkerEdgeColor','r')
[rho,pval] = corr([HA.pedmidas_grade HA.MCA1_HAaSx],'Type','Spearman');
title(sprintf('pedmidas grade, Spearmans Rho = %1.2f, p = %1.1e',[rho(1,2) pval(1,2)]))
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XTick = 0:3; ax.XTickLabel = {'none','mild','moderate','severe'}; ylabel('MCA1 factor loading'); ax.XLim = [-1 4];

% bootstrap analysis of headache disability vs. MCA rho value

bootstat = bootstrp(6662,@(a,b)corr([a b],'Type','Spearman'),HA.pedmidas_grade,HA.MCA1_HAaSx);
fprintf('pedmidas vs. MCA1, Rho = %1.2f; 95CI %1.2f - %1.2f  \n',prctile(bootstat(:,2),[50 2.5 97.5]))

bootstat = bootstrp(6662,@(a,b)corr([a b],'Type','Spearman'),HA.pedmidas_grade,HA.MCA2_HAaSx);
fprintf('pedmidas vs. MCA2, Rho = %1.2f; 95CI %1.2f - %1.2f  \n',prctile(bootstat(:,2),[50 2.5 97.5]))

subplot(2,3,3)
hold on
cont = HA.MCA1_HAaSx(HA.p_current_ha_pattern=='cons_flare'|HA.p_current_ha_pattern=='cons_same');
int = HA.MCA1_HAaSx(HA.p_current_ha_pattern=='episodic');
jtr = (rand(size(cont))-0.5).*0.1;
plot(2*ones(size(cont))+jtr,cont,'.k','MarkerSize',8)
errorbar(2,prctile(cont,50),abs(diff(prctile(cont,[25 50]))),abs(diff(prctile(cont,[50 75]))),'+r','MarkerSize',12,'MarkerEdgeColor','r')
jtr = (rand(size(int))-0.5).*0.1;
plot(ones(size(int))+jtr,int,'.k','MarkerSize',8)
errorbar(1,prctile(int,50),abs(diff(prctile(int,[25 50]))),abs(diff(prctile(int,[50 75]))),'+r','MarkerSize',12,'MarkerEdgeColor','r')
cohD = (mean(cont)-mean(int))./std(HA.MCA1_HAaSx);
title(sprintf('continuous vs. MCA1, Cohens D = %1.2f',cohD))
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XTick = 1:2; ax.XTickLabel = {'intermittent','continuous'}; ylabel('MCA1 factor loading'); ax.XLim = [0 3];

% Cohen's D 95CI for continuous HA vs. MCA
bootstat = bootstrp(6662,@(a)[mean(a.MCA1_HAaSx(a.p_current_ha_pattern=='cons_flare'|a.p_current_ha_pattern=='cons_same')) mean(a.MCA1_HAaSx(a.p_current_ha_pattern=='episodic')) std(a.MCA1_HAaSx)],HA(HA.MCA1_HAaSx~=min(HA.MCA1_HAaSx),:));
bootstat2 = (bootstat(:,1)-bootstat(:,2))./bootstat(:,3);
fprintf('continuous vs. MCA1, Cohens D = %1.2f; 95CI %1.2f - %1.2f  \n',prctile(bootstat2,[50 2.5 97.5]))


bootstat = bootstrp(6662,@(a)[mean(a.MCA2_HAaSx(a.p_current_ha_pattern=='cons_flare'|a.p_current_ha_pattern=='cons_same')) mean(a.MCA2_HAaSx(a.p_current_ha_pattern=='episodic')) std(a.MCA2_HAaSx)],HA);
bootstat2 = (bootstat(:,1)-bootstat(:,2))./bootstat(:,3);
fprintf('continuous vs. MCA2, Cohens D = %1.2f; 95CI %1.2f - %1.2f  \n',prctile(bootstat2,[50 2.5 97.5]))


subplot(2,3,4)
hold on
female = HA.MCA1_HAaSx(HA.gender==1);
male = HA.MCA1_HAaSx(HA.gender==2);
jtr = (rand(size(male))-0.5).*0.1;
plot(1*ones(size(male))+jtr,male,'.k','MarkerSize',8)
errorbar(1,prctile(male,50),abs(diff(prctile(male,[25 50]))),abs(diff(prctile(male,[50 75]))),'+r','MarkerSize',12,'MarkerEdgeColor','r')
jtr = (rand(size(female))-0.5).*0.1;
plot(2*ones(size(female))+jtr,female,'.k','MarkerSize',8)
errorbar(2,prctile(female,50),abs(diff(prctile(female,[25 50]))),abs(diff(prctile(female,[50 75]))),'+r','MarkerSize',12,'MarkerEdgeColor','r')
cohD = (mean(female)-mean(male))./std(HA.MCA1_HAaSx);
title(sprintf('male vs. female, Cohens D = %1.2f',cohD))
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XTick = 1:2; ax.XTickLabel = {'male','female'}; ylabel('MCA1 factor loading'); ax.XLim = [0 3];

% Cohen's D for sex assigned at birth vs. MCA
bootstat = bootstrp(6662,@(a)[mean(a.MCA1_HAaSx(a.gender==1)) mean(a.MCA1_HAaSx(a.gender==2)) std(a.MCA1_HAaSx)],HA);
bootstat2 = (bootstat(:,1)-bootstat(:,2))./bootstat(:,3);
fprintf('sex vs. MCA1, Cohens D = %1.2f; 95CI %1.2f - %1.2f  \n',prctile(bootstat2,[50 2.5 97.5]))

bootstat = bootstrp(6662,@(a)[mean(a.MCA2_HAaSx(a.gender==1)) mean(a.MCA2_HAaSx(a.gender==2)) std(a.MCA2_HAaSx)],HA);
bootstat2 = (bootstat(:,1)-bootstat(:,2))./bootstat(:,3);
fprintf('sex vs. MCA2, Cohens D = %1.2f; 95CI %1.2f - %1.2f  \n',prctile(bootstat2,[50 2.5 97.5]))

subplot(2,3,5)
hold on
plot(HA.age,HA.MCA1_HAaSx,'.k','MarkerSize',8)
lsline
[r,pval] = corr([HA.age HA.MCA1_HAaSx]);
title(sprintf('Age, Pearsons R = %1.2f, p = %1.1e',[r(1,2) pval(1,2)]))
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('MCA1 factor loading'); xlabel('Age');

% R 95CI for age vs. MCA
bootstat = bootstrp(6662,@corr,HA.age,HA.MCA1_HAaSx);
fprintf('age vs. MCA1, Pearsons R = %1.2f; 95CI %1.2f - %1.2f  \n',prctile(bootstat,[50 2.5 97.5]))

bootstat = bootstrp(6662,@corr,HA.age,HA.MCA2_HAaSx);
fprintf('age vs. MCA2, Pearsons R = %1.2f; 95CI %1.2f - %1.2f  \n',prctile(bootstat,[50 2.5 97.5]))

M = HA(HA.gender==2,:);
bootstat = bootstrp(height(M),@corr,M.age,M.MCA1_HAaSx);
fprintf('age (male) vs. MCA1, Pearsons R = %1.2f; 95CI %1.2f - %1.2f  \n',prctile(bootstat,[50 2.5 97.5]))

F = HA(HA.gender==1,:);
bootstat = bootstrp(height(F),@corr,F.age,F.MCA1_HAaSx);
fprintf('age (female) vs. MCA1, R = %1.2f; 95CI %1.2f - %1.2f  \n',prctile(bootstat,[50 2.5 97.5]))

subplot(2,3,6)
hold on
none3 = HA.MCA1_HAaSx(HA.freq_bad==1);
min3 = HA.MCA1_HAaSx(HA.freq_bad==2);
mild3 = HA.MCA1_HAaSx(HA.freq_bad==3);
mod3 = HA.MCA1_HAaSx(HA.freq_bad==4);
modsev3 = HA.MCA1_HAaSx(HA.freq_bad==5);
sev3 = HA.MCA1_HAaSx(HA.freq_bad==6);
Vsev3 = HA.MCA1_HAaSx(HA.freq_bad==7);
Esev3 = HA.MCA1_HAaSx(HA.freq_bad==8);
jtr = (rand(size(none3))-0.5).*0.1;
plot(0*ones(size(none3))+jtr,none3,'.k','MarkerSize',8)
errorbar(0,prctile(none3,50),abs(diff(prctile(none3,[25 50]))),abs(diff(prctile(none3,[50 75]))),'+r','MarkerSize',12,'MarkerEdgeColor','r')
jtr = (rand(size(min3))-0.5).*0.1;
plot(ones(size(min3))+jtr,min3,'.k','MarkerSize',8)
errorbar(1,prctile(min3,50),abs(diff(prctile(min3,[25 50]))),abs(diff(prctile(min3,[50 75]))),'+r','MarkerSize',12,'MarkerEdgeColor','r')
jtr = (rand(size(mild3))-0.5).*0.1;
plot(2*ones(size(mild3))+jtr,mild3,'.k','MarkerSize',8)
errorbar(2,prctile(mild3,50),abs(diff(prctile(mild3,[25 50]))),abs(diff(prctile(mild3,[50 75]))),'+r','MarkerSize',12,'MarkerEdgeColor','r')
jtr = (rand(size(mod3))-0.5).*0.1;
plot(3*ones(size(mod3))+jtr,mod3,'.k','MarkerSize',8)
errorbar(3,prctile(mod3,50),abs(diff(prctile(mod3,[25 50]))),abs(diff(prctile(mod3,[50 75]))),'+r','MarkerSize',12,'MarkerEdgeColor','r')
jtr = (rand(size(modsev3))-0.5).*0.1;
plot(4*ones(size(modsev3))+jtr,modsev3,'.k','MarkerSize',8)
errorbar(4,prctile(modsev3,50),abs(diff(prctile(modsev3,[25 50]))),abs(diff(prctile(modsev3,[50 75]))),'+r','MarkerSize',12,'MarkerEdgeColor','r')
jtr = (rand(size(sev3))-0.5).*0.1;
plot(5*ones(size(sev3))+jtr,sev3,'.k','MarkerSize',8)
errorbar(5,prctile(sev3,50),abs(diff(prctile(sev3,[25 50]))),abs(diff(prctile(sev3,[50 75]))),'+r','MarkerSize',12,'MarkerEdgeColor','r')
jtr = (rand(size(Vsev3))-0.5).*0.1;
plot(6*ones(size(Vsev3))+jtr,Vsev3,'.k','MarkerSize',8)
errorbar(6,prctile(Vsev3,50),abs(diff(prctile(Vsev3,[25 50]))),abs(diff(prctile(Vsev3,[50 75]))),'+r','MarkerSize',12,'MarkerEdgeColor','r')
jtr = (rand(size(Esev3))-0.5).*0.1;
plot(7*ones(size(Esev3))+jtr,Esev3,'.k','MarkerSize',8)
errorbar(7,prctile(Esev3,50),abs(diff(prctile(Esev3,[25 50]))),abs(diff(prctile(Esev3,[50 75]))),'+r','MarkerSize',12,'MarkerEdgeColor','r')
[rho,pval] = corr([HA.freq_bad(~isnan(HA.freq_bad)) HA.MCA1_HAaSx(~isnan(HA.freq_bad))],'Type','Spearman');
title(sprintf('pedmidas grade, Spearmans Rho = %1.2f, p = %1.1e',[rho(1,2) pval(1,2)]))
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XTick = 0:7; ax.XTickLabel = {'never','1/month','1-3/month','1/week','2-3/week','>3/week','daily','always'}; ylabel('MCA1 factor loading'); ax.XLim = [-1 8];


bootstat = bootstrp(6662,@(a,b)corr([a b],'Type','Spearman'),HA.freq_bad(~isnan(HA.freq_bad)),HA.MCA1_HAaSx(~isnan(HA.freq_bad)));
fprintf('HA freq vs. MCA1, Rho = %1.2f; 95CI %1.2f - %1.2f  \n',prctile(bootstat(:,2),[50 2.5 97.5]))

bootstat = bootstrp(6662,@(a,b)corr([a b],'Type','Spearman'),HA.freq_bad(~isnan(HA.freq_bad)),HA.MCA2_HAaSx(~isnan(HA.freq_bad)));
fprintf('HA freq vs. MCA2, Rho = %1.2f; 95CI %1.2f - %1.2f  \n',prctile(bootstat(:,2),[50 2.5 97.5]))

% lmFitMCA1 = fitlm(HA,'MCA1_HAaSx ~ severity_grade + pedmidas_grade + gender + age');


% Box plots
figure
subplot(2,3,1)
boxplot(HA.MCA1_HAaSx,HA.severity_grade)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XTick = 1:3; ax.XTickLabel = {'mild','moderate','severe'}; ylabel('MCA1 factor loading'); ax.XLim = [0 4];

subplot(2,3,2)
boxplot(HA.MCA1_HAaSx,HA.pedmidas_grade)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XTick = 1:4; ax.XTickLabel = {'none','mild','moderate','severe'}; ylabel('MCA1 factor loading'); ax.XLim = [0 5];

subplot(2,3,3)
boxplot(HA.MCA1_HAaSx,HA.freq_bad)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XTick = 1:8; ax.XTickLabel = {'never','1/month','1-3/month','1/week','2-3/week','>3/week','daily','always'}; ylabel('MCA1 factor loading'); ax.XLim = [0 9];

subplot(2,3,4)
hold on
hold on
plot(HA.age(HA.gender==1),HA.MCA1_HAaSx(HA.gender==1),'.k','MarkerSize',8)
lsline
plot(HA.age(HA.gender==2),HA.MCA1_HAaSx(HA.gender==2),'.g','MarkerSize',8)
lsline
[r,pval] = corr([HA.age HA.MCA1_HAaSx]);
title(sprintf('Age, Pearsons R = %1.2f, p = %1.1e',[r(1,2) pval(1,2)]))
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('MCA1 factor loading'); xlabel('Age');

subplot(2,3,5:6)
hold on
bins = 6:1:18;
Fmca1 = zeros(length(bins)-1,1);
Mmca1 = zeros(length(bins)-1,1);
for x = 1:length(bins)-1
    Fmca1(x,:) = median(HA.MCA1_HAaSx(HA.gender==1 & HA.age>=bins(x) & HA.age<bins(x+1)));
    Mmca1(x,:) = median(HA.MCA1_HAaSx(HA.gender==2 & HA.age>=bins(x) & HA.age<bins(x+1)));
end
plot(bins(1:end-1)+0.5,Fmca1,'-ob')
plot(bins(1:end-1)+0.5,Mmca1,'-or')
title(sprintf('Age, Pearsons R = %1.2f, p = %1.1e',[r(1,2) pval(1,2)]))
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('MCA1 factor loading'); xlabel('Age');

%% Primary headache diagnosis - migraine vs. probable migraine vs. tension-type headache

figure
p = [5 50 95];
migr_mca1 = prctile(HA.MCA1_HAaSx(HA.p_migraine_ichd==1 & HA.p_dx_overall_cat==1 & HA.p_dx_overall_pheno<4),p);
pmigr_mca1 = prctile(HA.MCA1_HAaSx(HA.p_migraine_ichd==0 & HA.p_dx_overall_cat==1 & HA.p_dx_overall_pheno<4),p);
tth_mca1 = prctile(HA.MCA1_HAaSx(HA.p_dx_overall_cat==1 & (HA.p_dx_overall_pheno==5|HA.p_dx_overall_pheno==6)),p);

migr_mca2 = prctile(HA.MCA2_HAaSx(HA.p_migraine_ichd==1 & HA.p_dx_overall_cat==1 & HA.p_dx_overall_pheno<4),p);
pmigr_mca2 = prctile(HA.MCA2_HAaSx(HA.p_migraine_ichd==0 & HA.p_dx_overall_cat==1 & HA.p_dx_overall_pheno<4),p);
tth_mca2 = prctile(HA.MCA2_HAaSx(HA.p_dx_overall_cat==1 & (HA.p_dx_overall_pheno==5|HA.p_dx_overall_pheno==6)),p);

errorbar(migr_mca1(2),migr_mca2(2),abs(diff(migr_mca2(1:2))),abs(diff(migr_mca2(2:3))),abs(diff(migr_mca1(1:2))),abs(diff(migr_mca1(2:3))),'ok')
hold on
errorbar(pmigr_mca1(2),pmigr_mca2(2),abs(diff(pmigr_mca2(1:2))),abs(diff(pmigr_mca2(2:3))),abs(diff(pmigr_mca1(1:2))),abs(diff(pmigr_mca1(2:3))),'ob')
errorbar(tth_mca1(2),tth_mca2(2),abs(diff(tth_mca2(1:2))),abs(diff(tth_mca2(2:3))),abs(diff(tth_mca1(1:2))),abs(diff(tth_mca1(2:3))),'oc')
plot([-6 6],[0 0],'--c')
plot([0 0],[-4 4],'--c')
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; title('MCA: migraine vs. probable migraine vs. tension-type headache')


migr_mca1 = HA.MCA1_HAaSx(HA.p_migraine_ichd==1 & HA.p_dx_overall_cat==1 & HA.p_dx_overall_pheno<4);
pmigr_mca1 = HA.MCA1_HAaSx(HA.p_migraine_ichd==0 & HA.p_dx_overall_cat==1 & HA.p_dx_overall_pheno<4);
tth_mca1 = HA.MCA1_HAaSx(HA.p_dx_overall_cat==1 & (HA.p_dx_overall_pheno==5|HA.p_dx_overall_pheno==6));
pth_mca1 = HA.MCA1_HAaSx(HA.p_dx_overall_cat==6);
ndph_mca1 = HA.MCA1_HAaSx(HA.p_dx_overall_cat==2 | HA.p_dx_overall_cat==3);
tac_mca1 = HA.MCA1_HAaSx(HA.p_dx_overall_cat==1 & HA.p_dx_overall_pheno==4);

migr_mca2 = HA.MCA2_HAaSx(HA.p_migraine_ichd==1 & HA.p_dx_overall_cat==1 & HA.p_dx_overall_pheno<4);
pmigr_mca2 = HA.MCA2_HAaSx(HA.p_migraine_ichd==0 & HA.p_dx_overall_cat==1 & HA.p_dx_overall_pheno<4);
tth_mca2 = HA.MCA2_HAaSx(HA.p_dx_overall_cat==1 & (HA.p_dx_overall_pheno==5|HA.p_dx_overall_pheno==6));
pth_mca2 = HA.MCA2_HAaSx(HA.p_dx_overall_cat==6);
ndph_mca2 = HA.MCA2_HAaSx(HA.p_dx_overall_cat==2 | HA.p_dx_overall_cat==3);
tac_mca2 = HA.MCA2_HAaSx(HA.p_dx_overall_cat==1 & HA.p_dx_overall_pheno==4);

figure
subplot(2,1,1)
plot(migr_mca1+rand(size(migr_mca1))*0.05,migr_mca2+rand(size(migr_mca2))*0.05,'.k','MarkerSize',12)
hold on
plot(pmigr_mca1+rand(size(pmigr_mca1))*0.05,pmigr_mca2+rand(size(pmigr_mca2))*0.05,'.b','MarkerSize',12)
plot(tth_mca1+rand(size(tth_mca1))*0.05,tth_mca2+rand(size(tth_mca2))*0.05,'.c','MarkerSize',12)
plot(tac_mca1+rand(size(tac_mca1))*0.05,tac_mca2+rand(size(tac_mca2))*0.05,'.y','MarkerSize',12)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; title('MCA: primary headache types')

subplot(2,1,2)
plot(pth_mca1+rand(size(pth_mca1))*0.05,pth_mca2+rand(size(pth_mca2))*0.05,'.r','MarkerSize',12)
hold on
plot(ndph_mca1+rand(size(ndph_mca1))*0.05,ndph_mca2+rand(size(ndph_mca2))*0.05,'.m','MarkerSize',12)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; title('MCA: other headache types')


% Compare ICHD3 diagnosis
ichd_mca1 = [migr_mca1;pmigr_mca1;tth_mca1;pth_mca1;ndph_mca1];
ichd_mca2 = [migr_mca2;pmigr_mca2;tth_mca2;pth_mca2;ndph_mca2];
ichd = [repmat({'mig'},length(migr_mca1),1);repmat({'pmg'},length(pmigr_mca1),1);repmat({'tth'},length(tth_mca1),1);repmat({'pth'},length(pth_mca1),1);repmat({'ndp'},length(ndph_mca1),1)];

[d,p,stats] = manova1([ichd_mca1 ichd_mca2],ichd);
gplotmatrix(ichd_mca1,ichd_mca2,ichd)

% ridgeline plot of MCA1 and MCA2 by ICHD diagnosis
figure
subplot(6,2,1)
pd = fitdist(tth_mca1,'Kernel','Bandwidth',0.6);
y = pdf(pd,x);
plot(x,y)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('TTH, MCA1'); ax.YLim = [0 0.5]; ax.XLim = [-6 10];

subplot(6,2,3)
pd = fitdist(pmigr_mca1,'Kernel','Bandwidth',0.6); 
y = pdf(pd,x);
plot(x,y)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('probable migraine, MCA1'); ax.YLim = [0 0.5]; ax.XLim = [-6 10];

subplot(6,2,5)
pd = fitdist(migr_mca1,'Kernel','Bandwidth',0.6);
x = -6:0.1:11.5;
y = pdf(pd,x);
plot(x,y)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('migraine, MCA1'); ax.YLim = [0 0.5]; ax.XLim = [-6 10];

subplot(6,2,7)
pd = fitdist(pth_mca1,'Kernel','Bandwidth',0.6);
y = pdf(pd,x);
plot(x,y)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('PTH, MCA1'); ax.YLim = [0 0.5]; ax.XLim = [-6 10];

subplot(6,2,9)
pd = fitdist(ndph_mca1,'Kernel','Bandwidth',0.6);
y = pdf(pd,x);
plot(x,y)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('NDPH/new onset, MCA1'); ax.YLim = [0 0.5]; ax.XLim = [-6 10];

subplot(6,2,11)
pd = fitdist(tac_mca1,'Kernel','Bandwidth',0.6);
y = pdf(pd,x);
plot(x,y)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('TAC, MCA1'); ax.YLim = [0 0.5]; ax.XLim = [-6 10];



subplot(6,2,2)
pd1 = fitdist(tth_mca2,'Kernel','Bandwidth',0.25);
pd2 = fitdist(tth_mca2(tth_mca1~=min(HA.MCA1_HAaSx)),'Kernel','Bandwidth',0.25);
x = -5:0.1:5;
y = pdf(pd1,x);
y2 = pdf(pd2,x);
y2 = y2.*(length(tth_mca2(tth_mca1~=min(HA.MCA1_HAaSx)))./length(tth_mca2)); % normalize pdf2 based on the total proportional number to pdf1
plot(x,y,'--')
hold on
plot(x,y2)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('TTH, MCA2'); ax.YLim = [0 1.3]; ax.XLim = [-5 5];

subplot(6,2,4)
pd1 = fitdist(pmigr_mca2,'Kernel','Bandwidth',0.25);
pd2 = fitdist(pmigr_mca2(pmigr_mca1~=min(HA.MCA1_HAaSx)),'Kernel','Bandwidth',0.25);
x = -5:0.1:5;
y = pdf(pd1,x);
y2 = pdf(pd2,x);
y2 = y2.*(length(pmigr_mca2(pmigr_mca1~=min(HA.MCA1_HAaSx)))./length(pmigr_mca2)); % normalize pdf2 based on the total proportional number to pdf1
plot(x,y,'--')
hold on
plot(x,y2)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('probable migraine, MCA2'); ax.YLim = [0 1.3]; ax.XLim = [-5 5];

subplot(6,2,6)
pd1 = fitdist(migr_mca2,'Kernel','Bandwidth',0.25);
pd2 = fitdist(migr_mca2(migr_mca1~=min(HA.MCA1_HAaSx)),'Kernel','Bandwidth',0.25);
x = -5:0.1:5;
y = pdf(pd1,x);
y2 = pdf(pd2,x);
y2 = y2.*(length(migr_mca2(migr_mca1~=min(HA.MCA1_HAaSx)))./length(migr_mca2)); % normalize pdf2 based on the total proportional number to pdf1
plot(x,y,'--')
hold on
plot(x,y2)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('Migraine, MCA2'); ax.YLim = [0 1.3]; ax.XLim = [-5 5];

subplot(6,2,8)
pd1 = fitdist(pth_mca2,'Kernel','Bandwidth',0.25);
pd2 = fitdist(pth_mca2(pth_mca1~=min(HA.MCA1_HAaSx)),'Kernel','Bandwidth',0.25);
x = -5:0.1:5;
y = pdf(pd1,x);
y2 = pdf(pd2,x);
y2 = y2.*(length(pth_mca2(pth_mca1~=min(HA.MCA1_HAaSx)))./length(pth_mca2)); % normalize pdf2 based on the total proportional number to pdf1
plot(x,y,'--')
hold on
plot(x,y2)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('PTH, MCA2'); ax.YLim = [0 1.3]; ax.XLim = [-5 5];

subplot(6,2,10)
pd1 = fitdist(ndph_mca2,'Kernel','Bandwidth',0.25);
pd2 = fitdist(ndph_mca2(ndph_mca1~=min(HA.MCA1_HAaSx)),'Kernel','Bandwidth',0.25);
x = -5:0.1:5;
y = pdf(pd1,x);
y2 = pdf(pd2,x);
y2 = y2.*(length(ndph_mca2(ndph_mca1~=min(HA.MCA1_HAaSx)))./length(ndph_mca2)); % normalize pdf2 based on the total proportional number to pdf1
plot(x,y,'--')
hold on
plot(x,y2)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('NDPH/new onset, MCA2'); ax.YLim = [0 1.3]; ax.XLim = [-5 5];

subplot(6,2,12)
pd1 = fitdist(tac_mca2,'Kernel','Bandwidth',0.25);
pd2 = fitdist(tac_mca2(tac_mca1~=min(HA.MCA1_HAaSx)),'Kernel','Bandwidth',0.25);
x = -5:0.1:5;
y = pdf(pd1,x);
y2 = pdf(pd2,x);
y2 = y2.*(length(tac_mca2(tac_mca1~=min(HA.MCA1_HAaSx)))./length(tac_mca2)); % normalize pdf2 based on the total proportional number to pdf1
plot(x,y,'--')
hold on
plot(x,y2)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('TAC, MCA2'); ax.YLim = [0 1.3]; ax.XLim = [-5 5];





% histogram plot of MCA1 and MCA2 by ICHD diagnosis
figure
subplot(6,1,1)
x = -5:0.25:12;
histogram(migr_mca1,x)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('migraine, MCA1');

subplot(6,1,2)
histogram(pmigr_mca1,x)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('probable migraine, MCA1');

subplot(6,1,3)
histogram(tth_mca1,x)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('TTH, MCA1');

subplot(6,1,4)
histogram(pth_mca1,x)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('PTH, MCA1');

subplot(6,1,5)

histogram(ndph_mca1,x)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('NDPH/new onset, MCA1');

subplot(6,1,6)
histogram(tac_mca1,x)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('TAC, MCA1');




figure
subplot(6,1,1)
x = -5:0.25:5;
histogram(migr_mca2,x)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('migraine, MCA2');

subplot(6,1,2)
histogram(pmigr_mca2,x)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('probable migraine, MCA2');

subplot(6,1,3)
histogram(tth_mca2,x)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('TTH, MCA2');

subplot(6,1,4)
histogram(pth_mca2,x)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('PTH, MCA2');

subplot(6,1,5)
histogram(ndph_mca2,x)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('NDPH/new onset, MCA2');

subplot(6,1,6)
histogram(tac_mca2,x)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ylabel('TAC, MCA2');




% Histogram on total number of symptoms reported by ICHD
sxN = sum(haASx,2);
migr_sxN = sxN(HA.p_migraine_ichd==1 & HA.p_dx_overall_cat==1 & HA.p_dx_overall_pheno<4);
pmigr_sxN = sxN(HA.p_migraine_ichd==0 & HA.p_dx_overall_cat==1 & HA.p_dx_overall_pheno<4);
tth_sxN = sxN(HA.p_dx_overall_cat==1 & (HA.p_dx_overall_pheno==5|HA.p_dx_overall_pheno==6));
pth_sxN = sxN(HA.p_dx_overall_cat==6);
ndph_sxN = sxN(HA.p_dx_overall_cat==2 | HA.p_dx_overall_cat==3);
tac_sxN = sxN(HA.p_dx_overall_cat==1 & HA.p_dx_overall_pheno==4);

figure
bins = 0:1:13;
subplot(6,1,1)
histogram(migr_sxN,bins)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off';
title('migraine')

subplot(6,1,2)
histogram(pmigr_sxN,bins)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off';
title('probable migraine')

subplot(6,1,3)
histogram(tth_sxN,bins)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off';
title('TTH')

subplot(6,1,4)
histogram(pth_sxN,bins)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off';
title('PTH')

subplot(6,1,5)
histogram(ndph_sxN,bins)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off';
title('NDPH')

subplot(6,1,6)
histogram(tac_sxN,bins)
ax = gca; ax.TickDir = 'out'; ax.Box = 'off';
title('TAC')

% Determine what percentage of probable migraine would be included in
% migraine diagnosis if dizziness was added as a symptom
HA.dizzy = sum([HA.p_assoc_sx_oth_sx___spinning HA.p_assoc_sx_oth_sx___balance==1 HA.p_assoc_sx_oth_sx___lighthead==0],2);
pmigr = HA(HA.p_migraine_ichd==0 & HA.p_dx_overall_cat==1 & HA.p_dx_overall_pheno<4,:);
pmigr_assoc = pmigr(isnan(pmigr.p_mig_d_nausea_vomit) & isnan(pmigr.p_mig_d_photo_phono) & pmigr.p_mig_b==1 & pmigr.p_mig_c==1,:);
pmigr_assoc_dizzy = pmigr_assoc(pmigr_assoc.dizzy>0,:);
pmigr_assoc_lighthead = pmigr_assoc(pmigr_assoc.p_assoc_sx_oth_sx___lighthead==1 & pmigr_assoc.p_assoc_sx_oth_sx___balance==0 & pmigr_assoc.p_assoc_sx_oth_sx___spinning==0,:);
pmigr_assoc_vertigo = pmigr_assoc((pmigr_assoc.p_assoc_sx_oth_sx___spinning==1 | pmigr_assoc.p_assoc_sx_oth_sx___balance==1) & pmigr_assoc.p_assoc_sx_oth_sx___lighthead==0,:);
pmigr_assoc_lhv = pmigr_assoc((pmigr_assoc.p_assoc_sx_oth_sx___spinning==1 | pmigr_assoc.p_assoc_sx_oth_sx___balance==1) & pmigr_assoc.p_assoc_sx_oth_sx___lighthead==1,:);
pmigr_assoc_notdizzy = pmigr_assoc(pmigr_assoc.dizzy==0,:);
