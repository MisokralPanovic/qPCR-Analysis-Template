.libPaths('C:/r_packages')

#########
# change 'a549_1',
#########

# data analysis ----
boxplot(Value_norm~Condition, a549_1)
plot(lm(Value_norm~Condition, a549_1))

#### test normality -----------------
shapiro.test(a549_1$Value_norm[1:3])
shapiro.test(a549_1$Value_norm[4:6])
shapiro.test(a549_1$Value_norm[7:9])
shapiro.test(a549_1$Value_norm[10:12])

plot(residuals(lm(Value_norm~Condition, a549_1)))
shapiro.test(residuals(lm(Value_norm~Condition, a549_1)))

#### equality of variance --------------

# if distribution normal
bartlett.test(Value_norm~Condition, a549_1)

# if non normal distribution
library(car)
leveneTest(Value_norm~Condition, a549_1)

#### getting p values ----------------

# for normal distribution and variance
library(DTK)
anova(lm(Value_norm~Condition, a549_1))
TukeyHSD(aov(Value_norm~Condition, a549_1))

# for non normal distribution but equal variance - multiple comparison
library(dunn.test)
kruskal.test(Value_norm~Condition, a549_1)
dunn.test(a549_1$Value_norm, 
          a549_1$Condition, 
          altp=T)

# for non normal distribution but equal variance - single comparison
t.test(Value_norm~Condition, 
       data= a549_1, 
       alternative='two.sided',
       var.equal=T)

# for normal distribution but non equal variance - multiple comparison
library(userfriendlyscience)
oneway.test(Value_norm~Condition, 
            a549_1, 
            var.equal = F)
oneway(a549_1$Value_norm, 
       as.factor(a549_1$Condition), 
       posthoc = 'games-howell')

# for normal distribution but non equal variance single comparison
t.test(Value_norm~Condition, 
       data=a549_1, 
       alternative='two.sided',
       var.equal=F)
