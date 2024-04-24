# Random slope and intercept simulations

set.seed(1)
# Genotype variation in the control (n = 50)
rel.c <- rnorm(50, 1, sd = 0.2)
# Genetic variation of tolerance, i.e. slopes (var = sd^2)
# The population effect is set to zero for simplicity
gen.var <- rnorm(50, 0, sd = 0.30)
# Genotype means for damage treatment
#rel.d <- rel.c + sort(gen.var)
rel.d <- rel.c + gen.var
# Set within genotype:treatment combination sample size
n.in = 20
# Generate within genotype:treatment data with sigma = 0.3
dat.con <- rnorm(50*n.in, rel.c, sd = 0.3)
dat.treat <- rnorm(50*n.in, rel.d, sd = 0.3)
# Give genotype id
gen.id <- rep(1:50, times = n.in)
# Put together data frame
dat <- data.frame( "genotype" = factor(c(gen.id, gen.id)), "density" = rep(c("con","treat"), each = 50*n.in),
                   "jday" = c(dat.con, dat.treat) )

unique(dat$treat)
table(dummy(dat$treat))

summary(mod.lmer)

summary(mod.lmer2)

mod.nest <- lme(jday ~ site, random=~1|genotype/site, data = phen_flower_kin) 
# In lme4
mod.nest.lmer <- lmer(jday ~ site + (1 | genotype)+ (1 | genotype:site), data = phen_flower_kin)

plot(predict(mod.nest),predict(mod.nest.lmer))

# This model can also be defined as estimating a variance for each group (control, treat) instead of estimating the differences as in the model above
mod.2var <- lme(jday ~ density, random=list(genotype = pdDiag(~ density-1)), data = phen_flower_kin) 
hist(ranef(mod.2var)[,2] - ranef(mod.2var)[,1])

# Only random intercept model
mod.null <- lme(rel.fit ~ treat, random = ~ 1|gen, data = dat)

set.seed(1)
# Genotype variation in the control (n = 50)
rel.c <- rnorm(50, 1, sd = 0.2)
# Genetic variation of tolerance, i.e. slopes (var = sd^2)
# The population effect is set to zero for simplicity
gen.var <- rnorm(50, 0, sd = 0.30)
# Genotype means for damage treatment
#rel.d <- rel.c + sort(gen.var)
rel.d <- rel.c + gen.var
# Set within genotype:treatment combination sample size
n.in = 20
# Generate within genotype:treatment data with sigma = 0.3
dat.con <- rnorm(50*n.in, rel.c, sd = 0.3)
dat.treat <- rnorm(50*n.in, rel.d, sd = 0.3)
# Give genotype id
gen.id <- rep(1:50, times = n.in)
# Put together data frame
dat <- data.frame( "gen" = factor(c(gen.id, gen.id)), "treat" = rep(c("con","treat"), each = 50*n.in),
                   "rel.fit" = c(dat.con, dat.treat) )

library(lme4)

anova(mod, mod2)
plot(predict(mod), predict(mod2))

# Genotype variation in the control (n = 50)
rel.c <- rnorm(50, 1, sd = 1)
# Genetic variation of tolerance, i.e. slopes (var = sd^2)
# The population effect is set to zero for simplicity
gen.var <- rnorm(50, 0, sd = 0.3)
# Genotype means for damage treatment
#rel.d <- rel.c + sort(gen.var)
rel.d <- rel.c + gen.var
# Set within genotype:treatment combination sample size
n.in = 20
# Generate within genotype:treatment data with sigma = 0.3
dat.con <- rnorm(50*n.in, rel.c, sd = 0.3)
dat.treat <- rnorm(50*n.in, rel.d, sd = 0.3)
# Give genotype id
gen.id <- rep(1:50, times = n.in)
# Put together data frame
dat <- data.frame( "gen" = factor(c(gen.id, gen.id)), "treat" = rep(c("con","treat"), each = 50*n.in),
                   "response" = c(dat.con, dat.treat) )

dat_sub <- dat[sample(1:2000, 500, replace=F),]

dat_sub %>% 
  group_by(gen, treat) %>% 
  summarize(n=n())

mod <- lmer(response ~ treat + (1|gen) + (1|gen:treat), data = dat_sub)
mod2 <- lmer(response ~ treat + (1|gen:treat), data = dat_sub)

sjPlot::plot_model(mod, type = "pred", pred.type = "re", terms = c("treat", "gen"), ci.lvl = NA, colors = colors) + 
  geom_line() + theme(legend.position = "none")+ylim(-2,4)+ggtitle("Model predictions (sd = 0.3)") -> model

dat_sub %>% 
  group_by(gen, treat) %>% 
  summarize(mean = mean(response)) %>% 
  ggplot(aes(x = treat, y = mean, group = gen, color = gen)) +
  geom_point() + geom_line() +
  scale_color_manual(values = colors) + theme(legend.position = "none") +
  ylim(-2,4) + ggtitle("Raw data (sd = 0.3)") + ylab("response") +
  xlab("treatment")-> data


data + model -> plot_sd2

png("~/Desktop/model_simulations.png", height = 12, width = 8, res = 300, units = "in")
plot_sd1 / plot_sd2  / plot_sd3
dev.off()

phen_flower_kin %>% 
  group_by(genotype) %>% 
  summarize(mean = mean(jday)) %>% 
  pull(mean) %>% sd()

phen_flower_kin %>% 
  group_by(genotype, density) %>% 
  summarize(mean = mean(jday)) %>% 
  spread(density, mean) %>% 
  mutate(diff = lo - hi) %>% 
  pull(diff) %>% sd()

color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
colors <- sample(color, 92)
sjPlot::plot_model(mod, type = "pred", pred.type = "re", terms = c("treat", "gen"), ci.lvl = NA, colors = colors) + 
  geom_line() + theme(legend.position = "none")+ylim(-2,4)+ggtitle("Random slopes and intercept (sd = 0.3)") -> e
sjPlot::plot_model(mod2, type = "pred", pred.type = "re", terms = c("treat", "gen"), ci.lvl = NA, colors = colors) +
  geom_line() + theme(legend.position = "none")+ylim(-2,4)+ggtitle("Random slopes only (sd = 0.3)") -> f

png("~/Desktop/simulations.png", height = 12, width = 10, res = 300, units = "in")
(a+b)/(e+f)/(c+d)
dev.off()

BIC(mod2)

anova(mod, mod2)
plot(predict(mod), predict(mod2))

phen_flower_kin %>% 
  group_by(genotype) %>% 
  summarize(n = n()) %>% 
  ggplot(aes(x = n)) +
  geom_histogram()



