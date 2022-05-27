tmp = data.frame(coords.x = sim2.1c1$data$x.coord,
                 coords.y = sim2.1c1$data$y.coord,
                 land = sim2.1c1$data$land,
                 X1 = sim2.1c1$data$X.1,
                 X2 = sim2.1c2$data$X.1,
                 X3 = sim2.1c3$data$X.1,
                 eps1 = sim2.1c1$eps,
                 eps2 = sim2.1c2$eps,
                 eps3 = sim2.1c3$eps,
                 r1 = sim2.1c1$r1,
                 r2 = sim2.1c2$r1,
                 r3 = sim2.1c3$r1)
# check data
tmp %>%
  reshape2::melt(measure.vars = c("X1", "X2", "X3")) %>%
  ggplot(aes(x = coords.x, y = coords.y, fill = value, col = land)) +
  geom_tile() +
  facet_wrap(~variable, ncol = 1) +
  scale_fill_gradient2(low = "blue", high = "orange") +
  scale_color_manual(values = c("grey70", "grey20")) +
  labs(fill = "X")

# check varepsilon
tmp %>%
  reshape2::melt(measure.vars = c("eps1", "eps2", "eps3")) %>%
  ggplot(aes(x = coords.x, y = coords.y, fill = value)) +
  geom_tile() +
  facet_wrap(~variable, ncol = 1) +
  scale_fill_gradient2(low = "blue", high = "orange") +
  scale_color_manual(values = c("grey70", "grey20")) +
  labs(fill = "varepsilon")

tmp %>%
  reshape2::melt(measure.vars = c("r1", "r2", "r3")) %>%
  ggplot(aes(x = coords.x, y = coords.y, fill = value, col = land)) +
  geom_tile() +
  facet_wrap(~variable, ncol = 1) +
  scale_fill_gradient2(low = "blue", high = "orange") +
  scale_color_manual(values = c("grey70", "grey20")) +
  labs(fill = "R1")

# check varepsilon + R1
S1 = sim2.1c1$sim.table$S1[1]
S1 = .5
tmp %>%
  mutate(err1 = eps1 + S1*r1, err2 = eps2 + S1*r2, err3 = eps3 + S1*r3) %>%
  reshape2::melt(measure.vars = c("err1", "err2", "err3")) %>%
  ggplot(aes(x = coords.x, y = coords.y, fill = value)) +
  geom_tile() +
  facet_wrap(~variable, ncol = 1) +
  scale_fill_gradient2(low = "blue", high = "orange") +
  scale_color_manual(values = c("grey70", "grey20")) +
  labs(fill = "eF + eR")

# # check delta
# tmp %>%
#   reshape2::melt(measure.vars = c("del1", "del2", "del3")) %>%
#   ggplot(aes(x = coords.x, y = coords.y, fill = value)) +
#   geom_tile() +
#   facet_wrap(~variable, ncol = 1) +
#   labs(fill = "delta")


{
  print("Sim 1:")
  print(lapply(sim2.1c1$GLS, function(x)x$t.pval))
  print("Sim 2:")
  print(lapply(sim2.1c2$GLS, function(x)x$t.pval))
  print("Sim 3:")
  print(lapply(sim2.1c3$GLS, function(x)x$t.pval))
}


sin_2d = function(x, y, A = 1, Tx = 1/(n/2), Ty = 1/(n/2), px = 45, py = 45){
  A * sin( (2*pi/Tx)*x + px + (2*pi/Ty)*y + py )
}

sin_2d_TI= function(x, y, A = 1, Tx = 1/(n/2), Ty = 1/(n/2), px = 45, py = 45){
  A * (sin((2*pi/Tx)*x + px) + sin((2*pi/Ty)*y + py ))
}

tmp$T1 = sin_2d_TI(x = tmp$coords.x, y = tmp$coords.y, Tx = 1/(1/2), Ty = 1/(1/2))
tmp$T2 = sin_2d_TI(x = tmp$coords.x, y = tmp$coords.y, Tx = 1/(8/2), Ty = 1/(8/2))
tmp$T3 = sin_2d_TI(x = tmp$coords.x, y = tmp$coords.y, Tx = 1/(4/2), Ty = 1/(4/2))

tmp %>%
  reshape2::melt(measure.vars = c("T1", "T2", "T3")) %>%
  ggplot(aes(x = coords.x, y = coords.y, fill = value, col = land)) +
  geom_tile() +
  facet_wrap(~variable, ncol = 1) +
  scale_fill_gradient2(low = "blue", high = "orange") +
  scale_color_manual(values = c("grey70", "grey20")) +
  labs(fill = "T1")

# n <- 8
# x <- matrix(rep(1:50, each=50), 50, 50)/50
# y <- t(x)
#
# par(mfrow=c(1,2))
#
# image(sin_2d(x, y))
#
# image(sin_2d_TI(x, y))
