# clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()

VaRest = function(y, method, alpha) {
  # parameter settings
  n = length(y)
  h = 240
  lam = 0.96
  w = 1
  
  # RMA
  if (method == 1) {
    sigma = matrix(1, (n - h), 1) - 1
    tmp = cumsum(y * y)
    tmp1 = (tmp[(h + 1):n] - tmp[1:(n - h)])/h
    sigma = sqrt(((w * tmp1) * w))
  }
  
  # EMA
  if (method == 2) {
    sigma = matrix(1, (n - h), 1) - 1
    j = h
    expo = (h - 1):0
    while (j < n) {
      j = j + 1
      tmp = (lam^expo)^0.5 * y[(j - h):(j - 1)]
      tmp1 = sum(tmp * tmp)
      sigma[j - h] = sqrt(tmp1 * (1 - lam))
    }
  }
  
  # Selecting alpha
  if (alpha == 0.01) {
    qf = qnorm(alpha, 0, 1)
    VaR = qf * sigma
    VaR = cbind(VaR, (-VaR))
  } else if (alpha == 0.05) {
    qf = qnorm(alpha, 0, 1)
    VaR = qf * sigma
    VaR = cbind(VaR, (-VaR))
  }
  
}

# Main computation
x1 = read.table("CMS.txt")
x = x1[1:nrow(x1), 1]
y = diff(log(x))
h = 240

# Option 1_1=RMA of 1%, Option 2_1=EMA of 1%, Option 1_2=RMA of 5%, Option 2_2=RMA of
# 5%
opt1_1 = VaRest(y, 1, 0.01)
opt2_1 = VaRest(y, 2, 0.01)
opt1_2 = VaRest(y, 1, 0.05)
opt2_2 = VaRest(y, 2, 0.05)

# Plots

# Subplots
A = c(1, 2)
layout(A)

# 1.Plots for estimated VaR of 1%
time = seq(h + 1, length(x) - 1)
plot(time, opt1_1[, 1], col = "green", type = "l", lty = "longdash", ylim = c(0.15, -0.15), 
     main = "VaR TimePlot (alpha=1%)", xlab = "Time", ylab = "Returns")
lines(time, opt1_1[, 2], col = "green", lty = "longdash")
lines(time, opt2_1[, 1], col = "blue", lty = "solid")
lines(time, opt2_1[, 2], col = "blue", lty = "solid")

points(241:(length(x) - 1), y[241:(length(x) - 1)], col = "black", pch = 4)

exceed = matrix(0, length(241:length(y)))
l = 1
for (i in 241:length(y)) {
  if ((opt1_1[i - 240, 2] < y[i]) || (y[i] < opt1_1[i - 240, 1])) {
    exceed[l] = y[i]
    points(i, exceed[l], col = "red", pch = 12)
    l = l + 1
  }
}

# 2.Plots for estimated VaR of 5%
time = seq(h + 1, length(x) - 1)
plot(time, opt1_2[, 1], col = "green", type = "l", lty = "longdash", ylim = c(0.15, -0.15), 
     main = "VaR TimePlot (alpha=5%)", xlab = "Time", ylab = "Returns")
lines(time, opt1_2[, 2], col = "green", lty = "longdash")
lines(time, opt2_2[, 1], col = "blue", lty = "solid")
lines(time, opt2_2[, 2], col = "blue", lty = "solid")

points(241:(length(x) - 1), y[241:(length(x) - 1)], col = "black", pch = 4)

exceed = matrix(0, length(241:length(y)))
l = 1
for (i in 241:length(y)) {
  if ((opt1_2[i - 240, 2] < y[i]) || (y[i] < opt1_2[i - 240, 1])) {
    exceed[l] = y[i]
    points(i, exceed[l], col = "red", pch = 12)
    l = l + 1
  }
}