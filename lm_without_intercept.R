# This is a simple script showcasing how to compute a standard linear model without an intercept.
# Which, bizarrely enough, is not very intuitive.

# xvals is vector containing the regressor, yvals vector for regressand

xvals = c(3,4,2)
yvals = c(5,7,3)

# just standard lm for sake of comparison
lmwithintercept <- lm(yvals ~ xvals)

# the intercept to use is set by offset, making a vector containing only 0 means lm function will go through origin
lmwithoutintercept <- lm(yvals~ 0 + xvals, offset=rep(0,length(xvals)))