# check that 0.1 + 0.2 = 0.3
# this will fail because of floating point precision issues
# add(0.1, 0.2) == 0.3
isapprox(add(0.1, 0.2), 0.3) # checks if in floating point precision