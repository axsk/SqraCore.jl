using SqraCore

sqra(rand(10))
sqra(rand(5,5,5), beta=3)

sqra(x->x^2, rand(5))
sqra(x->x^2, rand(5,5))
