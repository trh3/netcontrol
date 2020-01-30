###Scratch space

X = c(1:5)
u = -5:-1

u = c(5, -5,5,-5,5)

(X-u)

u.bounds = c(10,10,10,10,10)
l.bounds = c(1,1,1,1,1)

u.lb = pmin(u.bounds-X,u)
 pmax(l.bounds-X, pmin(u.bounds-X,u))