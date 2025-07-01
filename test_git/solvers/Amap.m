function AV = Amap(X,V)

XtV = X'*V;

AV = XtV + XtV';
