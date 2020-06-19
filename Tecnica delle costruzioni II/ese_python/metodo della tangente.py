def quadrato(x):
	return x*x

def d_quadrato(x):
	return 2*x

def tangenti(x_k, scarto, f, df):
    #controllo sulla vicinanza del risultato allo zero
	while f(x_k) > scarto:
		x_k = x_k - ( f(x_k) / df(x_k) )
		tangenti(x_k, scarto, f, df)
	return x_k

init_val = 10.0
scarto = 0.01

x_k = tangenti(init_val, scarto, quadrato, d_quadrato)
print(x_k)
print(quadrato(x_k))
print("E' la soluzione corretta?"),
ans = quadrato(x_k) <= scarto
print(ans)