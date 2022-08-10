import random
from decimal import *

getcontext().prec = 500

def predict(s, x, n):
    s = Decimal(s)
    x = Decimal(x)
    n = Decimal(n)
    f = ( 1+s-x ) * (s**-n)
    return (x**n) * f - 2*(x-1)**n * f + max(0, x-2)**n * f


def sample(s, x, n, tries):
    n_success = 0
    for _ in range(tries):
        x_max = 0
        x_min = s + 1
        for _ in range(n):
            x_curr = random.randrange(s)
            x_max = max(x_max, x_curr)
            x_min = min(x_min, x_curr)
        if x == 1 + x_max - x_min:
            n_success += 1
    return Decimal(n_success) / Decimal(tries)



def main(N=10, max_s=1000, max_n_diff=100, round_to=7, num_steps=10):
    for _ in range(N):
        s = random.randrange(max_s) + 1
        n = s + max_n_diff//2 - random.randrange(max_n_diff//2)
        print("")
        print("s", "n", "x", "pred", "sampled", "error", sep="\t")
        for x in range(1, s+1, max(1,s//num_steps)):
            t = 4
            error = 1
            p = predict(s, x, n)
            while abs(error) >= 1e-5 and t < 16:
                sm = sample(s, x, n, 2**t)
                t += 1
                if t >= 18 and abs(p - sm) >= abs(error):
                    print("warning: error =", error, "for t = ", t, "new error", p - sm)
                error = p - sm
            print(s, n, x, round(p, round_to), round(sm, round_to), round( p-sm, round_to*round_to), sep="\t")

main()