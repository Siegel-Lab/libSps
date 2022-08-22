import random
import numpy as np

def predict(num_total, num_distinct, n):
    result = 0
    for _ in range(n):
        num_draws = 0
        for i in range(num_distinct):
            draws_required = np.random.geometric(p=(num_total - i) / num_total)
            #print("xxx", num_total, num_distinct, draws_required, i)
            num_draws += draws_required
        if num_distinct < num_total:
            result += num_draws + (np.random.geometric(p=(num_total - num_distinct) / num_total) - 1 )/ 2
        else:
            result += num_draws
    return result / n


def sample(num_total, num_distinct, n):
    result = 0
    for _ in range(n):
        def f():
            num_draws_now = 0
            num_draws_start = 1
            val_hit = set()

            while len(val_hit) <= num_distinct:
                if len(val_hit) >= num_total:
                    return num_draws_now
                curr_draw = random.randrange(num_total)
                val_hit.add(curr_draw)
                num_draws_now += 1
                if len(val_hit) < num_distinct:
                    num_draws_start = num_draws_now + 1
            return ( num_draws_now + num_draws_start ) // 2
        result += f()
    return result / n




def main(N=10, max_s=1000, max_n_diff=100, round_to=2, num_steps=10):
    for _ in range(N):
        s = random.randrange(max_s) + 1
        print("")
        print("s", "x", "pred", "sampled", "error", sep="\t")
        for x in range(1, s+1, max(1,s//num_steps)):
            t = 4
            error = 100
            while abs(error) >= 10 and t < 18:
                p = predict(s, x, 2**t + 100)
                sm = sample(s, x, 2**t + 100)
                t += 1
                if t >= 15 and abs(p - sm) >= abs(error):
                    pass
                    #print("warning: error =", error, "for t = ", t, "new error", p - sm)
                error = p - sm
            print(s, x, round(p, round_to), round(sm, round_to), round( p-sm, round_to*round_to), sep="\t")

main()