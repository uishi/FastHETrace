import math

def next_power_of_2(n):

    # decrement n (to handle cases when n itself
    # is a power of 2)
    n = n - 1

    # do till only one bit is left
    while n & n - 1:
        n = n & n - 1  # unset rightmost bit

    # n is now a power of two (less than n)

    # return next power of 2
    return n << 1

def is_pow_of_2(n):
    return (n & (n-1) == 0) and n != 0

# Return n1 and n2
def bsgs_pow_of_2(n):
    if not is_pow_of_2(n):
        raise Exception("Must be power of 2")

    if n <= 4:
        raise Exception("No need to use BSGS for n <= 4")
    n1 = math.ceil(math.sqrt(n)) # real
    # a divisor of n: O(\sqrt{n})
    n1 = next_power_of_2(n1) # int rounded
    n2 = int(n / n1)
    if not is_pow_of_2(n1):
        raise Exception("n1 must be power of 2")

    if not is_pow_of_2(n2):
        raise Exception("n2 must be power of 2")

    return n1, n2
