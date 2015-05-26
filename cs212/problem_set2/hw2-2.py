#------------------
# User Instructions
#
# Hopper, Kay, Liskov, Perlis, and Ritchie live on
# different floors of a five-floor apartment building.
#
# Hopper does not live on the top floor.
# Kay does not live on the bottom floor.
# Liskov does not live on either the top or the bottom floor.
# Perlis lives on a higher floor than does Kay.
# Ritchie does not live on a floor adjacent to Liskov's.
# Liskov does not live on a floor adjacent to Kay's.
#
# Where does everyone live?
#
# Write a function floor_puzzle() that returns a list of
# five floor numbers denoting the floor of Hopper, Kay,
# Liskov, Perlis, and Ritchie.

import itertools
import time

#profiling purpose
def c(sequence):
    c.starts += 1
    for item in sequence:
        c.items += 1
        yield item

#c is a generator function
#generator function often goes with next() / for statement
def instrument_fn(fn, *args):
    c.starts, c.items = 0, 0
    result = fn(*args)
    print '%s got %s with %5d iters over %7d items' % (
            fn.__name__, result, c.starts, c.items)

def floor_puzzle():
    bottom = 1
    top = 5
    for Hopper, Kay, Liskov, Perlis, Ritchie in itertools.permutations(range(1,6)):
        if Kay != bottom and \
        Liskov != top and Liskov != bottom and \
        Perlis > Kay and \
        Ritchie != Liskov + 1 and Ritchie != Liskov -1 and \
        Liskov != Kay + 1 and Liskov != Kay -1:
            return [Hopper, Kay, Liskov, Perlis, Ritchie]

# *args, packing and unpacking arguments
def timecall(fn, *args): # packing
    t0 = time.clock()
    result = fn(*args)  # unpacking
    t1 = time.clock()
    return t1-t0, result

def timecalls(n, fn, *args):
    """call function n times and return min, avg, and max time"""
    if isinstance(n, int):
        times = [timecall(fn, *args)[0] for _ in range(n)]
    elif isinstance(n, float):
        times = []
        while sum(times) < n:
            times.append(timecall(fn, *args)[0])
    return min(times), average(times), max(times)

def average(numbers):
    return sum(numbers)/float(len(numbers))

print timecalls(50, floor_puzzle,)
