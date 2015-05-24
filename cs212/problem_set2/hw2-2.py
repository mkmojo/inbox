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
